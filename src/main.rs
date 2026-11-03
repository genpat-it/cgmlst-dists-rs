//! cgmlst-dists-rs — ultra-fast pairwise Hamming distance matrix for cgMLST
//! allele call tables. A Rust port of cgmlst-dists-py: same distance semantics
//! (differences at loci where both alleles are present/non-zero), same chewBBACA
//! input handling (strip a leading `INF-`, treat the missing char as 0), and the
//! same output layouts (full / lower-tri / upper-tri, to a file or stdout).

use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::time::Instant;

use clap::Parser;
use rayon::prelude::*;

mod gpu;

#[derive(Parser, Debug)]
#[command(name = "cgmlst-dists", version, about = "Pairwise cg-MLST Hamming distance matrix (Rust)")]
struct Args {
    /// Input allele-call TSV (chewBBACA format). Use '-' for stdin.
    #[arg(short = 'i', long)]
    input: String,

    /// Output file [default: stdout]
    #[arg(short = 'o', long)]
    output: Option<String>,

    /// Number of threads [default: all cores]
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,

    /// Output matrix layout
    #[arg(short = 'm', long, value_enum, default_value_t = MatrixFormat::Full)]
    matrix_format: MatrixFormat,

    /// Input field separator
    #[arg(short = 'd', long, default_value = "\t")]
    input_sep: String,

    /// Output field separator
    #[arg(short = 'D', long, default_value = "\t")]
    output_sep: String,

    /// Name for the top-left index cell
    #[arg(short = 'x', long, default_value = "cgmlst-dists")]
    index_name: String,

    /// Character used for missing data
    #[arg(short = 'n', long, default_value = "-")]
    missing_char: String,

    /// Do not strip the chewBBACA `INF-` prefix
    #[arg(short = 'r', long, default_value_t = false)]
    skip_input_replacements: bool,

    /// Count chewBBACA `*N` novel alleles (inferred, not added to the schema) as
    /// alleles, compared by their text within each locus. Default: missing, as
    /// in tseemann cgmlst-dists. The `*N` numbering is local to one chewBBACA
    /// run, so use this only on tables from a single run.
    #[arg(short = 'N', long, default_value_t = false)]
    novel_as_allele: bool,

    /// Keep only loci with at least this % of non-missing calls (0-100)
    #[arg(short = 'L', long)]
    locus_completeness: Option<f64>,

    /// Keep only samples with at least this % of non-missing calls (0-100)
    #[arg(short = 'S', long)]
    sample_completeness: Option<f64>,

    /// gzip-compress the output (also auto-enabled when -o ends with .gz).
    /// Default is uncompressed plain text.
    #[arg(short = 'z', long, default_value_t = false)]
    gzip: bool,

    /// Stop counting beyond this distance (early-exit; 0 = no cap)
    #[arg(short = 'X', long, default_value_t = 0)]
    max_dist: i32,

    /// Skip the up-front memory feasibility check
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,

    /// Compute on the GPU (build with `--features cuda`). Same output as the CPU path.
    #[arg(short = 'g', long, default_value_t = false)]
    gpu: bool,

    /// CUDA device index for --gpu
    #[arg(long, default_value_t = 0)]
    gpu_device: usize,

    /// Rows per GPU stripe for --gpu [default: 0 = sized from free VRAM]
    #[arg(long, default_value_t = 0, hide = true)]
    gpu_stripe: usize,

    /// Suppress progress/log messages on stderr
    #[arg(short = 's', long, default_value_t = false)]
    silent: bool,
}

#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq)]
enum MatrixFormat {
    Full,
    #[value(name = "lower-tri")]
    LowerTri,
    #[value(name = "upper-tri")]
    UpperTri,
}

/// Parse a single allele cell into an integer: strip one leading `INF-`, trim
/// whitespace, missing char -> 0, and anything that is not a plain non-negative
/// integer -> 0. Alleles are read as u64 so chewBBACA hashed profiles (CRC32,
/// up to 2^32-1) and any other large integer id are kept, not lost as missing.
/// Returns `Err` for an integer that does not fit in u64.
#[inline]
fn parse_allele(cell: &str, skip_repl: bool, missing: &str) -> Result<u64, ()> {
    let mut s = cell.trim();
    if !skip_repl {
        if let Some(rest) = s.strip_prefix("INF-") {
            s = rest;
        }
    }
    if s == missing || s.is_empty() {
        return Ok(0);
    }
    if s.bytes().all(|b| b.is_ascii_digit()) {
        s.parse::<u64>().map_err(|_| ())
    } else {
        Ok(0)
    }
}

/// A `*N` novel allele: the number after the `*`, or None if `cell` is not one.
#[inline]
fn parse_novel(cell: &str) -> Option<u64> {
    let n = cell.trim().strip_prefix('*')?;
    if !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) {
        n.parse().ok()
    } else {
        None
    }
}

/// An allele that cannot be stored directly as a u32 id.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Special {
    /// integer id above u32::MAX
    Big(u64),
    /// `*N` novel allele (with --novel-as-allele)
    Novel(u64),
}

/// Re-encode loci holding special alleles: each such locus gets dense ids 1..k
/// (0 stays missing), which preserves equality within the locus — the only
/// thing the distance depends on. `special` lists (row, locus, allele) for every
/// cell that could not be stored as is (its slot in `calls` is a placeholder).
fn remap_special_loci(calls: &mut [u32], n_loci: usize, special: &[(usize, usize, Special)]) {
    use std::collections::HashMap;
    let mut by_locus: HashMap<usize, Vec<(usize, Special)>> = HashMap::new();
    for &(row, k, v) in special {
        by_locus.entry(k).or_default().push((row, v));
    }
    let n = calls.len() / n_loci.max(1);
    // An integer id that fit in u32 and a special allele are different keys.
    #[derive(PartialEq, Eq, Hash)]
    enum Key {
        Id(u32),
        Special(Special),
    }
    for (k, cells) in by_locus {
        let mut col: Vec<Key> = (0..n).map(|i| Key::Id(calls[i * n_loci + k])).collect();
        for (row, v) in cells {
            col[row] = Key::Special(v);
        }
        let mut ids: HashMap<Key, u32> = HashMap::new();
        for (i, v) in col.into_iter().enumerate() {
            calls[i * n_loci + k] = if v == Key::Id(0) {
                0
            } else {
                let next = ids.len() as u32 + 1;
                *ids.entry(v).or_insert(next)
            };
        }
    }
}

/// Hamming distance: count loci where alleles differ and both are present
/// (non-zero). Optional early-exit at `max_dist`.
#[inline]
fn distance(a: &[u32], b: &[u32], max_dist: i32) -> i32 {
    let mut d: i32 = 0;
    for k in 0..a.len() {
        let x = unsafe { *a.get_unchecked(k) };
        let y = unsafe { *b.get_unchecked(k) };
        if x != y && x != 0 && y != 0 {
            d += 1;
            if max_dist > 0 && d >= max_dist {
                return max_dist;
            }
        }
    }
    d
}

/// Symmetric lookup into the upper-triangle storage: distance between i and j.
#[inline]
fn at(upper: &[Vec<i32>], i: usize, j: usize) -> i32 {
    if i == j {
        0
    } else if j > i {
        upper[i][j - i - 1]
    } else {
        upper[j][i - j - 1]
    }
}

fn read_input(path: &str) -> io::Result<String> {
    let mut buf = String::new();
    if path == "-" {
        io::stdin().read_to_string(&mut buf)?;
    } else {
        File::open(path)?.read_to_string(&mut buf)?;
    }
    Ok(buf)
}

/// Linux: available memory in bytes (MemAvailable from /proc/meminfo).
fn available_memory_bytes() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in s.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            let kb: u64 = rest.trim().trim_end_matches("kB").trim().parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

/// Write the header + matrix rows to `out`. Rows are formatted into byte
/// buffers in parallel (chunked to bound memory) and written in order, so the
/// per-row integer formatting (the CPU bottleneck at scale) is spread across
/// cores. Works with any writer, including a parallel gzip stream.
fn write_matrix<W: Write>(
    out: &mut W,
    ids: &[&str],
    upper: &[Vec<i32>],
    n: usize,
    osep: &[u8],
    index_name: &str,
    mfmt: MatrixFormat,
) {
    out.write_all(index_name.as_bytes()).unwrap();
    for id in ids {
        out.write_all(osep).unwrap();
        out.write_all(id.as_bytes()).unwrap();
    }
    out.write_all(b"\n").unwrap();

    const CHUNK: usize = 1024;
    let mut start = 0usize;
    while start < n {
        let end = (start + CHUNK).min(n);
        let bufs: Vec<Vec<u8>> = (start..end)
            .into_par_iter()
            .map(|i| {
                let mut buf: Vec<u8> = Vec::with_capacity((n + 1) * 7);
                let mut ib = itoa::Buffer::new();
                buf.extend_from_slice(ids[i].as_bytes());
                for j in 0..n {
                    buf.extend_from_slice(osep);
                    let v = match mfmt {
                        MatrixFormat::Full => at(upper, i, j),
                        MatrixFormat::LowerTri => if j <= i { at(upper, i, j) } else { 0 },
                        MatrixFormat::UpperTri => if j >= i { at(upper, i, j) } else { 0 },
                    };
                    buf.extend_from_slice(ib.format(v).as_bytes());
                }
                buf.push(b'\n');
                buf
            })
            .collect();
        for b in &bufs {
            out.write_all(b).unwrap();
        }
        start = end;
    }
}

#[cfg(feature = "cuda")]
fn gpu_upper(calls: &[u32], n: usize, n_loci: usize, max_dist: i32, args: &Args) -> Vec<Vec<i32>> {
    if gpu::device_count() <= args.gpu_device {
        eprintln!("error: --gpu: CUDA device {} not available", args.gpu_device);
        std::process::exit(1);
    }
    match gpu::upper_triangle(calls, n, n_loci, max_dist, args.gpu_device, args.gpu_stripe) {
        Ok((upper, secs)) => {
            if !args.silent {
                eprintln!("GPU pipeline: {:.2}s", secs);
            }
            upper
        }
        Err(e) => {
            eprintln!("error: --gpu: {e}");
            std::process::exit(1);
        }
    }
}

#[cfg(not(feature = "cuda"))]
fn gpu_upper(_: &[u32], _: usize, _: usize, _: i32, _: &Args) -> Vec<Vec<i32>> {
    eprintln!("error: --gpu needs a build with `cargo build --release --features cuda`");
    std::process::exit(1);
}

fn main() {
    let args = Args::parse();
    if args.gpu && !cfg!(feature = "cuda") {
        eprintln!("error: --gpu needs a build with `cargo build --release --features cuda`");
        std::process::exit(1);
    }
    let t_start = Instant::now();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let sep = args.input_sep.as_bytes();
    if sep.len() != 1 {
        eprintln!("error: input separator must be a single character");
        std::process::exit(1);
    }
    let sep_ch = sep[0] as char;

    // ---- Load ----
    let t_load = Instant::now();
    let text = match read_input(&args.input) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("error reading {}: {}", args.input, e);
            std::process::exit(1);
        }
    };

    let mut lines = text.lines();
    let header = match lines.next() {
        Some(h) => h,
        None => {
            eprintln!("error: empty input");
            std::process::exit(1);
        }
    };
    let mut n_loci = header.split(sep_ch).count() - 1;

    let data_lines: Vec<&str> = lines.filter(|l| !l.is_empty()).collect();
    let mut n = data_lines.len();

    let skip_repl = args.skip_input_replacements;
    let missing = args.missing_char.as_str();
    let novel_as_allele = args.novel_as_allele;

    // Values that fit in u32 (all of them in practice, CRC32 hashes included)
    // go straight into the u32 matrix; larger ones (and `*N` novel alleles with
    // --novel-as-allele) are collected and their loci re-encoded afterwards, so
    // the common case pays nothing extra.
    let parsed: Vec<(&str, Vec<u32>, Vec<(usize, Special)>)> = data_lines
        .par_iter()
        .map(|line| {
            let mut it = line.split(sep_ch);
            let id = it.next().unwrap_or("");
            let mut calls = Vec::with_capacity(n_loci);
            let mut big = Vec::new();
            for (k, cell) in it.enumerate() {
                let v = parse_allele(cell, skip_repl, missing).unwrap_or_else(|_| {
                    eprintln!("error: sample {}: allele '{}' is an integer too large (> {})", id, cell.trim(), u64::MAX);
                    std::process::exit(1);
                });
                let special = if v > u32::MAX as u64 {
                    Some(Special::Big(v))
                } else if v == 0 && novel_as_allele {
                    parse_novel(cell).map(Special::Novel)
                } else {
                    None
                };
                match special {
                    Some(sp) => {
                        big.push((k, sp));
                        calls.push(u32::MAX);
                    }
                    None => calls.push(v as u32),
                }
            }
            calls.resize(n_loci, 0);
            (id, calls, big)
        })
        .collect();

    let mut ids: Vec<&str> = parsed.iter().map(|(id, _, _)| *id).collect();
    let mut calls = vec![0u32; n * n_loci];
    let mut big: Vec<(usize, usize, Special)> = Vec::new();
    for (row, (_, c, b)) in parsed.iter().enumerate() {
        calls[row * n_loci..(row + 1) * n_loci].copy_from_slice(&c[..n_loci]);
        big.extend(b.iter().filter(|&&(k, _)| k < n_loci).map(|&(k, v)| (row, k, v)));
    }
    if !big.is_empty() {
        remap_special_loci(&mut calls, n_loci, &big);
    }
    drop(parsed);
    if !args.silent {
        eprintln!("Loaded {} samples x {} loci in {:.2}s", n, n_loci, t_load.elapsed().as_secs_f64());
    }

    // ---- Completeness filtering (present = non-zero call) ----
    if let Some(thr) = args.locus_completeness {
        let keep: Vec<usize> = (0..n_loci)
            .into_par_iter()
            .filter(|&k| {
                let present = (0..n).filter(|&i| calls[i * n_loci + k] != 0).count();
                (present as f64 / n as f64) * 100.0 >= thr
            })
            .collect();
        if keep.len() != n_loci {
            let mut nc = vec![0u32; n * keep.len()];
            nc.par_chunks_mut(keep.len()).enumerate().for_each(|(i, row)| {
                for (dst, &k) in keep.iter().enumerate() {
                    row[dst] = calls[i * n_loci + k];
                }
            });
            calls = nc;
            n_loci = keep.len();
        }
        if !args.silent {
            eprintln!("Locus completeness >= {}%: kept {} loci", thr, n_loci);
        }
    }

    if let Some(thr) = args.sample_completeness {
        let keep: Vec<usize> = (0..n)
            .into_par_iter()
            .filter(|&i| {
                let present = calls[i * n_loci..(i + 1) * n_loci].iter().filter(|&&v| v != 0).count();
                (present as f64 / n_loci as f64) * 100.0 >= thr
            })
            .collect();
        if keep.len() != n {
            let mut nc = vec![0u32; keep.len() * n_loci];
            let mut nids = Vec::with_capacity(keep.len());
            for (dst, &i) in keep.iter().enumerate() {
                nc[dst * n_loci..(dst + 1) * n_loci].copy_from_slice(&calls[i * n_loci..(i + 1) * n_loci]);
                nids.push(ids[i]);
            }
            calls = nc;
            ids = nids;
            n = keep.len();
        }
        if !args.silent {
            eprintln!("Sample completeness >= {}%: kept {} samples", thr, n);
        }
    }

    // ---- Memory feasibility guard (upper-triangle storage) ----
    let tri = (n as u64) * (n.saturating_sub(1) as u64) / 2;
    let need_bytes = tri.saturating_mul(4);
    if !args.force {
        if let Some(avail) = available_memory_bytes() {
            if need_bytes as f64 > 0.85 * avail as f64 {
                eprintln!(
                    "\nERROR: not enough memory.\n  {} samples -> upper-triangle matrix ~{:.1} GiB, available ~{:.1} GiB.\n  Aborting before computation. Use --force to try anyway, or reduce samples.",
                    n,
                    need_bytes as f64 / 1024f64.powi(3),
                    avail as f64 / 1024f64.powi(3)
                );
                std::process::exit(1);
            }
        }
    }

    // ---- Compute upper triangle (jagged, parallel over rows) ----
    let t_calc = Instant::now();
    let max_dist = args.max_dist;
    let upper: Vec<Vec<i32>> = if args.gpu {
        gpu_upper(&calls, n, n_loci, max_dist, &args)
    } else {
        // upper[i] holds distances to j = i+1 .. n-1 (length n-1-i).
        (0..n)
            .into_par_iter()
            .map(|i| {
                let ai = &calls[i * n_loci..(i + 1) * n_loci];
                let mut row = Vec::with_capacity(n - 1 - i.min(n - 1));
                for j in (i + 1)..n {
                    let aj = &calls[j * n_loci..(j + 1) * n_loci];
                    row.push(distance(ai, aj, max_dist));
                }
                row
            })
            .collect()
    };
    if !args.silent {
        eprintln!("Computed {} x {} matrix in {:.2}s", n, n, t_calc.elapsed().as_secs_f64());
    }

    // ---- Write ----
    let t_save = Instant::now();
    let osep: Vec<u8> = args.output_sep.as_bytes().to_vec();
    let gzip = args.gzip || args.output.as_deref().map_or(false, |p| p.ends_with(".gz"));
    let file_writer: Box<dyn Write + Send> = match &args.output {
        Some(p) => Box::new(File::create(p).unwrap_or_else(|e| {
            eprintln!("error creating {}: {}", p, e);
            std::process::exit(1);
        })),
        None => Box::new(io::stdout()),
    };
    // Optional gzip (when -z/--gzip or the output name ends with .gz): the matrix
    // is highly redundant and compresses ~4x, so far fewer bytes cross a slow
    // Docker bind mount. gzp compresses blocks across all cores (standard .gz,
    // readable by gunzip/zcat). Default output is uncompressed plain text.
    let mfmt = args.matrix_format;
    if gzip {
        use gzp::{deflate::Gzip, par::compress::ParCompressBuilder, ZWriter};
        let mut z = ParCompressBuilder::<Gzip>::new()
            .compression_level(gzp::Compression::new(3))
            .from_writer(file_writer);
        write_matrix(&mut z, &ids, &upper, n, &osep, &args.index_name, mfmt);
        z.finish().unwrap();
    } else {
        let mut out = BufWriter::with_capacity(16 * 1024 * 1024, file_writer);
        write_matrix(&mut out, &ids, &upper, n, &osep, &args.index_name, mfmt);
        out.flush().unwrap();
    }

    if !args.silent {
        eprintln!(
            "Wrote output in {:.2}s. Total {:.2}s",
            t_save.elapsed().as_secs_f64(),
            t_start.elapsed().as_secs_f64()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_semantics() {
        let p = |c: &str| parse_allele(c, false, "-").unwrap();
        assert_eq!(p("5"), 5);
        assert_eq!(p("INF-42"), 42); // strip one INF-
        assert_eq!(p("INF-INF-42"), 0); // double -> not numeric -> 0
        assert_eq!(p("-"), 0); // missing
        assert_eq!(p(" 7 "), 7); // trimmed
        assert_eq!(p(""), 0);
        assert_eq!(p("abc"), 0);
        assert_eq!(p("12x"), 0);
        assert_eq!(p("LNF"), 0); // chewBBACA classes -> missing
        assert_eq!(p("*486"), 0); // novel (non-numeric) -> missing, as tseemann
        // with --skip-input-replacements, INF- is left in place -> not numeric -> 0
        assert_eq!(parse_allele("INF-9", true, "-").unwrap(), 0);
        // custom missing char
        assert_eq!(parse_allele("NA", false, "NA").unwrap(), 0);
    }

    #[test]
    fn parse_large_integers() {
        let p = |c: &str| parse_allele(c, false, "-").unwrap();
        // CRC32 hashes above i32::MAX were silently read as missing (0)
        assert_eq!(p("2147483648"), 2_147_483_648);
        assert_eq!(p("3735928559"), 3_735_928_559);
        assert_eq!(p("4294967295"), u32::MAX as u64);
        assert_eq!(p("INF-4294967295"), u32::MAX as u64);
        assert_eq!(p("4294967296"), 4_294_967_296);
        assert_eq!(p("18446744073709551615"), u64::MAX);
        assert!(parse_allele("18446744073709551616", false, "-").is_err());
    }

    #[test]
    fn hashed_alleles_distance() {
        let h = |c: &str| parse_allele(c, false, "-").unwrap() as u32;
        let a = [h("4294967295"), h("2147483648"), h("3000000000"), h("-"), h("123")];
        let b = [h("4294967295"), h("2147483649"), h("3000000000"), h("4000000000"), h("124")];
        // equal near 2^32-1, differ by one just above 2^31, equal, missing on a, differ
        assert_eq!(distance(&a, &b, 0), 2);
        assert_eq!(distance(&a, &a, 0), 0);
    }

    #[test]
    fn remap_values_above_u32() {
        // 2 rows x 3 loci; locus 1 holds values > u32::MAX
        let big = [(0usize, 1usize, Special::Big(5_000_000_000)), (1, 1, Special::Big(5_000_000_000))];
        let mut calls = vec![7, u32::MAX, 0, 7, u32::MAX, 9];
        remap_special_loci(&mut calls, 3, &big);
        assert_eq!(distance(&calls[0..3], &calls[3..6], 0), 0); // equal big values
        let big = [(0usize, 1usize, Special::Big(5_000_000_000)), (1, 1, Special::Big(5_000_000_001))];
        let mut calls = vec![7, u32::MAX, 0, 7, u32::MAX, 9];
        remap_special_loci(&mut calls, 3, &big);
        assert_eq!(distance(&calls[0..3], &calls[3..6], 0), 1); // different big values
        // a big value and a u32 value in the same locus must stay distinct
        let big = [(0usize, 0usize, Special::Big((u32::MAX as u64) + 1))];
        let mut calls = vec![u32::MAX, 1, u32::MAX, 1];
        remap_special_loci(&mut calls, 2, &big);
        assert_eq!(distance(&calls[0..2], &calls[2..4], 0), 1);
        let mut calls = vec![0, 1, 0, 1]; // missing stays missing
        remap_special_loci(&mut calls, 2, &[(1, 0, Special::Big(6_000_000_000))]);
        assert_eq!(calls[0], 0);
        assert_eq!(distance(&calls[0..2], &calls[2..4], 0), 0);
    }

    #[test]
    fn novel_alleles() {
        assert_eq!(parse_novel("*486"), Some(486));
        assert_eq!(parse_novel(" *7 "), Some(7));
        for c in ["486", "*", "*x", "*12a", "LNF", "INF-3", ""] {
            assert_eq!(parse_novel(c), None, "{c}");
        }
        // rows: [*5, 7, *5], [*5, 7, *6]; locus 1 is a plain id
        let sp = [
            (0usize, 0usize, Special::Novel(5)),
            (0, 2, Special::Novel(5)),
            (1, 0, Special::Novel(5)),
            (1, 2, Special::Novel(6)),
        ];
        let mut calls = vec![u32::MAX, 7, u32::MAX, u32::MAX, 7, u32::MAX];
        remap_special_loci(&mut calls, 3, &sp);
        assert_eq!(distance(&calls[0..3], &calls[3..6], 0), 1); // same *5, different *5/*6
        // a novel *5 and an integer allele 5 at the same locus are different
        let mut calls = vec![u32::MAX, 5];
        remap_special_loci(&mut calls, 1, &[(0, 0, Special::Novel(5))]);
        assert_eq!(distance(&calls[0..1], &calls[1..2], 0), 1);
    }

    #[test]
    fn hamming_semantics() {
        assert_eq!(distance(&[1, 2, 3], &[1, 9, 3], 0), 1); // one difference
        assert_eq!(distance(&[1, 2, 3], &[1, 2, 3], 0), 0); // identical
        assert_eq!(distance(&[1, 2, 3], &[4, 5, 6], 0), 3); // all differ
        assert_eq!(distance(&[0, 2, 3], &[1, 2, 3], 0), 0); // missing on a -> ignored
        assert_eq!(distance(&[1, 2, 3], &[1, 2, 0], 0), 0); // missing on b -> ignored
        assert_eq!(distance(&[0, 0, 0], &[1, 2, 3], 0), 0); // all missing on one side
    }

    #[test]
    fn max_dist_cap() {
        // early-exit stops counting at the cap
        assert_eq!(distance(&[1, 2, 3, 4], &[5, 6, 7, 8], 2), 2);
        assert_eq!(distance(&[1, 2, 3, 4], &[5, 6, 7, 8], 0), 4); // 0 = no cap
    }
}
