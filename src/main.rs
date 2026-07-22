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

/// Parse a single allele cell into an integer, matching the Python semantics:
/// strip one leading `INF-`, trim whitespace, missing char -> 0, and anything
/// that is not a plain non-negative integer -> 0.
#[inline]
fn parse_allele(cell: &str, skip_repl: bool, missing: &str) -> i32 {
    let mut s = cell.trim();
    if !skip_repl {
        if let Some(rest) = s.strip_prefix("INF-") {
            s = rest;
        }
    }
    if s == missing || s.is_empty() {
        return 0;
    }
    if s.bytes().all(|b| b.is_ascii_digit()) {
        s.parse::<i32>().unwrap_or(0)
    } else {
        0
    }
}

/// Hamming distance: count loci where alleles differ and both are present
/// (non-zero). Optional early-exit at `max_dist`.
#[inline]
fn distance(a: &[i32], b: &[i32], max_dist: i32) -> i32 {
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

fn main() {
    let args = Args::parse();
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

    let parsed: Vec<(&str, Vec<i32>)> = data_lines
        .par_iter()
        .map(|line| {
            let mut it = line.split(sep_ch);
            let id = it.next().unwrap_or("");
            let mut calls = Vec::with_capacity(n_loci);
            for cell in it {
                calls.push(parse_allele(cell, skip_repl, missing));
            }
            calls.resize(n_loci, 0);
            (id, calls)
        })
        .collect();

    let mut ids: Vec<&str> = parsed.iter().map(|(id, _)| *id).collect();
    let mut calls = vec![0i32; n * n_loci];
    for (row, (_, c)) in parsed.iter().enumerate() {
        calls[row * n_loci..(row + 1) * n_loci].copy_from_slice(c);
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
            let mut nc = vec![0i32; n * keep.len()];
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
            let mut nc = vec![0i32; keep.len() * n_loci];
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
    // upper[i] holds distances to j = i+1 .. n-1 (length n-1-i).
    let upper: Vec<Vec<i32>> = (0..n)
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
        .collect();
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
        assert_eq!(parse_allele("5", false, "-"), 5);
        assert_eq!(parse_allele("INF-42", false, "-"), 42); // strip one INF-
        assert_eq!(parse_allele("INF-INF-42", false, "-"), 0); // double -> not numeric -> 0
        assert_eq!(parse_allele("-", false, "-"), 0); // missing
        assert_eq!(parse_allele(" 7 ", false, "-"), 7); // trimmed
        assert_eq!(parse_allele("", false, "-"), 0);
        assert_eq!(parse_allele("abc", false, "-"), 0);
        assert_eq!(parse_allele("12x", false, "-"), 0);
        // with --skip-input-replacements, INF- is left in place -> not numeric -> 0
        assert_eq!(parse_allele("INF-9", true, "-"), 0);
        // custom missing char
        assert_eq!(parse_allele("NA", false, "NA"), 0);
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
