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

    /// Stop counting beyond this distance (early-exit; 0 = no cap)
    #[arg(short = 'X', long, default_value_t = 0)]
    max_dist: i32,

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
    // Only accept pure ASCII digits (mirrors `^[0-9]+$`); otherwise 0.
    if s.bytes().all(|b| b.is_ascii_digit()) {
        s.parse::<i32>().unwrap_or(0)
    } else {
        0
    }
}

/// Hamming distance between two allele profiles: count loci where the alleles
/// differ and both are present (non-zero). Optional early-exit at `max_dist`.
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

fn read_input(path: &str) -> io::Result<String> {
    let mut buf = String::new();
    if path == "-" {
        io::stdin().read_to_string(&mut buf)?;
    } else {
        File::open(path)?.read_to_string(&mut buf)?;
    }
    Ok(buf)
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
    let n_loci = header.split(sep_ch).count() - 1;

    // Collect data lines, then parse rows in parallel.
    let data_lines: Vec<&str> = lines.filter(|l| !l.is_empty()).collect();
    let n = data_lines.len();

    let skip_repl = args.skip_input_replacements;
    let missing = args.missing_char.as_str();

    let parsed: Vec<(usize, &str, Vec<i32>)> = data_lines
        .par_iter()
        .enumerate()
        .map(|(row, line)| {
            let mut it = line.split(sep_ch);
            let id = it.next().unwrap_or("");
            let mut calls = Vec::with_capacity(n_loci);
            for cell in it {
                calls.push(parse_allele(cell, skip_repl, missing));
            }
            // Pad/truncate defensively to n_loci
            calls.resize(n_loci, 0);
            (row, id, calls)
        })
        .collect();

    let ids: Vec<&str> = parsed.iter().map(|(_, id, _)| *id).collect();
    // Flat row-major matrix of calls for cache-friendly access.
    let mut calls = vec![0i32; n * n_loci];
    for (row, _, c) in &parsed {
        calls[row * n_loci..(row + 1) * n_loci].copy_from_slice(c);
    }
    drop(parsed);
    let load_secs = t_load.elapsed().as_secs_f64();
    if !args.silent {
        eprintln!("Loaded {} samples x {} loci in {:.2}s", n, n_loci, load_secs);
    }

    // ---- Compute (upper triangle, parallel over rows) ----
    let t_calc = Instant::now();
    let max_dist = args.max_dist;
    let mut dist = vec![0i32; n * n];
    // Fill upper triangle: dist[i*n + j] for j > i. Each row is written by one
    // task; distances are symmetric so the lower triangle is read from the
    // upper one at output time (no separate mirror pass, no aliasing).
    dist.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
        let ai = &calls[i * n_loci..(i + 1) * n_loci];
        for j in (i + 1)..n {
            let aj = &calls[j * n_loci..(j + 1) * n_loci];
            row[j] = distance(ai, aj, max_dist);
        }
    });
    let calc_secs = t_calc.elapsed().as_secs_f64();
    if !args.silent {
        eprintln!("Computed {} x {} matrix in {:.2}s", n, n, calc_secs);
    }

    // symmetric lookup: value at (i, j)
    let at = |i: usize, j: usize| -> i32 {
        if i == j {
            0
        } else if j > i {
            dist[i * n + j]
        } else {
            dist[j * n + i]
        }
    };

    // ---- Write ----
    let t_save = Instant::now();
    let osep = args.output_sep.as_bytes();
    let writer: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(File::create(p).unwrap_or_else(|e| {
            eprintln!("error creating {}: {}", p, e);
            std::process::exit(1);
        })),
        None => Box::new(io::stdout()),
    };
    let mut out = BufWriter::with_capacity(16 * 1024 * 1024, writer);

    // Header
    out.write_all(args.index_name.as_bytes()).unwrap();
    for id in &ids {
        out.write_all(osep).unwrap();
        out.write_all(id.as_bytes()).unwrap();
    }
    out.write_all(b"\n").unwrap();

    let mut ibuf = itoa::Buffer::new();
    for i in 0..n {
        out.write_all(ids[i].as_bytes()).unwrap();
        for j in 0..n {
            out.write_all(osep).unwrap();
            let v = match args.matrix_format {
                MatrixFormat::Full => at(i, j),
                MatrixFormat::LowerTri => {
                    if j <= i {
                        at(i, j)
                    } else {
                        0
                    }
                }
                MatrixFormat::UpperTri => {
                    if j >= i {
                        at(i, j)
                    } else {
                        0
                    }
                }
            };
            out.write_all(ibuf.format(v).as_bytes()).unwrap();
        }
        out.write_all(b"\n").unwrap();
    }
    out.flush().unwrap();
    let save_secs = t_save.elapsed().as_secs_f64();

    if !args.silent {
        eprintln!(
            "Wrote output in {:.2}s. Total {:.2}s (load {:.2} / calc {:.2} / save {:.2})",
            save_secs,
            t_start.elapsed().as_secs_f64(),
            load_secs,
            calc_secs,
            save_secs
        );
    }
}
