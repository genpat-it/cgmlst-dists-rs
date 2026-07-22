# cgmlst-dists-rs

[![Crates.io](https://img.shields.io/crates/v/cgmlst-dists-rs.svg?logo=rust)](https://crates.io/crates/cgmlst-dists-rs)
[![CI](https://github.com/genpat-it/cgmlst-dists-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/genpat-it/cgmlst-dists-rs/actions/workflows/ci.yml)
[![Docker](https://img.shields.io/badge/docker-ghcr.io-blue?logo=docker)](https://github.com/genpat-it/cgmlst-dists-rs/pkgs/container/cgmlst-dists-rs)
[![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)

Ultra-fast pairwise **cgMLST Hamming distance matrix** from chewBBACA allele
call tables — a Rust reimplementation of
[cgmlst-dists-py](https://github.com/genpat-it/cgmlst-dists-py), producing
**byte-identical output** while running several times faster.

The distance kernel **auto-vectorizes to AVX-512 SIMD** (16 loci compared per
instruction), which is the main reason it outperforms the original C
`cgmlst-dists` and the Python implementation.

## Why it's fast

- **SIMD Hamming kernel**: the inner loop (`differ && both present`) is written
  so LLVM vectorizes it to packed compare + mask + packed-add (AVX-512 on
  capable CPUs). Verified in the disassembly (`vpcmpeqd`/`vptestmd`/`vpaddd`).
- **rayon** work-stealing parallelism over the (imbalanced) triangular workload.
- **Upper-triangle storage** with symmetric output lookup — ~half the memory of
  a full `N×N` matrix (e.g. ~6.7 GB for 50k samples instead of ~10 GB).
- Fast buffered integer output (`itoa`).

## Install

### Docker (no Rust toolchain needed)

```bash
docker pull ghcr.io/genpat-it/cgmlst-dists-rs
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd):/data" \
  ghcr.io/genpat-it/cgmlst-dists-rs -i /data/chewbbaca.tab -o /data/distances.tsv
```

The image (Debian slim) is built with portable AVX2 SIMD (`x86-64-v3`), so it
runs on any x86-64 Linux host from ~2013 onward.

### Cargo (crates.io)

```bash
cargo install cgmlst-dists-rs
```

### From source

```bash
git clone https://github.com/genpat-it/cgmlst-dists-rs.git
cd cgmlst-dists-rs
RUSTFLAGS="-C target-cpu=native" cargo build --release
# binary: target/release/cgmlst-dists
```

`-C target-cpu=native` enables the widest SIMD your CPU supports (e.g. AVX-512)
and is recommended for best performance when building for your own machine.

## Usage

```bash
cgmlst-dists -i chewbbaca.tab -o distances.tsv
cgmlst-dists -i chewbbaca.tab > distances.tsv        # stdout
cgmlst-dists -i in.tab -m lower-tri -o dist.tsv       # triangular
```

| flag | meaning |
|------|---------|
| `-i, --input` | input allele-call TSV (`-` for stdin) |
| `-o, --output` | output file [default: stdout] |
| `-t, --threads` | threads [default: all cores] |
| `-m, --matrix-format` | `full` \| `lower-tri` \| `upper-tri` |
| `-d, --input-sep` / `-D, --output-sep` | field separators [tab] |
| `-x, --index-name` | top-left cell name [`cgmlst-dists`] |
| `-n, --missing-char` | missing-data character [`-`] |
| `-r, --skip-input-replacements` | do not strip the chewBBACA `INF-` prefix |
| `-L, --locus-completeness` / `-S, --sample-completeness` | keep loci/samples with ≥ this % non-missing calls |
| `-X, --max-dist` | stop counting beyond this distance (early-exit; 0 = no cap) |
| `-f, --force` | skip the up-front memory feasibility check |
| `-s, --silent` | suppress progress on stderr |

### Large outputs / slow Docker bind mounts

At large scale the output matrix is huge (~12 GB for 50k samples). Row
formatting is parallelized, so on a local disk writing is fast (~5 s). If you
run in Docker and write through a **bind mount** (which can be slow, especially
on Docker Desktop), gzip the output to move ~4× fewer bytes:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd):/data" \
  ghcr.io/genpat-it/cgmlst-dists-rs -i /data/in.tab -o /data/distances.tsv.gz
# read it back with: zcat distances.tsv.gz
```

Gzip is opt-in (a `.gz` output name, or `-z/--gzip`); the default is plain text.

**Distance semantics** (identical to cgmlst-dists-py): distance = number of loci
where the two alleles differ and **both are present** (non-missing). The
chewBBACA `INF-` prefix is stripped and missing calls are treated as absent.

## Benchmarks

Measured on an 80-core box, computing the **full** matrix (identical output on
all tools; tseemann run with `-x` disabled so it computes the same thing).

| dataset | **cgmlst-dists-rs** | tseemann `cgmlst-dists` (C) | cgmlst-dists-py (CPU) | cgmlst-dists-py (GPU, L4) |
|---------|--------------------:|----------------------------:|----------------------:|--------------------------:|
| 5 000 × 2 000 | **0.41 s** | 3.1 s | 5.0 s | — |
| 50 000 × 5 000 | **~100 s** | 610 s | 782 s | 411 s |

Single-threaded, 2 000 × 2 000: **0.63 s** vs **8.5 s** for the C tool (~13×) —
the C inner loop stays scalar (an early-`return` inside the loop plus per-element
`abs()` block auto-vectorization), while this one runs AVX-512.

## Correctness

Output is verified **byte-identical** to both `cgmlst-dists-py` and
`tseemann/cgmlst-dists` (with its distance cap disabled) on clean data and on
data containing `INF-` prefixes and missing calls.

## Credits

- Original C tool: [tseemann/cgmlst-dists](https://github.com/tseemann/cgmlst-dists)
- Python implementation: [genpat-it/cgmlst-dists-py](https://github.com/genpat-it/cgmlst-dists-py)

## License

GPL-3.0-only
