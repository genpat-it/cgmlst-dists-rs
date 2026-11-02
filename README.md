# cgmlst-dists-rs

[![Crates.io](https://img.shields.io/crates/v/cgmlst-dists-rs.svg?logo=rust)](https://crates.io/crates/cgmlst-dists-rs)
[![CI](https://github.com/genpat-it/cgmlst-dists-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/genpat-it/cgmlst-dists-rs/actions/workflows/ci.yml)
[![Docker](https://img.shields.io/badge/docker-ghcr.io-blue?logo=docker)](https://github.com/genpat-it/cgmlst-dists-rs/pkgs/container/cgmlst-dists-rs)
[![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)

Ultra-fast pairwise **cgMLST Hamming distance matrix** from chewBBACA allele
call tables — a Rust port of
[cgmlst-dists](https://github.com/tseemann/cgmlst-dists) by Torsten Seemann (via
[cgmlst-dists-py](https://github.com/genpat-it/cgmlst-dists-py)), producing
**byte-identical output** while running several times faster.

> **Independent port.** This is not the original tool and is not affiliated with
> or endorsed by its author. `cgmlst-dists` is GPL-3.0, and so is this port.

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

### Bioconda

```bash
conda install -c bioconda cgmlst-dists-rs
```

In the conda package the binary is named **`cgmlst-dists-rs`**, so it can sit in
the same environment as the original `cgmlst-dists`.

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
| `-z, --gzip` | gzip the output (also auto-enabled when `-o` ends with `.gz`) |
| `-t, --threads` | threads [default: all cores] |
| `-m, --matrix-format` | `full` \| `lower-tri` \| `upper-tri` |
| `-d, --input-sep` / `-D, --output-sep` | field separators [tab] |
| `-x, --index-name` | top-left cell name [`cgmlst-dists`] |
| `-n, --missing-char` | missing-data character [`-`] |
| `-r, --skip-input-replacements` | do not strip the chewBBACA `INF-` prefix |
| `-N, --novel-as-allele` | count `*N` novel alleles as alleles (default: missing) |
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

Allele ids are read as 64-bit unsigned integers, so **chewBBACA hashed profiles**
(CRC32 hashes up to 4 294 967 295) work. Non-numeric tokens — chewBBACA classes
(`LNF`, `PLOT3`, `ASM`, …) and `*N` novel alleles — are treated as missing, as in
tseemann `cgmlst-dists`.

`*N` (e.g. `*486`) is an allele chewBBACA inferred but did not add to the schema.
With `-N/--novel-as-allele` these are compared as alleles within each locus (same
`*N` = same allele), which recovers real differences. The `*N` numbering is local
to one chewBBACA run, so use it only on a table from a single run.

## Benchmarks

Full matrix on an 80-core x86-64 box with AVX-512; both tools produce
md5-identical output. tseemann `cgmlst-dists` 0.6.0 is built with its Makefile
(`-Ofast -fopenmp`) and run with `-j` = core count and `-x 100000` (distance cap
disabled, so it computes the same thing). Best of three.

| dataset | threads | **cgmlst-dists-rs** | tseemann `cgmlst-dists` (C) |
|---------|--------:|--------------------:|----------------------------:|
| 2 099 × 1 748 (real Listeria, hashed) | 1 | **0.72 s** | 13.2 s |
| 2 099 × 1 748 (real Listeria, hashed) | 80 | **0.06 s** | 1.0 s |
| 20 000 × 3 000 (synthetic) | 80 | **5.2 s** | 106 s |

The C inner loop stays scalar (an early `return` inside the loop plus
per-element `abs()` block auto-vectorization), while this one runs AVX-512.

### Reproduce

No large datasets are shipped; generate them:

```bash
python3 scripts/generate_profiles.py -n 20000 -l 3000 -o bench.tsv
python3 scripts/generate_profiles.py -n 20000 -l 3000 --hashed -o bench_hashed.tsv
time cgmlst-dists -i bench.tsv -o rs.tsv
time cgmlst-dists -i bench_hashed.tsv -o rs_hashed.tsv
md5sum rs.tsv rs_hashed.tsv     # identical: same profiles, CRC32 ids
```

## Correctness

Output is verified **byte-identical** to both `cgmlst-dists-py` and
`tseemann/cgmlst-dists` (with its distance cap disabled) on clean data and on
data containing `INF-` prefixes and missing calls.
On chewBBACA hashed profiles it matches `tseemann/cgmlst-dists` on all pairs;
cgmlst-dists-py ≤ 0.1.7 mis-reads hashes above 2^31 (fixed in its next release).

## Credits

- Original C tool: [tseemann/cgmlst-dists](https://github.com/tseemann/cgmlst-dists)
- Python implementation: [genpat-it/cgmlst-dists-py](https://github.com/genpat-it/cgmlst-dists-py)

## License

GPL-3.0-only
