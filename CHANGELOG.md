# Changelog

## [0.1.3] - 2026-11-02

### Fixed
- **chewBBACA hashed profiles (CRC32 allele ids) were read wrong.** Alleles were
  parsed as i32: every hash above 2^31-1 (about half of a hashed profile) failed
  the parse and was silently treated as missing, collapsing distances (Listeria
  Pasteur cgMLST, 2099 genomes: 282 instead of 1196 for one pair; 2 199 325 of
  2 201 851 pairs wrong). Alleles are now parsed as u64 and stored as u32 (loci
  with larger ids are re-encoded per locus), so the SIMD kernel and speed are
  unchanged. Output is md5-identical to tseemann `cgmlst-dists` 0.6.0 and to
  cgmlst-dists-py 0.1.8 on these data. Integers beyond u64 are an error instead
  of missing.

### Added
- `-N/--novel-as-allele`: count chewBBACA `*N` novel alleles as alleles, compared
  within each locus (default unchanged: missing). Use only on tables from a
  single chewBBACA run.
- `scripts/generate_profiles.py`: synthetic (optionally hashed) profiles to
  benchmark and verify without shipped datasets.
- GitHub Pages site and a Bioconda recipe (`recipes/bioconda/`, binary named
  `cgmlst-dists-rs` in the package).
