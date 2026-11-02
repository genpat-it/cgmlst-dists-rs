#!/usr/bin/env python3
"""Generate a synthetic chewBBACA-style allele-call table, to benchmark or
verify cgmlst-dists without shipping large datasets.

Profiles are clonal (samples derive from a few founder profiles with random
allele changes), like surveillance data. --hashed writes CRC32 allele ids as in
chewBBACA hashed profiles: about half are above 2^31, which is what exposes
32-bit signed readers. The distance matrix is the same with or without
--hashed for the same --seed, so the two outputs must be identical.

    python3 scripts/generate_profiles.py -n 20000 -l 3000 -o bench.tsv
    python3 scripts/generate_profiles.py -n 20000 -l 3000 --hashed -o bench_hashed.tsv
"""
import argparse
import random
import zlib


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-n", "--samples", type=int, default=1000)
    ap.add_argument("-l", "--loci", type=int, default=1500)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--founders", type=int, default=50, help="number of founder profiles")
    ap.add_argument("--mutation", type=float, default=0.05, help="per-locus allele change rate")
    ap.add_argument("--missing", type=float, default=0.02, help="per-cell missing rate")
    ap.add_argument("--hashed", action="store_true", help="write CRC32 allele ids")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    founders = [[rng.randint(1, 49) for _ in range(a.loci)] for _ in range(a.founders)]

    def allele(k, v):
        return str(zlib.crc32(f"locus{k}:{v}".encode())) if a.hashed else str(v)

    with open(a.output, "w") as f:
        f.write("FILE\t" + "\t".join(f"locus{k}" for k in range(a.loci)) + "\n")
        for i in range(a.samples):
            base = founders[rng.randrange(a.founders)]
            row = []
            for k, v in enumerate(base):
                if rng.random() < a.mutation:
                    v += rng.randint(1, 29)
                row.append("-" if rng.random() < a.missing else allele(k, v))
            f.write(f"sample{i}\t" + "\t".join(row) + "\n")


if __name__ == "__main__":
    main()
