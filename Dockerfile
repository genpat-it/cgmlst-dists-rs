# ---- build stage (glibc) ----
FROM rust:1-slim AS build
WORKDIR /src

# C toolchain for gzp's native deps (libz-sys / libdeflate)
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential cmake \
 && rm -rf /var/lib/apt/lists/*

# target-cpu=x86-64-v3 enables AVX2 (portable to any x86-64 CPU since ~2013) so
# the SIMD Hamming kernel stays vectorized in the shipped image. It is
# intentionally NOT `native` (that would bake in the build machine's CPU).
COPY Cargo.toml ./
COPY src ./src
RUN RUSTFLAGS="-C target-cpu=x86-64-v3" cargo build --release \
 && cp target/release/cgmlst-dists /cgmlst-dists

# ---- runtime stage ----
# Debian slim (glibc). The binary statically bundles its C deps (zlib/libdeflate
# via gzp), so no extra runtime packages are needed.
FROM debian:stable-slim
COPY --from=build /cgmlst-dists /usr/local/bin/cgmlst-dists
ENTRYPOINT ["cgmlst-dists"]
