# ---- build stage: static musl binary with portable AVX2 SIMD ----
FROM rust:1-slim AS build
WORKDIR /src

# musl target -> fully static binary (runs on any x86-64 linux, no glibc needed).
# target-cpu=x86-64-v3 enables AVX2 (portable to any CPU since ~2013) so the
# SIMD Hamming kernel stays vectorized in the shipped image. It is intentionally
# NOT `native` (that would bake in the build machine's CPU and crash elsewhere).
RUN rustup target add x86_64-unknown-linux-musl

COPY Cargo.toml ./
COPY src ./src
RUN RUSTFLAGS="-C target-cpu=x86-64-v3" \
    cargo build --release --target x86_64-unknown-linux-musl \
 && cp target/x86_64-unknown-linux-musl/release/cgmlst-dists /cgmlst-dists

# ---- runtime stage: tiny image with just the static binary ----
FROM alpine:3
COPY --from=build /cgmlst-dists /usr/local/bin/cgmlst-dists
ENTRYPOINT ["cgmlst-dists"]
