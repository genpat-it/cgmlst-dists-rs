# ---- build stage (glibc) ----
FROM rust:1-slim AS build
WORKDIR /src

# target-cpu=x86-64-v3 enables AVX2 (portable to any x86-64 CPU since ~2013) so
# the SIMD Hamming kernel stays vectorized in the shipped image. It is
# intentionally NOT `native` (that would bake in the build machine's CPU).
COPY Cargo.toml Cargo.lock build.rs ./
COPY src ./src
RUN RUSTFLAGS="-C target-cpu=x86-64-v3" cargo build --release --locked \
 && cp target/release/cgmlst-dists /cgmlst-dists

# ---- runtime stage ----
FROM debian:stable-slim
COPY --from=build /cgmlst-dists /usr/local/bin/cgmlst-dists
ENTRYPOINT ["cgmlst-dists"]
