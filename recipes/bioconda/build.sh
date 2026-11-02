#!/bin/bash
set -ex

# Portable build: no -C target-cpu=native in a package (it would crash on older CPUs)
cargo install --no-track --locked --root "$PREFIX" --path .

# Bundle third-party licenses
cargo-bundle-licenses --format yaml --output THIRDPARTY.yml
