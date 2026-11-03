// Builds the CUDA kernel (cuda/hamming.cu) only with `--features cuda`.
// Needs nvcc; CUDA_HOME (default /usr/local/cuda) locates it and libcudart.
// CUDA_ARCH overrides the target GPU (default: `native`, the GPU in this box).
use std::{env, path::PathBuf, process::Command};

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    if env::var_os("CARGO_FEATURE_CUDA").is_none() {
        return;
    }
    println!("cargo:rerun-if-changed=cuda/hamming.cu");
    println!("cargo:rerun-if-env-changed=CUDA_HOME");
    println!("cargo:rerun-if-env-changed=CUDA_ARCH");
    let cuda = PathBuf::from(env::var("CUDA_HOME").unwrap_or_else(|_| "/usr/local/cuda".into()));
    let arch = env::var("CUDA_ARCH").unwrap_or_else(|_| "native".into());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let obj = out.join("hamming.o");
    let st = Command::new(cuda.join("bin/nvcc"))
        .args(["-O3", "-c", "-Xcompiler", "-fPIC", &format!("-arch={arch}"), "cuda/hamming.cu", "-o"])
        .arg(&obj)
        .status()
        .expect("nvcc not found (set CUDA_HOME)");
    assert!(st.success(), "nvcc failed");
    let st = Command::new("ar").arg("crs").arg(out.join("libhamming.a")).arg(&obj).status().unwrap();
    assert!(st.success(), "ar failed");
    let lib = cuda.join("lib64");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=hamming");
    println!("cargo:rustc-link-search=native={}", lib.display());
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib.display());
}
