// Builds the CUDA kernel (cuda/hamming.cu) only with `--features cuda`.
// Needs nvcc; CUDA_HOME (default /usr/local/cuda) locates it. The CUDA runtime
// is linked statically, so the binary needs only the NVIDIA driver at run time.
// CUDA_ARCH picks the target GPUs: `native` (default, the GPU in this box) or a
// comma-separated list of compute capabilities such as `75,80,86,89,90`; the
// last one is also embedded as PTX so newer GPUs can JIT-compile it.
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
    let mut arch_flags = Vec::new();
    if arch == "native" {
        arch_flags.push("-arch=native".to_string());
    } else {
        let caps: Vec<&str> = arch.split(',').map(|c| c.trim().trim_start_matches("sm_")).filter(|c| !c.is_empty()).collect();
        assert!(!caps.is_empty(), "CUDA_ARCH is empty");
        for c in &caps {
            arch_flags.push(format!("-gencode=arch=compute_{c},code=sm_{c}"));
        }
        let last = caps[caps.len() - 1];
        arch_flags.push(format!("-gencode=arch=compute_{last},code=compute_{last}"));
    }
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let obj = out.join("hamming.o");
    let st = Command::new(cuda.join("bin/nvcc"))
        .args(["-O3", "-c", "-Xcompiler", "-fPIC"])
        .args(&arch_flags)
        .args(["cuda/hamming.cu", "-o"])
        .arg(&obj)
        .status()
        .expect("nvcc not found (set CUDA_HOME)");
    assert!(st.success(), "nvcc failed");
    let st = Command::new("ar").arg("crs").arg(out.join("libhamming.a")).arg(&obj).status().unwrap();
    assert!(st.success(), "ar failed");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=hamming");
    println!("cargo:rustc-link-search=native={}", cuda.join("lib64").display());
    println!("cargo:rustc-link-lib=static=cudart_static");
    for lib in ["stdc++", "dl", "rt", "pthread"] {
        println!("cargo:rustc-link-lib=dylib={lib}");
    }
}
