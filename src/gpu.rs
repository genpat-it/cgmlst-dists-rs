//! GPU backend (`--features cuda`, `--gpu`): the Hamming distance computed by a
//! hand-wired ReLU network (see cuda/relu_hamming.cu). Allele ids are first
//! re-encoded as dense per-locus ranks (0 = missing, 1..k), which preserves
//! equality within each locus — the only thing the distance depends on — and
//! keeps every value exact in fp32 on the GPU.

use rayon::prelude::*;

/// Dense per-locus ranks of a row-major `n x n_loci` call matrix, as u16.
/// Fails if a locus holds more than 65 535 distinct alleles.
#[cfg_attr(not(feature = "cuda"), allow(dead_code))]
pub fn dense_ranks(calls: &[u32], n: usize, n_loci: usize) -> Result<Vec<u16>, String> {
    let cols: Vec<Vec<u16>> = (0..n_loci)
        .into_par_iter()
        .map(|k| {
            let col: Vec<u32> = (0..n).map(|i| calls[i * n_loci + k]).collect();
            let mut uniq: Vec<u32> = col.iter().copied().filter(|&v| v != 0).collect();
            uniq.sort_unstable();
            uniq.dedup();
            if uniq.len() > u16::MAX as usize {
                return Err(format!("locus {} has {} distinct alleles (GPU limit 65535)", k, uniq.len()));
            }
            Ok(col
                .iter()
                .map(|&v| if v == 0 { 0 } else { uniq.binary_search(&v).unwrap() as u16 + 1 })
                .collect())
        })
        .collect::<Result<_, _>>()?;
    let mut x = vec![0u16; n * n_loci];
    x.par_chunks_mut(n_loci.max(1)).enumerate().for_each(|(i, row)| {
        for (k, v) in row.iter_mut().enumerate() {
            *v = cols[k][i];
        }
    });
    Ok(x)
}

#[cfg(feature = "cuda")]
mod ffi {
    use std::os::raw::{c_char, c_int, c_void};
    pub type StripeCb = extern "C" fn(*mut c_void, c_int, c_int, *const u16);
    extern "C" {
        pub fn rh_device_count() -> c_int;
        pub fn rh_run(
            x: *const u16, n: c_int, l: c_int, device: c_int, stripe_rows: c_int,
            cb: StripeCb, ctx: *mut c_void, gpu_secs: *mut f64, err: *mut c_char, errlen: c_int,
        ) -> c_int;
    }
}

#[cfg(feature = "cuda")]
pub fn device_count() -> usize {
    unsafe { ffi::rh_device_count() as usize }
}

#[cfg(feature = "cuda")]
struct Sink<'a> {
    upper: &'a mut [Vec<i32>],
    n: usize,
    max_dist: i32,
}

#[cfg(feature = "cuda")]
extern "C" fn on_stripe(ctx: *mut std::os::raw::c_void, r0: i32, rows: i32, stripe: *const u16) {
    let sink = unsafe { &mut *(ctx as *mut Sink) };
    let (n, cap, r0, rows) = (sink.n, sink.max_dist, r0 as usize, rows as usize);
    let stripe = unsafe { std::slice::from_raw_parts(stripe, rows * n) };
    sink.upper[r0..r0 + rows].par_iter_mut().enumerate().for_each(|(s, row)| {
        let i = r0 + s;
        *row = stripe[s * n + i + 1..(s + 1) * n]
            .iter()
            .map(|&d| if cap > 0 { (d as i32).min(cap) } else { d as i32 })
            .collect();
    });
}

/// Upper triangle (`upper[i][j-i-1]`, the layout of the CPU path) on the GPU.
/// Returns it with the GPU pipeline time in seconds.
#[cfg(feature = "cuda")]
pub fn upper_triangle(
    calls: &[u32], n: usize, n_loci: usize, max_dist: i32, device: usize, stripe_rows: usize,
) -> Result<(Vec<Vec<i32>>, f64), String> {
    if n_loci > u16::MAX as usize {
        return Err(format!("{} loci: the GPU path supports at most 65535", n_loci));
    }
    if n > i32::MAX as usize {
        return Err("too many samples".into());
    }
    let x = dense_ranks(calls, n, n_loci)?;
    let mut upper: Vec<Vec<i32>> = vec![Vec::new(); n];
    if n < 2 || n_loci == 0 {
        upper.iter_mut().enumerate().for_each(|(i, r)| *r = vec![0; n - 1 - i]);
        return Ok((upper, 0.0));
    }
    let mut sink = Sink { upper: &mut upper, n, max_dist };
    let mut err = vec![0 as std::os::raw::c_char; 512];
    let mut secs = 0f64;
    let rc = unsafe {
        ffi::rh_run(
            x.as_ptr(), n as i32, n_loci as i32, device as i32, stripe_rows as i32, on_stripe,
            &mut sink as *mut Sink as *mut _, &mut secs, err.as_mut_ptr(), err.len() as i32,
        )
    };
    if rc != 0 {
        let msg = unsafe { std::ffi::CStr::from_ptr(err.as_ptr()) }.to_string_lossy().into_owned();
        return Err(format!("CUDA error: {msg}"));
    }
    Ok((upper, secs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranks_preserve_equality_and_missing() {
        // 3 samples x 2 loci, CRC32-sized ids
        let calls = [4_294_967_295, 0, 7, 3_000_000_000, 4_294_967_295, 3_000_000_000];
        let x = dense_ranks(&calls, 3, 2).unwrap();
        assert_eq!(x[0], x[4]); // equal ids, same rank
        assert_ne!(x[0], x[2]);
        assert_eq!(x[1], 0); // missing stays 0
        assert_eq!(x[3], x[5]);
        assert!(x.iter().all(|&v| v <= 2));
    }

    #[cfg(feature = "cuda")]
    #[test]
    fn gpu_matches_cpu() {
        if device_count() == 0 {
            eprintln!("no CUDA device, skipped");
            return;
        }
        // deterministic LCG data; sizes straddle tile edges (64) and force several stripes
        let (n, l) = (517usize, 203usize);
        let mut s: u64 = 42;
        let mut rnd = || { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (s >> 33) as u32 };
        let calls: Vec<u32> = (0..n * l)
            .map(|_| { let r = rnd(); if r % 20 == 0 { 0 } else if r % 3 == 0 { u32::MAX - r % 5 } else { r % 4 } })
            .collect();
        for &(cap, stripe) in &[(0, 0), (0, 64), (5, 128)] {
            let (up, _) = upper_triangle(&calls, n, l, cap, 0, stripe).unwrap();
            for i in 0..n {
                for j in i + 1..n {
                    let cpu = crate::distance(&calls[i * l..(i + 1) * l], &calls[j * l..(j + 1) * l], cap);
                    assert_eq!(up[i][j - i - 1], cpu, "i={i} j={j} cap={cap} stripe={stripe}");
                }
            }
        }
    }
}
