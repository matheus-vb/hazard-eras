#![allow(dead_code)]

use core::slice;
use std::{
    slice,
    sync::atomic::{AtomicPtr, AtomicU64, Ordering},
    usize,
};

use array_init::array_init;

const NONE: u64 = 0;
const HE_MAX_THREADS: usize = 128;
const MAX_HES: usize = 5;
const CLPAD: usize = 128 / std::mem::size_of::<*mut AtomicU64>();
const HE_THRESHOLD_R: i64 = 0;

pub trait Node {
    fn new_era(&self);
    fn del_era(&self);
}

#[repr(align(128))]
struct AlignedAtomicU64(AtomicU64);

#[repr(align(128))]
struct AlignedHe([*mut AtomicU64; HE_MAX_THREADS]);

#[repr(align(128))]
struct AlignedRetiredList<T>([Vec<*mut T>; HE_MAX_THREADS * CLPAD]);

pub struct HazardEras<T>
where
    T: Node,
{
    max_hes: usize,
    max_threads: usize,
    era_clock: AlignedAtomicU64,
    he: AlignedHe,
    retired: AlignedRetiredList<T>,
}

impl<T> HazardEras<T>
where
    T: Node,
{
    pub fn new(max_hes: usize, max_threads: usize) -> Self {
        let era_clock = AlignedAtomicU64(AtomicU64::new(1));

        let he = AlignedHe(array_init(|_| {
            let he_counters = Box::new([const { AtomicU64::new(NONE) }; CLPAD * 2]);

            Box::into_raw(he_counters) as *mut AtomicU64
        }));

        let retired = AlignedRetiredList(array_init(|_| Vec::with_capacity(max_threads * max_hes)));

        HazardEras {
            max_hes,
            max_threads,
            era_clock,
            he,
            retired,
        }
    }

    #[inline]
    fn get_era(&self) -> u64 {
        self.era_clock.0.load(Ordering::SeqCst)
    }

    fn clear(&self, tid: usize) {
        assert!(tid < HE_MAX_THREADS, "Invalid thread id");

        let he_ptr = self.he.0[tid];

        let he_slice = unsafe { slice::from_raw_parts(he_ptr, CLPAD * 2) };

        for ihe in 0..self.max_hes {
            he_slice[ihe].store(NONE, Ordering::Release);
        }
    }

    fn get_protected(&self, index: usize, atom: &AtomicPtr<T>, tid: usize) -> *mut T {
        assert!(tid < HE_MAX_THREADS, "Invalid thread id");
        assert!(index < self.max_hes, "Invalid hazard era index");

        let he_ptr = self.he.0[tid];

        let he_slice = unsafe { slice::from_raw_parts(he_ptr, CLPAD * 2) };

        let mut prev_era = he_slice[index].load(Ordering::Relaxed);

        loop {
            let ptr = atom.load(Ordering::Acquire);
            let era = self.era_clock.0.load(Ordering::Acquire);

            if era == prev_era {
                return ptr;
            }

            he_slice[index].store(era, Ordering::Release);
            prev_era = era;
        }
    }

    fn protect_era_release(&self, index: usize, other: usize, tid: usize) {
        assert!(tid < HE_MAX_THREADS, "Invalid thread id");
        assert!(index < self.max_hes, "Invalid hazard era index");
        assert!(other < self.max_hes, "Invalid hazard era index");

        let he_ptr = self.he.0[tid];
        let he_slice = unsafe { slice::from_raw_parts(he_ptr, CLPAD * 2) };

        let era = he_slice[other].load(Ordering::Relaxed);

        if he_slice[index].load(Ordering::Relaxed) == era {
            return;
        }

        he_slice[index].store(era, Ordering::Release);
    }
}

impl<T> Drop for HazardEras<T>
where
    T: Node,
{
    fn drop(&mut self) {
        for &he_ptr in &self.he.0 {
            if !he_ptr.is_null() {
                unsafe {
                    let value = Box::from_raw(he_ptr);
                    drop(value);
                }
            }
        }

        for retired_list in &mut self.retired.0 {
            while let Some(retired_obj) = retired_list.pop() {
                unsafe {
                    let value = Box::from_raw(retired_obj);
                    drop(value);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        //let he = HazardEras::new(10, 12);
    }
}
