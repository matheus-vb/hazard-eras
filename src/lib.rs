#![allow(dead_code)]

use std::{
    sync::atomic::{AtomicU64, Ordering},
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {}
}
