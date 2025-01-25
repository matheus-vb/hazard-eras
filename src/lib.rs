#![allow(dead_code)]

use core::slice;
use std::{
    cell::UnsafeCell,
    sync::atomic::{fence, AtomicPtr, AtomicU64, Ordering},
    usize,
};

use array_init::array_init;

const NONE: u64 = 0;
const HE_MAX_THREADS: usize = 32;
const MAX_HES: usize = 5;
const CLPAD: usize = 32 / std::mem::size_of::<*mut AtomicU64>();
const HE_THRESHOLD_R: i64 = 0;

pub trait Node {
    fn get_new_era(&self) -> u64;
    fn set_new_era(&mut self, era: u64);

    fn get_del_era(&self) -> u64;
    fn set_del_era(&mut self, era: u64);
}

#[repr(align(128))]
struct AlignedAtomicU64(AtomicU64);

#[repr(align(128))]
struct AlignedHe([*mut AtomicU64; HE_MAX_THREADS]);

#[repr(align(128))]
pub struct ThreadLocalRetired<T>(pub UnsafeCell<Vec<*mut T>>);

#[repr(align(128))]
pub struct AlignedRetiredList<T>(pub [ThreadLocalRetired<T>; HE_MAX_THREADS * CLPAD]);

pub struct HazardEras<T>
where
    T: Node + Send + Sync + 'static,
{
    max_hes: usize,
    max_threads: usize,
    era_clock: AlignedAtomicU64,
    he: AlignedHe,
    retired: AlignedRetiredList<T>,
}

impl<T> HazardEras<T>
where
    T: Node + Send + Sync + 'static,
{
    pub fn new(max_hes: usize, max_threads: usize) -> Self {
        let era_clock = AlignedAtomicU64(AtomicU64::new(1));

        let he = AlignedHe(array_init(|_| {
            let he_counters = Box::new([const { AtomicU64::new(NONE) }; CLPAD * 2]);

            Box::into_raw(he_counters) as *mut AtomicU64
        }));

        let retired = AlignedRetiredList(array_init(|_| {
            ThreadLocalRetired(UnsafeCell::new(Vec::with_capacity(max_threads * max_hes)))
        }));

        HazardEras {
            max_hes,
            max_threads,
            era_clock,
            he,
            retired,
        }
    }

    #[inline]
    pub fn get_era(&self) -> u64 {
        self.era_clock.0.load(Ordering::SeqCst)
    }

    pub fn clear(&self, tid: usize) {
        assert!(tid < HE_MAX_THREADS, "Invalid thread id");

        let he_ptr = self.he.0[tid];

        let he_slice = unsafe { slice::from_raw_parts(he_ptr, CLPAD * 2) };

        for ihe in 0..self.max_hes {
            he_slice[ihe].store(NONE, Ordering::Release);
        }
    }

    pub fn get_protected(&self, index: usize, atom: &AtomicPtr<T>, tid: usize) -> *mut T {
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

            he_slice[index].store(era, Ordering::SeqCst);
            prev_era = era;
        }
    }

    pub fn protect_era_release(&self, index: usize, other: usize, tid: usize) {
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

    pub fn protect_ptr(
        &self,
        index: usize,
        atom: &AtomicPtr<T>,
        prev_era: &mut u64,
        tid: usize,
    ) -> *mut T {
        assert!(tid < HE_MAX_THREADS, "Invalid thread id");
        assert!(index < self.max_hes, "Invalid hazard era index");

        let ptr = atom.load(Ordering::Acquire);
        let era = self.era_clock.0.load(Ordering::SeqCst);

        if era != *prev_era {
            *prev_era = era;

            let he_ptr = self.he.0[tid];
            let he_slice = unsafe { slice::from_raw_parts(he_ptr, CLPAD * 2) };

            he_slice[index].store(era, Ordering::Relaxed);
            fence(Ordering::SeqCst);
        }

        ptr
    }

    pub fn retire(&self, ptr: *mut T, my_tid: usize) {
        let curr_era = self.era_clock.0.load(Ordering::SeqCst);

        unsafe {
            (*ptr).set_del_era(curr_era);
        }

        let local_vec = unsafe { &mut *self.retired.0[my_tid * CLPAD].0.get() };
        local_vec.push(ptr);

        if curr_era == self.era_clock.0.load(Ordering::SeqCst) {
            self.era_clock.0.fetch_add(1, Ordering::SeqCst);
        }

        let mut iret = 0;
        while iret < local_vec.len() {
            let obj = local_vec[iret];

            if self.can_delete(obj) {
                local_vec.swap_remove(iret);
                continue;
            }

            iret += 1;
        }
    }

    fn can_delete(&self, obj: *mut T) -> bool {
        for tid in 0..self.max_threads {
            for ihe in 0..self.max_hes {
                let he_ptr = self.he.0[tid];
                let he_slice = unsafe { slice::from_raw_parts(he_ptr, CLPAD * 2) };

                let era = he_slice[ihe].load(Ordering::Acquire);
                if era == NONE
                    || era < unsafe { &*obj }.get_new_era()
                    || era > unsafe { &*obj }.get_del_era()
                {
                    continue;
                }

                return false;
            }
        }

        true
    }
}

impl<T> Drop for HazardEras<T>
where
    T: Node + Send + Sync + 'static,
{
    fn drop(&mut self) {
        for &he_ptr in &self.he.0 {
            if !he_ptr.is_null() {
                unsafe {
                    let array_ptr = he_ptr as *mut [AtomicU64; CLPAD * 2];

                    let boxed_array = Box::from_raw(array_ptr);
                    drop(boxed_array);
                }
            }
        }

        for retired_list in &mut self.retired.0 {
            while let Some(retired_obj) = retired_list.0.get_mut().pop() {
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
    use std::sync::atomic::AtomicPtr;

    /// A test implementation of the `Node` trait.
    /// This struct allows us to manipulate `new_era` and `del_era` values for testing.
    struct TestNode {
        new_era: u64,
        del_era: u64,
    }

    impl Node for TestNode {
        fn get_new_era(&self) -> u64 {
            self.new_era
        }

        fn set_new_era(&mut self, era: u64) {
            self.new_era = era;
        }

        fn get_del_era(&self) -> u64 {
            self.del_era
        }

        fn set_del_era(&mut self, era: u64) {
            self.del_era = era;
        }
    }

    /// Test the initialization of HazardEras.
    #[test]
    fn test_new_hazard_eras() {
        let hazard_eras = HazardEras::<TestNode>::new(MAX_HES, HE_MAX_THREADS);

        // Check that era_clock is initialized to 1
        assert_eq!(hazard_eras.get_era(), 1);

        // Check that all hazard era pointers are non-null and initialized to NONE
        for &he_ptr in &hazard_eras.he.0 {
            assert!(!he_ptr.is_null(), "HE pointer should not be null");
            unsafe {
                let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
                for ihe in 0..MAX_HES {
                    assert_eq!(
                        he_slice[ihe].load(Ordering::Relaxed),
                        NONE,
                        "HE counter should be initialized to NONE"
                    );
                }
            }
        }

        // Check that each retired list has the correct capacity
        for retired in &hazard_eras.retired.0 {
            unsafe {
                let vec = &*retired.0.get();
                assert_eq!(
                    vec.capacity(),
                    MAX_HES * HE_MAX_THREADS,
                    "Retired list capacity mismatch"
                );
            }
        }
    }

    #[test]
    fn test_get_era() {
        let hazard_eras = HazardEras::<TestNode>::new(MAX_HES, HE_MAX_THREADS);
        assert_eq!(hazard_eras.get_era(), 1);

        // Increment era_clock and verify
        hazard_eras.era_clock.0.fetch_add(1, Ordering::SeqCst);
        assert_eq!(hazard_eras.get_era(), 2);
    }

    #[test]
    fn test_clear() {
        let hazard_eras = HazardEras::<TestNode>::new(MAX_HES, HE_MAX_THREADS);
        let tid = 0;

        // Initially, all hazard eras should be NONE
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            for ihe in 0..MAX_HES {
                assert_eq!(he_slice[ihe].load(Ordering::Relaxed), NONE);
            }
        }

        // Set hazard eras to a specific value
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts_mut(he_ptr, CLPAD * 2);
            for ihe in 0..MAX_HES {
                he_slice[ihe].store(42, Ordering::Relaxed);
            }
        }

        // Verify that hazard eras have been set
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            for ihe in 0..MAX_HES {
                assert_eq!(he_slice[ihe].load(Ordering::Relaxed), 42);
            }
        }

        // Clear hazard eras
        hazard_eras.clear(tid);

        // Verify that hazard eras have been reset to NONE
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            for ihe in 0..MAX_HES {
                assert_eq!(he_slice[ihe].load(Ordering::Relaxed), NONE);
            }
        }
    }

    #[test]
    fn test_get_protected() {
        let hazard_eras = HazardEras::<TestNode>::new(MAX_HES, HE_MAX_THREADS);
        let tid = 0;
        let index = 0;

        // Create a TestNode and convert it into a raw pointer
        let node = Box::new(TestNode {
            new_era: 1,
            del_era: 2,
        });
        let node_ptr = Box::into_raw(node);
        let atom = AtomicPtr::new(node_ptr);

        // Initially, era_clock is 1 and he[tid][index] is NONE (0)
        // get_protected should set he[tid][index] to 1 and return the pointer
        let protected_ptr = hazard_eras.get_protected(index, &atom, tid);
        assert_eq!(protected_ptr, node_ptr, "Protected pointer mismatch");

        // Verify that he[tid][index] has been updated to 1
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            assert_eq!(
                he_slice[index].load(Ordering::SeqCst),
                1,
                "HE counter should be updated to 1"
            );
        }

        // Call get_protected again; since era hasn't changed, it should return immediately
        let protected_ptr2 = hazard_eras.get_protected(index, &atom, tid);
        assert_eq!(
            protected_ptr2, node_ptr,
            "Protected pointer mismatch on second call"
        );

        // Change era_clock to 2
        hazard_eras.era_clock.0.store(2, Ordering::SeqCst);

        // Now, get_protected should update he[tid][index] to 2 and return the pointer
        let protected_ptr3 = hazard_eras.get_protected(index, &atom, tid);
        assert_eq!(
            protected_ptr3, node_ptr,
            "Protected pointer mismatch after era change"
        );

        // Verify that he[tid][index] has been updated to 2
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            assert_eq!(
                he_slice[index].load(Ordering::SeqCst),
                2,
                "HE counter should be updated to 2"
            );
        }

        // Clean up by converting the raw pointer back into a Box to prevent memory leaks
        unsafe {
            let value = Box::from_raw(node_ptr);
            drop(value);
        }
    }

    #[test]
    fn test_protect_era_release() {
        let hazard_eras = HazardEras::<TestNode>::new(MAX_HES, HE_MAX_THREADS);
        let tid = 0;
        let index = 0;
        let other = 1;

        // Initialize he[tid][other] to a specific era (e.g., 5)
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts_mut(he_ptr, CLPAD * 2);
            he_slice[other].store(5, Ordering::Relaxed);
            he_slice[index].store(3, Ordering::Relaxed);
        }

        // Call protect_era_release
        hazard_eras.protect_era_release(index, other, tid);

        // Verify that he[tid][index] has been updated to he[tid][other] (5)
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            assert_eq!(
                he_slice[index].load(Ordering::Relaxed),
                5,
                "HE counter should be updated to 5"
            );
        }

        // Call protect_era_release again; since he[tid][index] == he[tid][other], it should remain unchanged
        hazard_eras.protect_era_release(index, other, tid);

        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            assert_eq!(
                he_slice[index].load(Ordering::Relaxed),
                5,
                "HE counter should remain at 5"
            );
        }

        // Update he[tid][other] to a new era (10) and set he[tid][index] to a different value (5)
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts_mut(he_ptr, CLPAD * 2);
            he_slice[other].store(10, Ordering::Relaxed);
            he_slice[index].store(5, Ordering::Relaxed);
        }

        // Call protect_era_release again; it should update he[tid][index] to 10
        hazard_eras.protect_era_release(index, other, tid);

        // Verify that he[tid][index] has been updated to 10
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            assert_eq!(
                he_slice[index].load(Ordering::Relaxed),
                10,
                "HE counter should be updated to 10"
            );
        }
    }

    #[test]
    fn test_protect_ptr() {
        let hazard_eras = HazardEras::<TestNode>::new(MAX_HES, HE_MAX_THREADS);
        let tid = 0;
        let index = 0;

        // Create a TestNode and convert it into a raw pointer
        let node = Box::new(TestNode {
            new_era: 1,
            del_era: 2,
        });
        let node_ptr = Box::into_raw(node);
        let atom = AtomicPtr::new(node_ptr);

        // Initially, era_clock is 1 and prev_era is 0
        let mut prev_era = 0;

        // Call protect_ptr; since era != prev_era, it should update he[tid][index] to 1
        let protected_ptr = hazard_eras.protect_ptr(index, &atom, &mut prev_era, tid);
        assert_eq!(protected_ptr, node_ptr, "Protected pointer mismatch");
        assert_eq!(prev_era, 1, "Previous era should be updated to 1");

        // Verify that he[tid][index] has been updated to 1
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            assert_eq!(
                he_slice[index].load(Ordering::Relaxed),
                1,
                "HE counter should be updated to 1"
            );
        }

        // Change era_clock to 2
        hazard_eras.era_clock.0.store(2, Ordering::SeqCst);

        // Call protect_ptr again; since era != prev_era, it should update he[tid][index] to 2
        let protected_ptr2 = hazard_eras.protect_ptr(index, &atom, &mut prev_era, tid);
        assert_eq!(
            protected_ptr2, node_ptr,
            "Protected pointer mismatch after era change"
        );
        assert_eq!(prev_era, 2, "Previous era should be updated to 2");

        // Verify that he[tid][index] has been updated to 2
        unsafe {
            let he_ptr = hazard_eras.he.0[tid];
            let he_slice = slice::from_raw_parts(he_ptr, CLPAD * 2);
            assert_eq!(
                he_slice[index].load(Ordering::Relaxed),
                2,
                "HE counter should be updated to 2"
            );
        }

        // Clean up by converting the raw pointer back into a Box to prevent memory leaks
        unsafe {
            let value = Box::from_raw(node_ptr);
            drop(value);
        }
    }

    #[test]
    fn test_retire_and_can_delete() {
        let hazard_eras = HazardEras::<TestNode>::new(MAX_HES, HE_MAX_THREADS);
        let tid = 0;

        // Create a TestNode and convert it into a raw pointer
        let node = Box::new(TestNode {
            new_era: 1,
            del_era: 2,
        });
        let node_ptr = Box::into_raw(node);
        let atom = AtomicPtr::new(node_ptr);

        // Protect the node
        let _protected_ptr = hazard_eras.get_protected(0, &atom, tid);

        // Retire the node; this should set del_era to current_era (1)
        hazard_eras.retire(node_ptr, tid);

        // Verify that node_ptr is in the retired list
        unsafe {
            let retired_list = &*hazard_eras.retired.0[tid * CLPAD].0.get();
            assert!(
                retired_list.contains(&node_ptr),
                "Retired list should contain the node pointer"
            );
        }

        // Clear the hazard eras to remove protection
        hazard_eras.clear(tid);

        // Increment era_clock to exceed del_era (now era_clock = 3)
        hazard_eras.era_clock.0.store(3, Ordering::SeqCst);

        // Retire the node again; this should set del_era to current_era (3) and allow deletion
        hazard_eras.retire(node_ptr, tid);

        // Verify that node_ptr has been removed from the retired list
        unsafe {
            let retired_list = &*hazard_eras.retired.0[tid * CLPAD].0.get();
            assert!(
                !retired_list.contains(&node_ptr),
                "Retired list should no longer contain the node pointer"
            );
        }

        // Clean up by converting the raw pointer back into a Box to prevent memory leaks
        unsafe {
            let value = Box::from_raw(node_ptr);
            drop(value);
        }
    }
}
