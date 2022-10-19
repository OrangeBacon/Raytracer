use std::sync::atomic::{AtomicI32, AtomicI64, AtomicU32, AtomicU64, Ordering};

use crate::Number;

pub trait AtomicNumber<T: Number>: Default {
    fn from_float(v: T) -> Self;

    fn to_float(&self) -> T;

    fn add(&self, v: T);
}

impl AtomicNumber<i32> for AtomicI32 {
    fn from_float(v: i32) -> Self {
        AtomicI32::new(v)
    }

    fn to_float(&self) -> i32 {
        self.load(Ordering::SeqCst)
    }

    fn add(&self, v: i32) {
        self.fetch_add(v, Ordering::SeqCst);
    }
}

impl AtomicNumber<i64> for AtomicI64 {
    fn from_float(v: i64) -> Self {
        AtomicI64::new(v)
    }

    fn to_float(&self) -> i64 {
        self.load(Ordering::SeqCst)
    }

    fn add(&self, v: i64) {
        self.fetch_add(v, Ordering::SeqCst);
    }
}

#[derive(Debug, Default)]
pub struct AtomicF32 {
    bits: AtomicU32,
}
impl AtomicNumber<f32> for AtomicF32 {
    /// Create a new 32 bit atomic float
    fn from_float(f: f32) -> Self {
        Self {
            bits: AtomicU32::new(f.to_bits()),
        }
    }

    fn to_float(&self) -> f32 {
        f32::from_bits(self.bits.load(Ordering::SeqCst))
    }

    fn add(&self, v: f32) {
        self.bits
            .fetch_update(Ordering::SeqCst, Ordering::SeqCst, |f| {
                Some((f32::from_bits(f) + v).to_bits())
            })
            .unwrap(); // unwrap because there isn't anything to fail in the function
    }
}

#[derive(Debug, Default)]
pub struct AtomicF64 {
    bits: AtomicU64,
}
impl AtomicNumber<f64> for AtomicF64 {
    /// Create a new 32 bit atomic float
    fn from_float(f: f64) -> Self {
        Self {
            bits: AtomicU64::new(f.to_bits()),
        }
    }

    fn to_float(&self) -> f64 {
        f64::from_bits(self.bits.load(Ordering::SeqCst))
    }

    fn add(&self, v: f64) {
        self.bits
            .fetch_update(Ordering::SeqCst, Ordering::SeqCst, |f| {
                Some((f64::from_bits(f) + v).to_bits())
            })
            .unwrap(); // unwrap because there isn't anything to fail in the function
    }
}
