use std::num::Wrapping;

use crate::Number;

const STATE: Wrapping<u64> = Wrapping(0x853c49e6748fea9b);
const STREAM: Wrapping<u64> = Wrapping(0xda3e39cb94b95bdb);
const MULTIPLY: Wrapping<u64> = Wrapping(0x5851f42d4c957f2d);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Rng {
    state: Wrapping<u64>,
    inc: Wrapping<u64>,
}

impl Default for Rng {
    fn default() -> Self {
        Self {
            state: STATE,
            inc: STREAM,
        }
    }
}

impl Rng {
    /// Create a rng with the given seed
    pub fn new(sequence_index: u64) -> Self {
        let mut this = Self::default();
        this.set_sequence(sequence_index);
        this
    }

    /// Set the seed of the random number generator
    pub fn set_sequence(&mut self, index: u64) {
        self.state = Wrapping(0);
        self.inc = (Wrapping(index) << 1) | Wrapping(1);
        self.uniform_u32();
        self.state += STATE;
        self.uniform_u32();
    }

    /// Generate a new uniformly distributed u32 random number
    pub fn uniform_u32(&mut self) -> u32 {
        let old = self.state;
        self.state = old * MULTIPLY + self.inc;
        let shift = Wrapping((((old >> 18) ^ old) >> 27).0 as u32);
        let rot = Wrapping((old >> 59).0 as u32);
        ((shift >> rot.0 as usize) | (shift << ((!rot + Wrapping(1)) & Wrapping(31)).0 as usize)).0
    }

    /// Generate a new uniformly distributed random u32 in the range [0, limit)
    pub fn uniform_u32_limit(&mut self, limit: u32) -> u32 {
        let threshold = (!limit + 1) % limit;
        loop {
            let r = self.uniform_u32();
            if r >= threshold {
                return r % limit;
            }
        }
    }

    /// Generate a uniformly distributed floating point number, from 32 bits of
    /// random data.  The output will be in the range [0, 1).
    pub fn uniform_float<T: Number>(&mut self) -> T {
        let step = 2.0f32.powi(-32);
        let num = self.uniform_u32() as f32 * step;
        T::ONE_MINUS_EPSILON.min(T::cast(num))
    }

    /// Randomly shuffle a given range of data
    pub fn shuffle<T: Copy>(&mut self, data: &mut [T]) {
        for idx in (0..data.len()).rev() {
            data.swap(idx, self.uniform_u32_limit(data.len() as _) as _)
        }
    }
}
