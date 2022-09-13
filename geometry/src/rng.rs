use crate::{Float, Number};

/// Simple xor shift pseudo random number generator used to generate
/// test cases in a deterministic way.
pub struct Rng {
    state: <Float as Number>::Bits,
    counter: u32,
}

impl Rng {
    /// Create a new random number generator
    pub fn new(seed: <Float as Number>::Bits) -> Self {
        Self {
            state: seed,
            counter: seed as _,
        }
    }

    /// Generate a 32 bit random number
    #[cfg(not(feature = "double"))]
    pub fn float_bits(&mut self) -> u32 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        self.state = x;
        self.counter = self.counter.wrapping_add(36247);
        x.wrapping_add(self.counter)
    }

    /// Generate a 64 bit random number
    #[cfg(feature = "double")]
    pub fn float_bits(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        self.counter = self.counter.wrapping_add(36247);
        x.wrapping_add(self.counter as u64)
    }

    /// Generate a new floating point number in [0, 1)
    pub fn float(&mut self) -> Float {
        let num = self.float_bits();

        let float_bits = std::mem::size_of::<Float>() * 8;
        let multiplier = (-(float_bits as Float)).exp2();
        let f = (num as Float) * multiplier;

        (1.0 - Float::EPSILON).min(f)
    }
}
