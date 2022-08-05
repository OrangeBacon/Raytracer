use crate::Float;

/// Simple xor shift pseudo random number generator used to generate
/// test cases in a deterministic way.
pub struct Rng {
    state_32: u32,
    state_64: u64,
    counter: u32,
}

impl Rng {
    /// Create a new random number generator
    pub fn new(seed: u64) -> Self {
        let bytes = seed.to_le_bytes();
        let state_32 = {
            let mut data = [0; 4];
            data.copy_from_slice(&bytes[0..4]);
            u32::from_le_bytes(data)
        };
        Self {
            state_32,
            state_64: seed,
            counter: state_32,
        }
    }

    /// Generate a 32 bit random number
    pub fn u32(&mut self) -> u32 {
        let mut x = self.state_32;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        self.state_32 = x;
        self.counter = self.counter.wrapping_add(36247);
        x.wrapping_add(self.counter)
    }

    /// Generate a 64 bit random number
    pub fn u64(&mut self) -> u64 {
        let mut x = self.state_64;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state_64 = x;
        self.counter = self.counter.wrapping_add(36247);
        x.wrapping_add(self.counter as u64)
    }

    /// Generate a new floating point number in [0, 1)
    pub fn float(&mut self) -> Float {
        #[cfg(feature = "double")]
        let num = self.u64();

        #[cfg(not(feature = "double"))]
        let num = self.u32();

        let float_bits = std::mem::size_of::<Float>() * 8;
        let multiplier = (-(float_bits as Float)).exp2();
        let f = (num as Float) * multiplier;

        (1.0 - Float::EPSILON).min(f)
    }
}
