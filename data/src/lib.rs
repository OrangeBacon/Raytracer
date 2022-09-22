//! simple crate with no logic used so that the large amounts of data stored
//! directly in the binary are not repeatedly re-compiled

use proc::{radical_permutations, prime_sums};

pub mod sobol;

/// Statically calculated random permutations of the radical inverse function
/// for the first 1000 prime numbers
pub static RADICAL_INVERSE_PERMUTATIONS: &[usize] = &radical_permutations!();

/// The list of sums of prime numbers below the nth prime, i.e.
/// sum(primes less than 0th prime) = 0, sum(primes less than 1st prime) = 2, ...
pub static PRIME_SUMS: &[usize] = &prime_sums!(1000);