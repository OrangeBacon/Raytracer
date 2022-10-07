//! simple crate with no logic used so that the large amounts of data stored
//! directly in the binary are not repeatedly re-compiled

pub mod sobol;

use geometry::Rng;
use once_cell::sync::Lazy;
use proc::prime_sums;

/// Statically calculated random permutations of the radical inverse function
/// for the first 1000 prime numbers
pub static RADICAL_INVERSE_PERMUTATIONS: Lazy<Vec<u64>> = Lazy::new(|| {
    let mut rng = Rng::default();

    let primes = n_primes(1000);

    let mut size = 0;
    for prime in &primes {
        size += prime;
    }

    let mut index = 0;
    let mut permutations = vec![0; size as usize];
    for i in 0..primes.len() {
        for j in 0..primes[i] {
            permutations[(index + j) as usize] = j;
        }
        rng.shuffle_dims(
            &mut permutations[index as usize..(index + primes[i]) as usize],
            1,
        );

        index += primes[i];
    }

    permutations
});

/// Calculate the first COUNT prime numbers.
/// I.e. generates the numbers [2, 3, 5, 7, 11, ...]
pub fn n_primes(count: u64) -> Vec<u64> {
    let mut res = Vec::with_capacity(count as usize);
    let mut current = 2;

    'main: while res.len() < count as usize {
        for i in 2..=((current as f64).sqrt().ceil() as u64) {
            if current % i == 0 && current != i {
                current += 1;
                continue 'main;
            }
        }

        res.push(current);
        current += 1;
    }

    res
}

/// The list of sums of prime numbers below the nth prime, i.e.
/// sum(primes less than 0th prime) = 0, sum(primes less than 1st prime) = 2, ...
pub static PRIME_SUMS: &[usize] = &prime_sums!(1000);
