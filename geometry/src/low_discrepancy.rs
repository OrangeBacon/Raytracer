use proc::{primes, radical_inverse, scrambled_radical_inverse};

use crate::{Number, Rng};

/// Calculate the radical inverse of a number in the given base, works for any
/// base where the base can be chosen at compile time.
pub fn radical_inverse_base<const BASE: u64, T: Number>(mut num: u64) -> T {
    let inv = T::ONE / T::cast(BASE as f64);
    let mut digits = 0;
    let mut inv_base_n = T::ONE;

    while num > 0 {
        let next = num / BASE;
        let digit = num - next * BASE;
        digits = digits * BASE + digit;
        inv_base_n *= inv;
        num = next;
    }

    T::ONE_MINUS_EPSILON.min(T::cast(digits as f64) * inv_base_n)
}

/// Compute the base n radical inverse of a given number.
/// Will panic if the base index is greater than 1023.
pub fn radical_inverse<T: Number>(base_index: usize, num: u64) -> T {
    radical_inverse!(radical_inverse_base, base_index, num)
}

/// Compute 1/radical_inverse in the given base, i.e. the integer portion of the inverse.
/// The number of digits in the original value should be known to give the right result.
pub fn inverse_radical_inverse<const BASE: u64>(mut inverse: u64, n_digits: usize) -> u64 {
    let mut index = 0;
    for _ in 0..n_digits {
        let digit = inverse % BASE;
        inverse /= BASE;
        index = index * BASE + digit;
    }

    index
}

/// The first 1000 prime numbers
const PRIMES: [usize; 1000] = primes!(1000);

/// Compute the radical inverse permutation tables
pub fn radical_inverse_permutations(rng: &mut Rng) -> Vec<usize> {
    let mut size = 0;
    for prime in PRIMES {
        size += prime;
    }

    let mut index = 0;
    let mut permutations = Vec::with_capacity(size);
    for i in 0..PRIMES.len() {
        for j in 0..PRIMES[i] {
            permutations[index + j] = j;
        }
        rng.shuffle_dims(&mut permutations[index..index + PRIMES[i]], 1);

        index += PRIMES[i];
    }

    permutations
}

/// Compute the scrambled base n radical inverse of a given number.
/// Works for any base that can be provided at compile time.
pub fn scrambled_radical_inverse_base<const BASE: u64, T: Number>(
    mut num: u64,
    permutation: &[usize],
) -> T {
    let inv = T::ONE / T::cast(BASE as f64);
    let mut digits = 0;
    let mut inv_base_n = T::ONE;

    while num > 0 {
        let next = num / BASE;
        let digit = num - next * BASE;
        digits = digits * BASE + permutation[digit as usize] as u64;
        inv_base_n *= inv;
        num = next;
    }

    T::ONE_MINUS_EPSILON.min(
        inv_base_n
            * (T::cast(digits as f64) + inv * T::cast(permutation[0] as f64) / (T::ONE - inv)),
    )
}

/// Compute the scrambled base n radical inverse of a given number.
/// Will panic if the base index is greater than 1023.
pub fn scrambled_radical_inverse<T: Number>(
    base_index: usize,
    num: u64,
    permutation: &[usize],
) -> T {
    scrambled_radical_inverse!(scrambled_radical_inverse_base, base_index, num, permutation)
}
