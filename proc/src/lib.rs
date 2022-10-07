use proc_macro::TokenStream;

/// Calculate the radical inverse of a number.  Should be used through
/// [`geometry::radical_inverse`], not this proc_macro.
#[proc_macro]
pub fn radical_inverse(item: TokenStream) -> TokenStream {
    let mut input = item.into_iter();
    let function = input
        .next()
        .expect("Expected radical_index_base function name");
    input.next().expect("Expected comma");
    let base_index = input.next().expect("Expected base index variable name");
    input.next().expect("Expected comma");
    let number = input
        .next()
        .expect("Expected number to process variable name");

    let mut lines = Vec::with_capacity(1027);
    lines.push(format!("match {} {{", base_index));
    lines.push(format!(
        "0 => T::cast({}.reverse_bits() as f64) * (T::ONE / T::TWO.pow(T::cast(64))),",
        number
    ));

    for (idx, prime) in n_primes(1024).into_iter().enumerate().skip(1) {
        lines.push(format!(
            "{} => {}::<{}, _>({}),",
            idx, function, prime, number
        ));
    }

    lines.push("_ => panic!(\"Maximum base for radical inverse is 1023\")".into());
    lines.push("}".into());

    let lines = lines.join("");

    lines.parse().unwrap()
}

/// Calculate the first COUNT prime numbers.
/// I.e. generates the numbers [2, 3, 5, 7, 11, ...]
fn n_primes(count: usize) -> Vec<usize> {
    let mut res = Vec::with_capacity(count);
    let mut current = 2;

    'main: while res.len() < count {
        for i in 2..=((current as f64).sqrt().ceil() as usize) {
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

#[proc_macro]
pub fn scrambled_radical_inverse(item: TokenStream) -> TokenStream {
    let mut input = item.into_iter();
    let function = input
        .next()
        .expect("Expected scrambled_radical_index_base function name");
    input.next().expect("Expected comma");
    let base_index = input.next().expect("Expected base index variable name");
    input.next().expect("Expected comma");
    let number = input
        .next()
        .expect("Expected number to process variable name");
    input.next().expect("Expected comma");
    let permutation = input.next().expect("Expected permutation variable name");

    let mut lines = Vec::with_capacity(1027);
    lines.push(format!("match {} {{", base_index));

    for (idx, prime) in n_primes(1024).into_iter().enumerate() {
        lines.push(format!(
            "{} => {}::<{}, _>({}, {}),",
            idx, function, prime, number, permutation
        ));
    }

    lines.push("_ => panic!(\"Maximum base for scrambled radical inverse is 1023\")".into());
    lines.push("}".into());

    let lines = lines.join("");

    lines.parse().unwrap()
}

/// Generate a list of the sums of all prime numbers below each prime number.
/// 0 => primes below 2, 1 => primes below 3, 2 => primes below 5, ...
#[proc_macro]
pub fn prime_sums(item: TokenStream) -> TokenStream {
    let count = item
        .into_iter()
        .next()
        .expect("Expected number of primes to generate");
    let count = count
        .to_string()
        .parse::<usize>()
        .expect("Unable to parse count as number");

    let primes = n_primes(count);

    let prime_sums = primes
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, _)| primes[0..idx].iter().sum())
        .collect::<Vec<usize>>();

    format!("{:?}", prime_sums).parse().unwrap()
}
