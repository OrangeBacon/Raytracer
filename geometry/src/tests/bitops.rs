use crate::{is_pow_2, log2u32, log2usize, round_up_pow_2};

#[test]
fn log2() {
    for i in 0..32 {
        assert_eq!(i, log2u32(1 << i));
        assert_eq!(i as usize, log2usize(1 << i));
    }

    for i in 1..31 {
        assert_eq!(i, log2u32((1 << i) + 1));
        assert_eq!(i as usize, log2usize((1 << i) + 1));
    }

    for i in 0..std::mem::size_of::<usize>() {
        assert_eq!(i, log2usize(1 << i));
    }

    for i in 1..std::mem::size_of::<usize>() {
        assert_eq!(i, log2usize((1 << i) + 1));
    }
}

#[test]
fn pow2() {
    for i in 0..32 {
        let i = 1 << i;
        assert!(is_pow_2(i));

        if i > 1 {
            assert!(!is_pow_2(i + 1));
        }
        if i > 2 {
            assert!(!is_pow_2(i - 1));
        }
    }
}

#[test]
fn round_pow_2() {
    assert_eq!(round_up_pow_2(7), 8);

    for i in 1..(1 << 24) {
        if is_pow_2(i) {
            assert_eq!(round_up_pow_2(i), i);
        } else {
            assert_eq!(round_up_pow_2(i), 1 << (log2usize(i) + 1), "{i}");
        }
    }

    for i in 1..62 {
        let v = 1 << i;
        assert_eq!(v, round_up_pow_2(v));
        if v > 2 {
            assert_eq!(v, round_up_pow_2(v - 1));
        }
        assert_eq!(2 * v, round_up_pow_2(v + 1));
    }
}
