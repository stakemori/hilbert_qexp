use std::ops::{Mul, Add, Sub, Neg, Shr, MulAssign};

fn int_sqrt(x: u64) -> u64 {
    debug_assert!(x < 1 << 53);
    (x as f64).sqrt().floor() as u64
}

pub trait PowGen {
    fn set_one(&mut self);
    fn square(&mut self);
}

pub fn pow_mut<T>(res: &mut T, f: &T, a: usize)
where
    for<'a> T: MulAssign<&'a T>,
    T: Clone + PowGen,
{
    res.set_one();
    let s = format!("{:b}", a);
    let bts = s.into_bytes();
    let strs: Vec<char> = bts.iter().rev().map(|&i| i as char).collect();
    let mut tmp = f.clone();
    for &c in strs.iter() {
        if c == '0' {
            tmp.square();
        } else if c == '1' {
            T::mul_assign(res, &tmp);
            tmp.square();
        }
    }
}

pub fn prime_sieve(n: usize) -> Vec<usize> {
    let mut vec = Vec::with_capacity(n + 1);
    for _ in 0..(n + 1) {
        vec.push(true);
    }
    let bd = (n as f64).sqrt().floor() as usize;
    for i in 2..(bd + 1) {
        if vec[i] {
            let mut j = i * i;
            while j <= n {
                vec[j] = false;
                j += i;
            }
        }
    }
    let mut res = vec.iter()
        .enumerate()
        .filter(|&(_, &bl)| bl)
        .map(&|(i, _)| i as usize)
        .collect::<Vec<usize>>();
    res.remove(0);
    res.remove(0);
    res
}
