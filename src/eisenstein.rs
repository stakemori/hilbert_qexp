use elements::{HmfGen, UBounds};
use misc::prime_sieve;
use bignum::BigNumber;
use flint::fmpz::{Fmpz, FmpzFactor};
use flint::fmpq::Fmpq;
use flint::ulong_extras;
use gmp::mpz::Mpz;

fn norm_realquad(u: i64, v: i64, m: u64) -> i64 {
    if is_1mod4!(m) {
        (u * u - (m as i64) * v * v) >> 1
    } else {
        u * u - m as i64 * v * v
    }
}

fn content(u: i64, v: usize, m: u64) -> u64 {
    let d = ulong_extras::gcd(u.abs() as u64, v as u64);
    if is_1mod4!(m) && !is_even!(u - v as i64) {
        d >> 1
    } else {
        d
    }
}

fn disc(m: u64) -> u64 {
    if is_1mod4!(m) { m } else { 4 * m }
}

/// This assumes narrow class number of Q(sqrt{m}) is one.
pub fn eisenstein_series_from_lvals(
    k: u64,
    m: u64,
    l_val_num: &Fmpz,
    l_val_den: &Fmpz,
    prec: usize,
) -> HmfGen<Fmpq> {
    let mut res = HmfGen::<Fmpq>::new(m, prec);
    res.fcvec.fc_ref_mut(0, 0, 0).set_ui_g(1);
    let mut tmp_z = Fmpz::new();
    let mut fac = FmpzFactor::new();
    let mut tmp_mpz = Mpz::new();
    tmp_mpz.set_si(m as i64);
    let dsc = disc(m);
    v_u_bd_iter_non_const!((res.m, res.u_bds, v, u, bd) {
        let norm = norm_realquad(u, v as i64, m).abs();
        let d = content(u, v, m) as i64;
        tmp_z.set_si(norm);
        fac.factor_mut(&tmp_z);
        for &(ref p_z, n) in fac.to_vec().iter() {
            let p = p_z.to_slong().unwrap();
            let chi = Mpz::kronecker_si(&tmp_mpz, p);
        }
    });
    res
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn it_works() {
//         let a = Mpz::from_si(40);
//         for &p in prime_sieve(100).iter() {
//             println!("{}: {}", p, Mpz::kronecker_si(&a, p as i64));
//         }
//     }
// }
