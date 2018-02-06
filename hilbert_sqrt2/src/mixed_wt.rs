use hilbert_qexp::elements::{HmfGen, div_mut};
use hilbert_qexp::bignum::BigNumber;
use parallel_wt::*;
use flint::fmpq::Fmpq;

use std::cmp::max;
use std::ops::*;

pub fn star_op<T>(res: &mut HmfGen<T>, f: &HmfGen<T>)
where
    T: BigNumber,
{
    let (k1, k2) = f.weight.unwrap();
    v_u_bd_iter!((f.m, f.u_bds, v, u, bd) {
        res.fcvec.fc_ref_mut(v, u, bd).set_g(f.fcvec.fc_ref(v, -u, bd));
    });
    if !is_even!((k1 + k2) >> 1) {
        res.negate();
    }
    res.weight = Some((k2, k1));
}


fn proj_s9_part<T>(f: &HmfGen<T>, s9: &HmfGen<T>) -> HmfGen<T>
where
    T: BigNumber + Clone + ShrAssign<u64>,
    for<'a> T: SubAssign<&'a T>,
{
    let mut tmp = HmfGen::new(f.m, f.prec);
    star_op(&mut tmp, f);
    let g = f - &tmp;
    let mut tmp1 = HmfGen::new(f.m, f.prec);
    div_mut(&mut tmp1, &g, s9);
    tmp1 >>= 1_u64;
    tmp1
}


pub fn bracket_inner_prod<T>(f: &HmfGen<T>, g: &HmfGen<T>) -> HmfGen<T>
where
    T: BigNumber + Clone + ShrAssign<u64>,
    for<'a> T: AddAssign<&'a T>,
    for<'a> T: SubAssign<&'a T>,
    for<'a> T: From<&'a Fmpq>,
{
    let prec = max(f.prec, g.prec);
    let mut tmp = HmfGen::new(f.m, prec);
    star_op(&mut tmp, g);
    tmp *= f;
    let s9 = s9_form(prec);
    proj_s9_part(&tmp, &Into::<HmfGen<T>>::into(&s9))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_proj_s9() {
        let prec = 5;
        let s4 = s4_form(prec);
        let s6 = s6_form(prec);
        let mut s5 = s5_form(prec);
        let s9 = s9_form(prec);
        let f = &(&(s4.pow(2)) * &s6) + &(&s5 * &s9);
        let g = proj_s9_part(&f, &s9);
        s5.decrease_prec(prec - 1);
        assert_eq!(g, s5);
    }
}
