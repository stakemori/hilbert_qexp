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


fn proj_f9_part<T>(f: &HmfGen<T>, f9: &HmfGen<T>) -> HmfGen<T>
where
    T: BigNumber + Clone + MulAssign<u64> + ShrAssign<usize>,
    for<'a> T: SubAssign<&'a T>,
{
    let mut tmp = HmfGen::new(f.m, f.prec);
    star_op(&mut tmp, f);
    let g = &tmp + f;
    let mut tmp1 = HmfGen::new(f.m, f.prec);
    div_mut(&mut tmp1, &g, f9);
    tmp1 >>= 1;
    tmp1
}


pub fn bracket_inner_prod<T>(f: &HmfGen<T>, g: &HmfGen<T>) -> HmfGen<T>
where
    T: BigNumber + Clone + MulAssign<u64> + ShrAssign<usize>,
    for<'a> T: AddAssign<&'a T>,
    for<'a> T: SubAssign<&'a T>,
    for<'a> T: From<&'a Fmpq>,
{
    let prec = max(f.prec, g.prec);
    let mut tmp = HmfGen::new(f.m, prec);
    star_op(&mut tmp, g);
    tmp *= f;
    let s9 = s9_form(prec);
    proj_f9_part(&tmp, &Into::<HmfGen<T>>::into(&s9))
}
