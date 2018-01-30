use elements::HmfGen;
use bignum::{BigNumber, RealQuadElement};
use flint::fmpq::Fmpq;
use flint::fmpz::Fmpz;
use std::ops::{MulAssign, AddAssign};

fn diff_mut<T>(res: &mut HmfGen<T>, expt: (usize, usize), f: &HmfGen<T>)
where
    T: BigNumber + From<(i64, u64)> + Clone,
{
    res.set(f);
    let mut tmp = T::R::default();
    let mut tmp_t = T::new_g();
    let mut tmp_t1 = T::new_g();
    v_u_bd_iter!(
        (f.m, f.u_bds, v, u, bd)
        {
            tmp_t1.set_ui_g(1);
            if expt.0 > 0 {
                let a: T = From::from((u, v as u64));
                tmp_t.pow_mut(&a, expt.0);
                tmp_t1.mul_assign_g(&tmp_t, &mut tmp);
            }
            if expt.1 > 0 {
                let a: T = From::from((-u, v as u64));
                tmp_t.pow_mut(&a, expt.1);
                tmp_t1.mul_assign_g(&tmp_t, &mut tmp);
            }
            res.fcvec.fc_ref_mut(v, u, bd).mul_assign_g(&tmp_t1, &mut tmp);
        }
    );
}

fn diff_mul<T>(
    res: &mut HmfGen<T>,
    expt0: (usize, usize),
    expt1: (usize, usize),
    f: &HmfGen<T>,
    g: &HmfGen<T>,
) where
    T: BigNumber + From<(i64, u64)> + Clone,
    for<'a> T: AddAssign<&'a T>,
{
    let mut tmp = HmfGen::<T>::new(f.m, f.prec);
    diff_mut(res, expt0, f);
    diff_mut(&mut tmp, expt1, g);
    *res *= &tmp;
}

#[derive(Debug)]
pub struct NotHhmError {}

pub fn rankin_cohen<T>(n: usize, f: &HmfGen<T>, g: &HmfGen<T>) -> Result<HmfGen<T>, NotHhmError>
where
    T: BigNumber + From<(i64, u64)> + Clone + RealQuadElement<Fmpq>,
    for<'a> T: MulAssign<&'a Fmpq>,
    for<'a> T: AddAssign<&'a T>,
{
    assert_eq!(f.prec, g.prec);
    if !f.weight.is_none() && !g.weight.is_none() {
        let mut res = HmfGen::<T>::new(f.m, f.prec);
        let mut tmp = HmfGen::<T>::new(f.m, f.prec);

        let mut tmp_z = Fmpz::new();
        let mut tmp_z1 = Fmpz::new();
        let mut sgn = if is_even!(n) { 1 } else { -1 } as i64;

        let (k1, k2) = f.weight.unwrap();
        let (l1, l2) = g.weight.unwrap();

        for i in 0..(n + 1) {
            tmp_z.bi_uiui_mut((n + k2 - 1) as u64, (n - i) as u64);
            tmp_z1.bi_uiui_mut((n + l2 - 1) as u64, i as u64);
            tmp_z *= &tmp_z1;
            tmp_z *= sgn;
            diff_mul(&mut tmp, (0, i), (0, n - i), &f, &g);
            tmp *= &Into::<Fmpq>::into(&tmp_z);
            res += &tmp;
            sgn *= -1;
        }
        res.weight = Some((k1 + l1, k2 + l2 + 2 * n));
        Ok(res)
    } else {
        Err(NotHhmError {})
    }
}
