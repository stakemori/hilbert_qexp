use hilbert_qexp::diff_op::rankin_cohen;
use hilbert_qexp::elements::{HmfGen, relations_over_q, div_mut};
use hilbert_qexp::bignum::Sqrt2Q;
use parallel_wt::*;
use flint::fmpq::Fmpq;

/// Corresponds to s2^a * s4^b * s6^c where (a, b, c) = idx.
#[derive(Clone, Debug)]
pub struct MonomFormal {
    pub idx: (usize, usize, usize),
}

fn monom_s2_s4_s6(prec: usize, expts: (usize, usize, usize)) -> HmfGen<Fmpq> {
    let mut tmp = HmfGen::new(2, prec);
    let mut res = HmfGen::new(2, prec);
    let s2 = s2_form(prec);
    let s4 = s4_form(prec);
    let s6 = s6_form(prec);
    let (e2, e4, e6) = expts;
    res.pow_mut(&s2, e2);
    tmp.pow_mut(&s6, e6);
    if e6 > 0 {
        res *= &tmp;
    }
    tmp.pow_mut(&s4, e4);
    if e4 > 0 {
        res *= &tmp;
    }
    res
}

impl MonomFormal {
    pub fn to_form(&self, prec: usize) -> HmfGen<Fmpq> {
        monom_s2_s4_s6(prec, self.idx)
    }
}

/// Return a vector of (a, b, c) s.t. 2*a + 4*b + 6*c = k.
pub fn tpls_of_wt(k: usize) -> Vec<(usize, usize, usize)> {
    let mut res = Vec::new();
    let c_max = k / 6;
    for c in 0..(c_max + 1) {
        let b_max = (k - 6 * c) / 4;
        for b in 0..(b_max + 1) {
            let rem = k - (6 * c + 4 * b);
            if is_even!(rem) {
                res.push((rem >> 1, b, c));
            }
        }
    }
    res
}

pub fn monoms_of_s2_s4_s6(k: usize) -> Vec<MonomFormal> {
    tpls_of_wt(k)
        .into_iter()
        .map(|x| MonomFormal { idx: x })
        .collect()
}

pub fn r_elt_as_pol_over_q(f: &HmfGen<Fmpq>) -> Option<Vec<(MonomFormal, Fmpq)>> {
    let prec = f.prec;
    let monoms = monoms_of_s2_s4_s6(f.weight.unwrap().0);
    let mut forms: Vec<_> = monoms.iter().map(|x| x.to_form(prec)).collect();
    forms.insert(0, f.clone());
    let rels = relations_over_q(&forms);
    if rels.len() == 1 && !rels[0][0].is_zero() {
        let cfs: Vec<_> = rels[0]
            .iter()
            .skip(1)
            .map(|x| -&(x / &rels[0][0]))
            .collect();
        Some(monoms.into_iter().zip(cfs.into_iter()).collect())
    } else {
        None
    }
}


pub fn three_forms_a1_0(prec: usize) -> Vec<HmfGen<Sqrt2Q>> {
    let s2 = Into::<HmfGen<Sqrt2Q>>::into(&s2_form(prec));
    let s4 = Into::<HmfGen<Sqrt2Q>>::into(&s4_form(prec));
    let s6 = Into::<HmfGen<Sqrt2Q>>::into(&s6_form(prec));
    vec![
        rankin_cohen(1, &s2, &s4).unwrap(),
        rankin_cohen(1, &s2, &s6).unwrap(),
        rankin_cohen(1, &s4, &s6).unwrap(),
    ]
}

pub fn three_forms_a1_1(prec: usize) -> Vec<HmfGen<Sqrt2Q>> {
    let s2 = Into::<HmfGen<Sqrt2Q>>::into(&s2_form(prec));
    let s4 = Into::<HmfGen<Sqrt2Q>>::into(&s4_form(prec));
    let s6 = Into::<HmfGen<Sqrt2Q>>::into(&s6_form(prec));
    let s5 = Into::<HmfGen<Sqrt2Q>>::into(&s5_form(prec));
    vec![
        rankin_cohen(1, &s2, &s5).unwrap(),
        rankin_cohen(1, &s4, &s5).unwrap(),
        rankin_cohen(1, &s6, &s5).unwrap(),
    ]
}

pub fn div_by_s5(f: &HmfGen<Fmpq>) -> HmfGen<Fmpq> {
    let s5 = s5_form(f.prec);
    let mut res = HmfGen::new(2, f.prec);
    div_mut(&mut res, f, &s5);
    res
}

type Tuple3 = (usize, usize, usize);

pub fn mixed_weight_forms(
    df: usize,
    prec: usize,
    len: usize,
) -> Vec<(HmfGen<Sqrt2Q>, Tuple3, Tuple3)> {
    let mut num = 0;
    let mut res = Vec::with_capacity(len);
    for (i, m) in (2..).flat_map(monoms_of_s2_s4_s6).enumerate() {
        for n in (2..).flat_map(monoms_of_s2_s4_s6).take(if is_even!(df) {
            i + 1
        } else {
            i
        })
        {
            if num >= len {
                return res;
            }
            let f_m = m.to_form(prec);
            let f_n = n.to_form(prec);
            let f = rankin_cohen(df, &From::from(&f_m), &From::from(&f_n)).unwrap();
            if !f.is_zero() {
                num += 1;
                res.push((f, m.idx, n.idx));
            }
        }
    }
    res
}
