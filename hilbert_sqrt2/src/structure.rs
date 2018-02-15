use hilbert_qexp::diff_op::rankin_cohen;
use hilbert_qexp::elements::{HmfGen, relations_over_q, div_mut};
use hilbert_qexp::bignum::Sqrt2Q;
use hilbert_qexp::bignum::RealQuadElement;
use parallel_wt::*;
use mixed_wt::*;
use flint::fmpq::Fmpq;
use std::fs::File;
use serde;
use serde_pickle;
use std::io::Write;

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
    for c in (0..(c_max + 1)).rev() {
        let b_max = (k - 6 * c) / 4;
        for b in (0..(b_max + 1)).rev() {
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
    is_odd: bool,
    prec: usize,
    len: usize,
) -> Vec<(HmfGen<Sqrt2Q>, Tuple3, Tuple3)> {
    if is_odd {
        mixed_weight_forms_odd(df, prec, len)
    } else {
        mixed_weight_forms_even(df, prec, len)
    }
}

pub fn mixed_weight_forms_even(
    df: usize,
    prec: usize,
    len: usize,
) -> Vec<(HmfGen<Sqrt2Q>, Tuple3, Tuple3)> {
    let mut num = 0;
    let mut res = Vec::with_capacity(len);
    for (i, m) in (2..).flat_map(monoms_of_s2_s4_s6).enumerate() {
        for n in (2..).flat_map(monoms_of_s2_s4_s6).take(i + 1) {
            if num >= len {
                return res;
            }
            let f_m = m.to_form(prec);
            let f_n = n.to_form(prec);
            if let Ok(f) = rankin_cohen(df, &From::from(&f_m), &From::from(&f_n)) {
                if !f.is_zero() {
                    num += 1;
                    res.push((f, m.idx, n.idx));
                }
            } else {
                panic!();
            }
        }
    }
    res
}

pub fn mixed_weight_forms_odd(
    df: usize,
    prec: usize,
    len: usize,
) -> Vec<(HmfGen<Sqrt2Q>, Tuple3, Tuple3)> {
    let mut num = 0;
    let mut res = Vec::with_capacity(len);
    let s5 = s5_form(prec);
    for n in (2..).flat_map(monoms_of_s2_s4_s6) {
        if num >= len {
            return res;
        }
        let f_n = n.to_form(prec);

        if let Ok(f) = rankin_cohen(df, &From::from(&s5.clone()), &From::from(&f_n)) {
            if !f.is_zero() {
                num += 1;
                res.push((f, (0, 0, 0), n.idx));
            }
        } else {
            panic!();
        }
    }
    res
}

pub type ParaWtPolyQ = Vec<(MonomFormal, Fmpq)>;

/// If the bracket of f an g is of odd weight, this returns a modular form
/// divided by s5.
pub fn bracket_inner_prod_as_pol_over_q(
    f: &HmfGen<Sqrt2Q>,
    g: &HmfGen<Sqrt2Q>,
) -> Option<ParaWtPolyQ> {
    let h = bracket_inner_prod(f, g);
    if h.is_zero() {
        return Some(vec![]);
    }
    if !h.rt_part().is_zero() {
        None
    } else {
        let h_ir = h.ir_part();
        if is_even!(h_ir.weight.unwrap().0) {
            r_elt_as_pol_over_q(&h_ir)
        } else {
            let f = div_by_s5(&h_ir);
            r_elt_as_pol_over_q(&f)
        }
    }
}

pub fn brackets(forms: &[HmfGen<Sqrt2Q>]) -> Vec<ParaWtPolyQ> {
    forms
        .iter()
        .enumerate()
        .flat_map(|(i, f)| {
            forms
                .iter()
                .skip(i + 1)
                .map(|g| bracket_inner_prod_as_pol_over_q(f, g).unwrap())
                .collect::<Vec<_>>()
        })
        .collect()
}

pub fn save_as_pickle<T>(a: T, f: &mut File)
where
    T: serde::Serialize,
{
    let v = serde_pickle::to_vec(&a, false).unwrap();
    f.write(&v).unwrap();
}

pub fn save_brackets_for_candidates<'a, I>(vals_iter: I, len: usize)
where
    I: Iterator<Item = &'a (usize, u32)>,
{
    for i_parity in vals_iter {
        let i = i_parity.0;
        let parity = i_parity.1;
        println!("{}, {}", i, parity);
        let prec = (2 * i + 6) / 5 + 10;
        let forms_w_monoms = mixed_weight_forms(i, !is_even!(parity), prec, len);

        {
            let monoms = forms_w_monoms
                .iter()
                .map(|f_t| (f_t.1, f_t.2))
                .collect::<Vec<_>>();
            let ref mut monms_file =
                File::create(format!("./data/str{}_{}_monoms.sobj", i, parity)).unwrap();

            save_as_pickle(&monoms, monms_file);
        }

        {
            let weights = forms_w_monoms
                .iter()
                .map(|f| f.0.weight.unwrap().0)
                .collect::<Vec<_>>();
            println!("{:?}", weights);
            let ref mut weights_file =
                File::create(format!("./data/str{}_{}_weights.sobj", i, parity)).unwrap();
            save_as_pickle(&weights, weights_file);
        }

        {
            let forms = forms_w_monoms
                .clone()
                .into_iter()
                .map(|f| f.0)
                .collect::<Vec<_>>();
            let brs = brackets(&forms);
            let brs = brs.iter()
                .map(|br| {
                    br.iter()
                        .map(|x| (x.0.idx, x.1.clone()))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let ref mut br_file = File::create(format!("./data/str{}_{}_brs.sobj", i, parity))
                .unwrap();
            save_as_pickle(&brs, br_file);
        }
    }
}
