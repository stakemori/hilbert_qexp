use hilbert_qexp::diff_op::rankin_cohen;
use hilbert_qexp::elements::{HmfGen, relations_over_q, div_mut};
use hilbert_qexp::bignum::Sqrt2Q;
use parallel_wt::*;
use mixed_wt::*;
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

fn div_by_s5(f: &HmfGen<Fmpq>) -> HmfGen<Fmpq> {
    let s5 = s5_form(f.prec);
    let mut res = HmfGen::new(2, f.prec);
    div_mut(&mut res, f, &s5);
    res
}


#[cfg(test)]
mod tests {
    use super::*;
    use hilbert_qexp::bignum::RealQuadElement;

    #[test]
    fn br_a1_0() {
        let prec = 6;
        let forms = three_forms_a1_0(prec);
        let br01 = bracket_inner_prod(&forms[0], &forms[1]);
        assert!(br01.rt_part().is_zero());
        let f01 = div_by_s5(&br01.ir_part());
        // println!("{}", &f - &(&s2_form(prec -2) * 4));
        let pl01 = r_elt_as_pol_over_q(&f01);
        println!("{:?}", pl01);

        let br02 = bracket_inner_prod(&forms[0], &forms[2]);
        let br12 = bracket_inner_prod(&forms[1], &forms[2]);
        assert!(br02.rt_part().is_zero());
        assert!(br12.rt_part().is_zero());

        let f02 = div_by_s5(&br02.ir_part());
        let f12 = div_by_s5(&br12.ir_part());
        let pl02 = r_elt_as_pol_over_q(&f02);
        let pl12 = r_elt_as_pol_over_q(&f12);
        println!("{:?}", pl02);
        println!("{:?}", pl12);
    }
}
