use hilbert_qexp::diff_op::rankin_cohen;
use hilbert_qexp::elements::HmfGen;
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
