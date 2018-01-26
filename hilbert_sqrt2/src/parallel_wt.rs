use hilbert_qexp::eisenstein::eisenstein_series_from_lvals;
use hilbert_qexp::elements::{HmfGen, square_root_mut};
use flint::fmpq::Fmpq;
use flint::fmpz::Fmpz;

pub fn eisensten_series(k: u64, prec: usize) -> HmfGen<Fmpq> {
    assert!(2 <= k && k <= 10);
    let l_vals = [
        ("48", "1"),
        ("480", "11"),
        ("1008", "361"),
        ("960", "24611"),
        ("528", "2873041"),
    ];
    let (num_s, den_s) = l_vals[((k - 2) >> 1) as usize];
    let num = &Fmpz::from_str(&num_s, 10).unwrap();
    let den = &Fmpz::from_str(&den_s, 10).unwrap();
    eisenstein_series_from_lvals(k, 2, num, den, prec)
}

pub fn s4_form(prec: usize) -> HmfGen<Fmpq> {
    let g2 = eisensten_series(2, prec);
    let g4 = eisensten_series(4, prec);
    let mut res = &g2 * &g2;
    res -= &g4;
    res *= 11;
    res /= &Into::<Fmpq>::into((192 * 3, 1));
    res
}

pub fn s6_form(prec: usize) -> HmfGen<Fmpq> {
    let g2 = eisensten_series(2, prec);
    let g4 = eisensten_series(4, prec);
    let g6 = eisensten_series(6, prec);
    let mut res = g2.pow(3);
    res *= -25 * 49;
    res += &(&(&g2 * &g4) * (3 * 11 * 59));
    res -= &(&g6 * (2 * 19 * 19));
    res /= &From::from(449280);
    res
}

fn s5_squared(prec: usize) -> HmfGen<Fmpq> {
    let s4 = s4_form(prec);
    let s6 = s6_form(prec);
    let g2 = eisensten_series(2, prec);
    (&s6 * 4 + &s4 * &g2) * &s4
}

pub fn s5_form(prec: usize) -> HmfGen<Fmpq> {
    let mut res = HmfGen::<Fmpq>::new(2, prec + 1);
    let bd = res.u_bds.vec[1] as i64;
    res.fcvec.fc_ref_mut(1, 1, bd).set_si(1, 1);
    res.fcvec.fc_ref_mut(1, -1, bd).set_si(-1, 1);
    square_root_mut(&mut res, 1, &s5_squared(prec + 1));
    res
}

fn s9_squared(prec: usize) -> HmfGen<Fmpq> {
    let g2 = eisensten_series(2, prec);
    let s4 = s4_form(prec);
    let s6 = s6_form(prec);
    let mut res = &g2.pow(3) * &s6;
    res += &(&(&g2.pow(2) * &s4.pow(2)) * 4);
    res -= &(&(&g2 * &(&s4 * &s6)) * (32 * 9));
    res -= &(&s4.pow(3) * 1024);
    res -= &(&s6.pow(2) * (64 * 27));
    res * &s6
}

pub fn s9_form(prec: usize) -> HmfGen<Fmpq> {
    let f = s9_squared(prec + 1);
    let mut res = HmfGen::<Fmpq>::new(2, prec + 1);
    let bd = res.u_bds.vec[1] as i64;
    res.fcvec.fc_ref_mut(1, 0, bd).set_si(1, 1);
    square_root_mut(&mut res, 1, &f);
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_s5_squared() {
        let prec = 5;
        let f = s5_squared(prec);
        let s5 = s5_form(prec);
        assert_eq!(f, &s5 * &s5);
    }

    #[test]
    fn test_s9_squared() {
        let prec = 5;
        let f = s9_squared(prec);
        println!("{}", f);
        assert_eq!(
            f.diagonal_restriction(),
            vec![0, 0, 1, -1056, 270216, 4819328]
                .into_iter()
                .map(|x| Into::<Fmpq>::into((x, 1)))
                .collect::<Vec<_>>()
        );
    }
}
