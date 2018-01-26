extern crate hilbert_qexp;
extern crate flint;

mod square_root {
    use hilbert_qexp::eisenstein::eisenstein_series_from_lvals;
    use flint::fmpq::Fmpq;
    use hilbert_qexp::elements::{HmfGen, square_root_mut};
    use hilbert_qexp::misc::PowGen;

    fn eisenstein_sqrt5(k: usize, prec: usize) -> HmfGen<Fmpq> {
        assert!(k == 2 || k == 6 || k == 10);
        if k == 2 {
            eisenstein_series_from_lvals(
                2,
                5,
                &From::from(120_i64),
                &From::from(1_i64),
                prec,
            )
        } else if k == 6 {
            let mut res = eisenstein_series_from_lvals(
                6,
                5,
                &From::from(2520),
                &From::from(67),
                prec,
            );
            res *= 67;
            res
        } else {
            let mut res = eisenstein_series_from_lvals(
                10,
                5,
                &From::from(6600),
                &From::from(412751),
                prec,
            );
            res *= 412751;
            res
        }
    }

    fn g5_squared(prec: usize) -> HmfGen<Fmpq> {
        let e2 = eisenstein_sqrt5(2, prec);
        let e6 = eisenstein_sqrt5(6, prec);
        let mut res = eisenstein_sqrt5(10, prec);
        let mut tmp = HmfGen::new(5, prec);
        tmp.pow_mut(&e2, 2);
        let mut f10_2 = &tmp * &e6;
        tmp.square();
        tmp *= &e2;
        let mut f10_1 = tmp.clone();
        f10_2 *= &Into::<Fmpq>::into((-11465, 1));
        f10_1 *= &Into::<Fmpq>::into((355404, 1));
        res += &f10_2;
        res += &f10_1;
        res /= &Into::<Fmpq>::into((5443200000, 1));
        res
    }

    fn g5_normalized(prec: usize) -> HmfGen<Fmpq> {
        let mut res = HmfGen::<Fmpq>::new(5, prec + 1);
        let bd = res.u_bds.vec[1] as i64;
        res.fcvec.fc_ref_mut(1, 1, bd).set_si(1, 1);
        res.fcvec.fc_ref_mut(1, -1, bd).set_si(-1, 1);
        square_root_mut(&mut res, 1, &g5_squared(prec + 1));
        res
    }

    #[test]
    fn test_square_root() {
        let prec = 5;
        let g5 = g5_normalized(prec);
        let g10 = g5_squared(prec);
        assert_eq!(&g5 * &g5, g10);
    }
}
