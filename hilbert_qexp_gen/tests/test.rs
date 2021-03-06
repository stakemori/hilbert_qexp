extern crate flint;
extern crate hilbert_qexp;

mod sqrt5 {
    use hilbert_qexp::eisenstein::eisenstein_series_from_lvals;
    use flint::fmpq::Fmpq;
    use hilbert_qexp::elements::{square_root_mut, HmfGen};
    use hilbert_qexp::misc::PowGen;
    use hilbert_qexp::diff_op::rankin_cohen;
    use hilbert_qexp::bignum::Sqrt5Q;

    fn eisenstein_sqrt5(k: usize, prec: usize) -> HmfGen<Fmpq> {
        assert!(k == 2 || k == 6 || k == 10);
        if k == 2 {
            eisenstein_series_from_lvals(2, 5, &From::from(120_i64), &From::from(1_i64), prec)
        } else if k == 6 {
            let mut res =
                eisenstein_series_from_lvals(6, 5, &From::from(2520), &From::from(67), prec);
            res *= 67;
            res
        } else {
            let mut res =
                eisenstein_series_from_lvals(10, 5, &From::from(6600), &From::from(412751), prec);
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

    #[test]
    fn test_rankin_cohen() {
        let prec = 5;
        let g2: HmfGen<Sqrt5Q> = From::from(&eisenstein_sqrt5(2, prec));
        let mut f = rankin_cohen(2, &g2, &g2).unwrap();
        let a: Sqrt5Q = From::from((4320, -1440_i64));
        let b: Sqrt5Q = From::from((3, -1_i64));
        f /= &a;
        f *= &b;
        for x in f.diagonal_restriction().iter() {
            println!("{}", x);
        }
        println!("{}", f);
    }
}

mod eisen_sqrt2 {
    use hilbert_qexp::eisenstein::eisenstein_series_from_lvals;
    use flint::fmpq::Fmpq;
    use hilbert_qexp::elements::{relations_over_q, HmfGen};
    // use hilbert_qexp::misc::PowGen;

    fn eisenstein_sqrt2(k: u64, prec: usize) -> HmfGen<Fmpq> {
        assert!(k == 2 || k == 4);
        match k {
            2 => eisenstein_series_from_lvals(k, 2, &From::from(48), &From::from(1), prec),
            _ => eisenstein_series_from_lvals(k, 2, &From::from(480), &From::from(11), prec),
        }
    }

    fn f4(prec: usize) -> HmfGen<Fmpq> {
        let e2 = eisenstein_sqrt2(2, prec);
        let mut res = eisenstein_sqrt2(4, prec);
        res -= &(&e2 * &e2);
        res
    }

    #[test]
    fn diag_res() {
        let v = eisenstein_sqrt2(2, 5).diagonal_restriction();
        assert_eq!(
            v,
            vec![1, 240, 2160, 6720, 17520, 30240]
                .iter()
                .map(|&i| Into::<Fmpq>::into((i, 1)))
                .collect::<Vec<_>>()
        );

        let v = eisenstein_sqrt2(4, 5).diagonal_restriction();
        assert_eq!(
            v,
            vec![1, 480, 61920, 1050240, 7926240, 37500480]
                .iter()
                .map(|&i| Into::<Fmpq>::into((i, 1)))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_wt4_cusp_form() {
        let mut f = f4(5);
        f /= &From::from((-576, 11));
        assert!(f.diagonal_restriction().iter().all(|v| v.is_zero()));
    }

    #[test]
    fn test_relations() {
        let prec = 5;
        let f = eisenstein_sqrt2(2, prec);
        let forms = vec![f.clone(), &f * 2];
        let v = relations_over_q(&forms);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].len(), 2);
        assert_eq!(v[0][0], &v[0][1] * (-2));
    }
}
