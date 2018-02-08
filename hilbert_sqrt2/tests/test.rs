extern crate hilbert_sqrt2;
extern crate hilbert_qexp;

mod parallel {
    use hilbert_sqrt2::parallel_wt::*;

    #[test]
    fn test_s4() {
        let s4 = s4_form(5);
        assert_eq!(s4.weight, Some((4, 4)));
        println!("{}", s4);
    }

    #[test]
    fn test_s6() {
        let s6 = s6_form(5);
        assert_eq!(s6.weight, Some((6, 6)));
        println!("{}", s6);
    }

    #[test]
    fn test_s5() {
        let s5 = s5_form(5);
        assert_eq!(s5.weight, Some((5, 5)));
        println!("{}", s5);
    }

    #[test]
    fn test_s9() {
        let s9 = s9_form(5);
        assert_eq!(s9.weight, Some((9, 9)));
        println!("{}", s9);
    }
}

mod structure {
    use hilbert_sqrt2::structure::*;
    use hilbert_qexp::bignum::RealQuadElement;
    use hilbert_sqrt2::mixed_wt::*;

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
