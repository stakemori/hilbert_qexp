extern crate hilbert_sqrt2;

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
