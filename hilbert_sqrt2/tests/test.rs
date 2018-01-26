extern crate hilbert_sqrt2;

mod parallel {
    use hilbert_sqrt2::parallel_wt::*;

    #[test]
    fn test_s4() {
        let s4 = s4_form(5);
        println!("{}", s4);
    }

    #[test]
    fn test_s6() {
        let s6 = s6_form(5);
        println!("{}", s6);
    }

    #[test]
    fn test_s9() {
        let s9 = s9_form(5);
        println!("{}", s9);
    }
}
