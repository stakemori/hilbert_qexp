#[macro_use]
extern crate hilbert_qexp;
extern crate flint;
extern crate gmp;

pub mod parallel_wt;
pub mod mixed_wt;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
