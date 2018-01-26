extern crate libc;
extern crate flint;
extern crate gmp;

#[macro_use]
pub mod elements;
pub mod bignum;
pub mod eisenstein;
pub mod misc;
pub mod diff_op;
mod fcvec;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
