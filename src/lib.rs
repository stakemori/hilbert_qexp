extern crate libc;
extern crate flint;

#[macro_use]
pub mod elements;
pub mod bignum;
mod misc;
mod fcvec;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
