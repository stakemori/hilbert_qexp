extern crate libc;
extern crate flint;

pub mod bignum;
pub mod elements;
mod misc;
mod fcvec;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
