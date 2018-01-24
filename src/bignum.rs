use libc::{c_ulong, c_long};
use flint::fmpz::Fmpz;
use std::fmt;
use std::ops::{AddAssign, SubAssign, ShlAssign, ShrAssign, MulAssign};
use std;
use flint::fmpz_poly::FmpzPoly;
use flint::fmpq::Fmpq;

pub trait RealQuadElement<S> {
    fn rt_part(&self) -> S;
    fn ir_part(&self) -> S;
}

pub trait BigNumber {
    fn is_zero_g(&self) -> bool;
    fn from_ui_g(x: c_ulong) -> Self;
    fn from_si_g(x: c_long) -> Self;
    fn set_ui_g(&mut self, x: c_ulong);
    fn set_si_g(&mut self, x: c_long);
    fn set_g(&mut self, other: &Self);
    fn add_mut_g(&mut self, x: &Self, y: &Self);
    fn mul_mut_g(&mut self, x: &Self, y: &Self);
    fn sub_mut_g(&mut self, x: &Self, y: &Self);
    fn is_multiple_of_g(&self, x: &Self, tmpelt: &mut Self, tmp: &mut Fmpz) -> bool;
    fn new_g() -> Self;
    fn negate_g(&mut self);
    fn set_divexact_g(&mut self, x: &Self, tmp: &mut Fmpz);
    fn mul_assign_g(&mut self, other: &Self, tmp: &mut Fmpz);
    fn addmul_mut_g(&mut self, x: &Self, y: &Self, tmp: &mut Fmpz);
    fn submul_mut_g(&mut self, x: &Self, y: &Self, tmp: &mut Fmpz);
    fn square_g(&mut self, tmp_elt: &mut Self, tmp: &mut Fmpz);
    fn pow_mut(&mut self, f: &Self, a: usize)
    where
        Self: std::marker::Sized + Clone,
    {
        self.set_ui_g(1);
        let s = format!("{:b}", a);
        let bts = s.into_bytes();
        let strs: Vec<char> = bts.iter().rev().map(|&i| i as char).collect();
        let mut n = (*f).clone();
        let tmp_elt = &mut Self::new_g();
        let tmp = &mut Fmpz::new();
        for &c in strs.iter() {
            if c == '0' {
                n.square_g(tmp_elt, tmp);
            } else if c == '1' {
                self.mul_assign_g(&n, tmp);
                n.square_g(tmp_elt, tmp);
            }
        }

    }
}

impl BigNumber for Fmpq {
    fn is_zero_g(&self) -> bool {
        self.is_zero()
    }

    fn from_ui_g(x: c_ulong) -> Self {
        debug_assert!(x < ::std::i64::MAX as u64);
        From::from((x as c_long, 1))
    }

    fn from_si_g(x: c_long) -> Self {
        From::from((x, 1))
    }

    fn set_ui_g(&mut self, x: c_ulong) {
        self.set_ui(x, 1);
    }

    fn set_si_g(&mut self, x: c_long) {
        self.set_si(x, 1);
    }

    fn set_g(&mut self, other: &Self) {
        self.set(other);
    }

    fn add_mut_g(&mut self, x: &Self, y: &Self) {
        self.add_mut(x, y);
    }

    fn sub_mut_g(&mut self, x: &Self, y: &Self) {
        self.sub_mut(x, y);
    }

    fn mul_mut_g(&mut self, x: &Self, y: &Self) {
        self.mul_mut(x, y);
    }

    fn is_multiple_of_g(&self, x: &Self, _tmpelt: &mut Self, _tmp: &mut Fmpz) -> bool {
        if x.is_zero() { self.is_zero() } else { true }
    }

    fn new_g() -> Self {
        Fmpq::new()
    }

    fn set_divexact_g(&mut self, x: &Self, _tmp: &mut Fmpz) {
        debug_assert!(!x.is_zero());
        *self /= x;
    }

    fn mul_assign_g(&mut self, other: &Self, _tmp: &mut Fmpz) {
        *self *= other;
    }

    fn addmul_mut_g(&mut self, x: &Self, y: &Self, _tmp: &mut Fmpz) {
        self.addmul_mut(x, y);
    }

    fn submul_mut_g(&mut self, x: &Self, y: &Self, _tmp: &mut Fmpz) {
        self.submul_mut(x, y);
    }

    fn square_g(&mut self, _tmp_elt: &mut Self, _tmp: &mut Fmpz) {
        self.set_pow_si(2);
    }

    fn negate_g(&mut self) {
        self.negate();
    }
}

impl BigNumber for Fmpz {
    fn is_zero_g(&self) -> bool {
        self.is_zero()
    }

    fn from_ui_g(x: c_ulong) -> Self {
        Fmpz::from_ui(x)
    }

    fn from_si_g(x: c_long) -> Self {
        Fmpz::from_si(x)
    }

    fn set_ui_g(&mut self, x: c_ulong) {
        self.set_ui(x);
    }

    fn set_si_g(&mut self, x: c_long) {
        self.set_si(x);
    }

    fn set_g(&mut self, other: &Fmpz) {
        self.set(other)
    }

    fn add_mut_g(&mut self, x: &Fmpz, y: &Fmpz) {
        self.add_mut(x, y);
    }

    fn mul_mut_g(&mut self, x: &Fmpz, y: &Fmpz) {
        self.mul_mut(x, y);
    }

    fn sub_mut_g(&mut self, x: &Fmpz, y: &Fmpz) {
        self.sub_mut(x, y);
    }

    fn is_multiple_of_g(&self, x: &Fmpz, _tmpelt: &mut Self, _tmp: &mut Fmpz) -> bool {
        self.is_divisible(x)
    }

    fn new_g() -> Fmpz {
        Fmpz::new()
    }

    fn negate_g(&mut self) {
        self.negate()
    }

    fn set_divexact_g(&mut self, x: &Fmpz, _tmp: &mut Fmpz) {
        self.set_divexact(x);
    }

    fn mul_assign_g(&mut self, other: &Fmpz, _tmp: &mut Fmpz) {
        *self *= other;
    }

    fn addmul_mut_g(&mut self, x: &Fmpz, y: &Fmpz, _tmp: &mut Fmpz) {
        self.addmul_mut(x, y);
    }

    fn submul_mut_g(&mut self, x: &Fmpz, y: &Fmpz, _tmp: &mut Fmpz) {
        self.submul_mut(x, y);
    }

    fn square_g(&mut self, _tmp_elt: &mut Fmpz, _tmp: &mut Fmpz) {
        self.set_pow_ui(2);
    }
}

/// (rt + ir sqrt(5))/2
#[derive(Clone)]
pub struct Sqrt5Z {
    pub rt: Fmpz,
    pub ir: Fmpz,
}

impl<'a> From<&'a Fmpz> for Sqrt5Z {
    fn from(a: &Fmpz) -> Self {
        let mut rt = Fmpz::new();
        rt.set(&a);
        rt <<= 1;
        Self {
            rt: rt,
            ir: Fmpz::zero(),
        }
    }
}

impl RealQuadElement<Fmpz> for Sqrt5Z {
    fn rt_part(&self) -> Fmpz {
        self.rt.clone()
    }

    fn ir_part(&self) -> Fmpz {
        self.ir.clone()
    }
}

impl<'a> From<&'a (Fmpz, Fmpz)> for Sqrt5Z {
    fn from(tpl: &(Fmpz, Fmpz)) -> Self {
        let mut rt = Fmpz::new();
        let mut ir = Fmpz::new();
        rt.set(&tpl.0);
        ir.set(&tpl.1);
        Self { rt: rt, ir: ir }
    }
}

impl From<(i64, i64)> for Sqrt5Z {
    fn from(a: (i64, i64)) -> Self {
        let mut rt = Fmpz::new();
        let mut ir = Fmpz::new();
        rt.set_si(a.0);
        ir.set_si(a.1);
        Self { rt: rt, ir: ir }
    }
}

impl Sqrt5Z {
    pub fn conj_mut(&mut self) {
        self.ir.negate();
    }

    pub fn norm(&self, res: &mut Fmpz) {
        let &Self {
            rt: ref a,
            ir: ref b,
        } = self;
        res.mul_mut(b, b);
        *res *= -5 as c_long;
        res.addmul_mut(a, a);
        *res >>= 2;
    }

    pub fn from_sisi(rt: c_long, ir: c_long) -> Self {
        Self {
            rt: Fmpz::from_si(rt),
            ir: Fmpz::from_si(ir),
        }
    }

    // Return maximum positive integer `a` such that `self/a` is integral.
    pub fn content(&self) -> Fmpz {
        let a = Fmpz::gcd(&self.rt, &self.ir);
        let mut b = &self.rt - &self.ir;
        b /= &a;
        if !b.is_even() { a >> 1 } else { a }
    }
}

impl BigNumber for Sqrt5Z {
    fn is_zero_g(&self) -> bool {
        self.rt.is_zero() && self.ir.is_zero()
    }

    fn from_ui_g(x: c_ulong) -> Self {
        let mut a = Fmpz::from_ui(x);
        a <<= 1;
        Self {
            rt: a,
            ir: Fmpz::zero(),
        }
    }

    fn from_si_g(x: c_long) -> Self {
        let mut a = Fmpz::from_si(x);
        a <<= 1;
        Self {
            rt: a,
            ir: Fmpz::zero(),
        }
    }

    fn set_ui_g(&mut self, x: c_ulong) {
        self.rt.set_ui(x);
        self.rt <<= 1;
        self.ir.set_ui(0);
    }

    fn set_si_g(&mut self, x: c_long) {
        self.rt.set_si(x);
        self.rt <<= 1;
        self.ir.set_ui(0);
    }

    fn set_g(&mut self, other: &Self) {
        self.ir.set(&other.ir);
        self.rt.set(&other.rt);
    }

    fn add_mut_g(&mut self, x: &Self, y: &Self) {
        self.ir.add_mut(&x.ir, &y.ir);
        self.rt.add_mut(&x.rt, &y.rt);
    }

    fn sub_mut_g(&mut self, x: &Self, y: &Self) {
        self.ir.sub_mut(&x.ir, &y.ir);
        self.rt.sub_mut(&x.rt, &y.rt);
    }

    fn mul_mut_g(&mut self, x: &Self, y: &Self) {
        self.rt.mul_mut(&x.ir, &y.ir);
        self.rt *= 5 as c_ulong;
        self.rt.addmul_mut(&x.rt, &y.rt);
        self.ir.mul_mut(&x.rt, &y.ir);
        self.ir.addmul_mut(&x.ir, &y.rt);
        self.rt >>= 1;
        self.ir >>= 1;
    }

    fn new_g() -> Self {
        Self {
            rt: Fmpz::new(),
            ir: Fmpz::new(),
        }
    }

    fn negate_g(&mut self) {
        self.rt.negate();
        self.ir.negate();
    }

    fn set_divexact_g(&mut self, x: &Self, tmp: &mut Fmpz) {
        // (x.conj() * y/y.norm()).conj()
        self.conj_mut();
        self.mul_assign_g(&x, tmp);
        x.norm(tmp);
        self.ir.set_divexact(tmp);
        self.rt.set_divexact(tmp);
        self.conj_mut();
    }

    fn is_multiple_of_g(&self, x: &Self, tmpelt: &mut Self, tmp: &mut Fmpz) -> bool {
        tmpelt.set_g(self);
        tmpelt.conj_mut();
        tmpelt.mul_assign_g(&x, tmp);
        x.norm(tmp);
        if tmpelt.ir.is_divisible(&tmp) && tmpelt.rt.is_divisible(&tmp) {
            tmpelt.ir /= tmp as &Fmpz;
            tmpelt.rt /= tmp as &Fmpz;
            tmp.add_mut(&tmpelt.ir, &tmpelt.rt);
            // TODO: Use macro Fmpz_even_p.
            tmp.is_even()
        } else {
            false
        }
    }

    fn mul_assign_g(&mut self, other: &Self, tmp: &mut Fmpz) {
        let &mut Self {
            rt: ref mut a,
            ir: ref mut b,
        } = self;
        let &Self {
            rt: ref c,
            ir: ref d,
        } = other;
        tmp.set(b);
        *b *= c;
        b.addmul_mut(a, d);
        *b >>= 1;

        *tmp *= 5 as c_ulong;
        *a *= c;
        a.addmul_mut(tmp, d);
        *a >>= 1;
    }

    fn addmul_mut_g(&mut self, x: &Self, y: &Self, tmp: &mut Fmpz) {
        tmp.mul_mut(&x.ir, &y.ir);
        *tmp *= 5 as c_ulong;
        tmp.addmul_mut(&x.rt, &y.rt);
        *tmp >>= 1;
        self.rt += tmp as &Fmpz;
        tmp.mul_mut(&x.rt, &y.ir);
        tmp.addmul_mut(&x.ir, &y.rt);
        *tmp >>= 1;
        self.ir += tmp as &Fmpz;
    }

    fn submul_mut_g(&mut self, x: &Self, y: &Self, tmp: &mut Fmpz) {
        self.negate_g();
        self.addmul_mut_g(&x, &y, tmp);
        self.negate_g();
    }

    fn square_g(&mut self, tmp_elt: &mut Self, tmp: &mut Fmpz) {
        tmp_elt.set_g(self);
        self.mul_assign_g(tmp_elt, tmp);
    }
}

impl fmt::Display for Sqrt5Z {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.ir.is_zero() {
            write!(f, "{}", &self.rt >> 1)
        } else {
            let mut tmp_ply = FmpzPoly::new();
            if self.ir.is_even() {
                tmp_ply.set_coeff(&(&self.rt >> 1), 0);
                tmp_ply.set_coeff(&(&self.ir >> 1), 1);
                write!(f, "{}", tmp_ply)
            } else {
                tmp_ply.set_coeff(&(self.rt), 0);
                tmp_ply.set_coeff(&(&self.ir), 1);
                write!(f, "({})/2", tmp_ply)
            }
        }
    }
}

impl fmt::Debug for Sqrt5Z {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl<'a> AddAssign<&'a Sqrt5Z> for Sqrt5Z {
    fn add_assign(&mut self, other: &Sqrt5Z) {
        self.rt += &other.rt;
        self.ir += &other.ir;
    }
}

impl<'a> SubAssign<&'a Sqrt5Z> for Sqrt5Z {
    fn sub_assign(&mut self, other: &Sqrt5Z) {
        self.rt -= &other.rt;
        self.ir -= &other.ir;
    }
}

impl MulAssign<c_long> for Sqrt5Z {
    fn mul_assign(&mut self, other: c_long) {
        self.ir *= other;
        self.rt *= other;
    }
}

impl MulAssign<c_ulong> for Sqrt5Z {
    fn mul_assign(&mut self, other: c_ulong) {
        self.ir *= other;
        self.rt *= other;
    }
}

impl ShlAssign<u64> for Sqrt5Z {
    fn shl_assign(&mut self, other: u64) {
        self.rt <<= other;
        self.ir <<= other;
    }
}

impl ShrAssign<u64> for Sqrt5Z {
    fn shr_assign(&mut self, other: u64) {
        self.rt >>= other;
        self.ir >>= other;
    }
}

impl PartialEq for Sqrt5Z {
    fn eq(&self, other: &Self) -> bool {
        self.rt == other.rt && self.ir == other.ir
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_divexact_g() {
        let a = Sqrt5Z {
            rt: Fmpz::from_ui(1),
            ir: Fmpz::from_ui(3),
        };
        let mut b = Sqrt5Z {
            rt: Fmpz::from_ui(22),
            ir: Fmpz::from_si(0),
        };
        let mut tmp = Fmpz::new();
        b.set_divexact_g(&a, &mut tmp);
        assert_eq!(
            b,
            Sqrt5Z {
                rt: Fmpz::from_si(-1),
                ir: Fmpz::from_si(3),
            }
        );
    }

    #[test]
    fn test_set_fun() {
        let a = Sqrt5Z::from_si_g(3);
        let mut b = Sqrt5Z::new_g();
        b.set_si_g(3);
        assert_eq!(a, b);
    }

    #[test]
    fn test_addmul_mut() {
        let mut tmp = Fmpz::new();
        let a = Sqrt5Z::from_sisi(3, 5);
        let b = Sqrt5Z::from_sisi(7, 1);
        let mut res = Sqrt5Z::from_sisi(4, 6);
        res.addmul_mut_g(&a, &b, &mut tmp);
        assert_eq!(res, Sqrt5Z::from_sisi(27, 25));
    }

    #[test]
    fn test_mul_mut() {
        let mut res = Sqrt5Z::from_sisi(2, 4);
        let mut tmp = Fmpz::new();
        let mut a = Sqrt5Z {
            rt: Fmpz::from_si(5),
            ir: Fmpz::from_si(1),
        };
        let b = Sqrt5Z {
            rt: Fmpz::from_si(3),
            ir: Fmpz::from_si(7),
        };
        res.mul_mut_g(&a, &b);
        a.mul_assign_g(&b, &mut tmp);
        assert_eq!(res, a);
    }

    #[test]
    fn test_is_multiple_of() {
        let ref mut tmpelt = Sqrt5Z::new_g();
        let ref mut tmp = Fmpz::new();
        let a = Sqrt5Z::from_ui_g(10);
        let b = Sqrt5Z::from_ui_g(2);
        assert!(a.is_multiple_of_g(&b, tmpelt, tmp));
        let b = Sqrt5Z::from_ui_g(3);
        assert!(!a.is_multiple_of_g(&b, tmpelt, tmp));
        let a = Sqrt5Z {
            rt: Fmpz::from_si(2),
            ir: Fmpz::from_si(4),
        };
        let b = Sqrt5Z::from_ui_g(2);
        assert!(!a.is_multiple_of_g(&b, tmpelt, tmp));
        let a = Sqrt5Z::from_ui_g(33);
        let b = Sqrt5Z {
            rt: Fmpz::from_si(-1),
            ir: Fmpz::from_si(3),
        };
        assert!(a.is_multiple_of_g(&b, tmpelt, tmp));
    }

    #[test]
    fn test_submul_mut() {
        let mut a = Sqrt5Z::from_sisi(3, 3);
        let b = Sqrt5Z::from_sisi(1, 1);
        let c = Sqrt5Z::from_sisi(2, 4);
        let mut tmp = Fmpz::new();
        a.submul_mut_g(&b, &c, &mut tmp);
        assert_eq!(a, Sqrt5Z::from_si_g(-4));
    }

    #[test]
    fn test_pow() {
        let mut a = Sqrt5Z::new_g();
        let b = Sqrt5Z::from_sisi(3, 5);
        a.pow_mut(&b, 10);
        assert_eq!(a.rt_part().get_string(10), "322355827");
        assert_eq!(a.ir_part().get_string(10), "142989825");
    }

    #[test]
    fn test_fmt() {
        let a: Sqrt5Z = From::from((1, -1));
        println!("{}", a);
        let a: Sqrt5Z = From::from((1, 1));
        println!("{}", a);
        let a: Sqrt5Z = From::from((2, -2));
        println!("{}", a);
    }

    #[test]
    fn test_content() {
        let a: Sqrt5Z = From::from((3, 3));
        assert_eq!(a.content(), 3_i64);
        let a: Sqrt5Z = From::from((10, 4));
        assert_eq!(a.content(), 1_i64);
    }
}
