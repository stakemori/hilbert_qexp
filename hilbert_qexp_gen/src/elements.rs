use std;
use std::fmt;
use misc::{int_sqrt, pow_mut, PowGen};
use libc::{c_long, c_ulong};
use std::ops::{Add, AddAssign, DivAssign, Mul, MulAssign, Neg, ShlAssign, ShrAssign, Sub,
               SubAssign};
use std::cmp::min;
use bignum::{BigNumber, RealQuadElement};

use fcvec;
use flint::fmpz::Fmpz;
use flint::fmpq::Fmpq;
use flint::traits::*;
use flint::fmpz_mat::FmpzMat;
use flint::fmpq_mat::FmpqMat;

type Weight = Option<(usize, usize)>;
/// struct for hilbert modualr form over Q(sqrt(m))
/// this corresponds finite sum of the q-expansion of the form
/// Σ a(u, v) exp(2piTr 1/sqrt(m) (u + v * sqrt(m))/2 z) if m equiv 1 mod 4
/// Σ a(u, v) exp(2piTr 1/2sqrt(m) (u + v * sqrt(m)) z) otherwise
/// where v <= prec.
/// a(u, v) = fc[v][a], where a = u + u_bds[v]
#[derive(Debug, Clone)]
pub struct HmfGen<T> {
    pub prec: usize,
    pub fcvec: FcVec<T>,
    pub weight: Weight,
    // vth element of u_bds.vec is (sqrt(m) * v).floor()
    pub u_bds: UBounds,
    // Square free positive integer.
    pub m: u64,
}

#[macro_export]
macro_rules! is_1mod4 {
    ($expr: expr) => {($expr) & 0b11 == 1}
}

#[macro_export]
macro_rules! is_even {
    ($expr: expr) => {($expr) & 1 == 0}
}

#[macro_export]
macro_rules! u_iter {
    ($m: expr, $v: expr, $bd: ident) => {
        {
            (-$bd..($bd+1)).filter(|&x| $m & 0b11 != 1 || is_even!(x+$v))
        }
    }
}

#[macro_export]
macro_rules! v_u_bd_iter {
    (($m: expr, $u_bds: expr, $v: ident, $u: ident, $bd: ident) $body:expr) =>
    {
        for ($v, &$bd) in $u_bds.vec.iter().enumerate() {
            let $bd = $bd as i64;
            let v_i64 = $v as i64;
            let m = $m;
            for $u in u_iter!(m, v_i64, $bd) {
                $body
            }
        };
    }
}

macro_rules! v_u_bd_iter_non_const {
    (($m: expr, $u_bds: expr, $v: ident, $u: ident, $bd: ident) $body:expr) =>
    {
        for ($v, &$bd) in $u_bds.vec.iter().enumerate().skip(1) {
            let $bd = $bd as i64;
            let v_i64 = $v as i64;
            let m = $m;
            for $u in u_iter!(m, v_i64, $bd) {
                $body
            }
        };
    }
}

impl<T> fmt::Display for HmfGen<T>
where
    T: BigNumber + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut vec = Vec::new();
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            let a = self.fcvec.fc_ref(v, u, bd);
            if !a.is_zero_g() {
                vec.push(format!("({}, {}): {}", u, v, a));
            }
        }
        );
        write!(f, "{}", vec.join("\n"))
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct FcVec<T> {
    pub vec: Vec<Vec<T>>,
}

impl<'a, T, S> From<&'a FcVec<S>> for FcVec<T>
where
    for<'b> T: From<&'b S>,
{
    fn from(a: &FcVec<S>) -> FcVec<T> {
        Self {
            vec: a.vec
                .iter()
                .map(|v| v.iter().map(From::from).collect())
                .collect(),
        }
    }
}

impl<T, S> RealQuadElement<FcVec<S>> for FcVec<T>
where
    T: RealQuadElement<S>,
{
    fn ir_part(&self) -> FcVec<S> {
        let vec = self.vec
            .iter()
            .map(|v| v.iter().map(|x| x.ir_part()).collect())
            .collect();
        FcVec::<S> { vec: vec }
    }

    fn rt_part(&self) -> FcVec<S> {
        let vec = self.vec
            .iter()
            .map(|v| v.iter().map(|x| x.rt_part()).collect())
            .collect();
        FcVec::<S> { vec: vec }
    }
}

impl<T> FcVec<T>
where
    T: BigNumber,
{
    pub fn fc_ref(&self, v: usize, u: i64, bd: i64) -> &T {
        debug_assert!(u + bd >= 0);
        self.vec[v].get((u + bd) as usize).unwrap()
    }

    pub fn fc_ref_mut(&mut self, v: usize, u: i64, bd: i64) -> &mut T {
        debug_assert!(u + bd >= 0);
        self.vec[v].get_mut((u + bd) as usize).unwrap()
    }

    fn new(u_bds: &UBounds) -> FcVec<T> {
        let vec = u_bds
            .vec
            .iter()
            .map(|&bd| (0..(2 * bd + 1)).map(|_| T::from_ui_g(0)).collect())
            .collect();
        FcVec { vec: vec }
    }
}

#[derive(Debug, Clone)]
pub struct UBounds {
    pub vec: Vec<usize>,
}

impl UBounds {
    // m is a square free positive integer.
    pub fn new(m: u64, prec: usize) -> UBounds {
        let mut u_bds = Vec::new();
        for v in 0..(prec + 1) {
            u_bds.push(int_sqrt(m * v as u64 * v as u64) as usize);
        }
        UBounds { vec: u_bds }
    }

    pub fn take(&self, n: usize) -> UBounds {
        let v = self.vec.iter().cloned().take(n).collect();
        UBounds { vec: v }
    }
}

impl<T> PowGen for HmfGen<T>
where
    T: BigNumber + Clone,
    for<'a> T: AddAssign<&'a T>,
{
    fn set_one(&mut self) {
        v_u_bd_iter_non_const!((self.m, self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).set_ui_g(0);
        });
        self.fcvec.fc_ref_mut(0, 0, 0).set_ui_g(1);
        self.weight = Some((0, 0));
    }

    fn square(&mut self) {
        let f = self.clone();
        self.weight = weight_pow(self.weight, 2);
        let mut tmp = T::from_ui_g(0);
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            _mul_mut_tmp(&mut tmp, u, v, &f.fcvec, &f.fcvec, &self.u_bds, self.m);
            self.fcvec.fc_ref_mut(v, u, bd).set_g(&tmp);
            })
    }
}

fn weight_mul(a: Weight, b: Weight) -> Weight {
    a.and_then(|x| b.and_then(|y| Some((x.0 + y.0, x.1 + y.1))))
}

fn weight_add(a: Weight, b: Weight) -> Weight {
    a.and_then(|x| b.and_then(|y| if x == y { Some(x) } else { None }))
}

fn weight_pow(a: Weight, n: usize) -> Weight {
    a.and_then(|(k1, k2)| Some((k1 * n, k2 * n)))
}

pub fn weight_div(a: Weight, b: Weight) -> Weight {
    a.and_then(|x| b.and_then(|y| Some((x.0 - y.0, x.1 - y.1))))
}

impl<'a, T, S> From<&'a HmfGen<S>> for HmfGen<T>
where
    for<'b> T: From<&'b S>,
{
    fn from(f: &HmfGen<S>) -> Self {
        let fcvec = From::from(&f.fcvec);
        let u_bds = f.u_bds.clone();
        Self {
            fcvec: fcvec,
            weight: f.weight,
            prec: f.prec,
            u_bds: u_bds,
            m: f.m,
        }
    }
}

impl<T> From<HmfGen<Fmpz>> for HmfGen<T>
where
    for<'b> T: From<&'b Fmpz>,
    T: RealQuadElement<Fmpz>,
{
    fn from(f: HmfGen<Fmpz>) -> Self {
        From::from(&f)
    }
}

impl<T, S> RealQuadElement<HmfGen<S>> for HmfGen<T>
where
    T: RealQuadElement<S>,
{
    fn rt_part(&self) -> HmfGen<S> {
        HmfGen::<S> {
            weight: self.weight,
            prec: self.prec,
            fcvec: self.fcvec.rt_part(),
            u_bds: self.u_bds.clone(),
            m: self.m,
        }
    }

    fn ir_part(&self) -> HmfGen<S> {
        HmfGen::<S> {
            weight: self.weight,
            prec: self.prec,
            fcvec: self.fcvec.ir_part(),
            u_bds: self.u_bds.clone(),
            m: self.m,
        }
    }
}

impl<T> HmfGen<T>
where
    T: BigNumber,
{
    /// Return 0 q-expantion
    pub fn new(m: u64, prec: usize) -> HmfGen<T> {
        let u_bds = UBounds::new(m, prec);
        let fcvec = FcVec::new(&u_bds);
        HmfGen {
            weight: None,
            prec: prec,
            fcvec: fcvec,
            u_bds: u_bds,
            m: m,
        }
    }

    /// Decrease prec to prec.
    pub fn decrease_prec(&mut self, prec: usize) {
        assert!(self.prec >= prec);
        self.u_bds = self.u_bds.take(prec + 1);
        self.prec = prec;
    }

    pub fn one(m: u64, prec: usize) -> HmfGen<T> {
        let mut f = Self::new(m, prec);
        f.fcvec.fc_ref_mut(0, 0, 0).set_ui_g(1);
        f
    }

    /// set self = f1 + f2
    pub fn add_mut(&mut self, f1: &HmfGen<T>, f2: &HmfGen<T>) {
        let prec = min(f1.prec, f2.prec);
        self.decrease_prec(prec);
        self.weight = weight_add(f1.weight, f2.weight);
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            T::add_mut_g(
                self.fcvec.fc_ref_mut(v, u, bd),
                f1.fcvec.fc_ref(v, u, bd),
                f2.fcvec.fc_ref(v, u, bd),
            );
        })
    }

    pub fn is_divisible_by_const(&self, a: &T) -> bool {
        let mut tmpelt = T::new_g();
        let mut tmp = T::R::default();
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            if !self.fcvec.fc_ref(v, u, bd).is_multiple_of_g(
                a,
                &mut tmpelt,
                &mut tmp) {
                return false;
            }
        });
        true
    }

    /// set self = f1 - f2
    pub fn sub_mut(&mut self, f1: &HmfGen<T>, f2: &HmfGen<T>) {
        let prec = min(f1.prec, f2.prec);
        self.decrease_prec(prec);
        self.weight = weight_add(f1.weight, f2.weight);
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            T::sub_mut_g(
                self.fcvec.fc_ref_mut(v, u, bd),
                f1.fcvec.fc_ref(v, u, bd),
                f2.fcvec.fc_ref(v, u, bd),
            );
        })
    }

    /// set self = f1 * f2
    pub fn mul_mut(&mut self, f1: &HmfGen<T>, f2: &HmfGen<T>)
    where
        for<'a> T: AddAssign<&'a T>,
    {
        let mut tmp = T::from_ui_g(0);
        let prec = min(f1.prec, f2.prec);
        self.decrease_prec(prec);
        self.weight = weight_mul(f1.weight, f2.weight);
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            _mul_mut_tmp(&mut tmp, u, v, &f1.fcvec, &f2.fcvec, &self.u_bds, self.m);
            self.fcvec.fc_ref_mut(v, u, bd).set_g(&tmp);
        })
    }

    /// self = f * a
    pub fn mul_mut_by_const(&mut self, f: &HmfGen<T>, a: &T) {
        self.weight = f.weight;
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            T::mul_mut_g(self.fcvec.fc_ref_mut(v, u, bd), f.fcvec.fc_ref(v, u, bd), a)
            })
    }

    pub fn pow_mut(&mut self, f: &HmfGen<T>, a: usize)
    where
        T: Clone,
        for<'a> T: AddAssign<&'a T>,
    {
        if a == 0 {
            self.set_one();
        } else if a == 1 {
            self.set(&f);
        } else {
            pow_mut(self, f, a)
        };
        self.weight = weight_pow(f.weight, a);
    }

    pub fn pow(&self, a: usize) -> Self
    where
        T: Clone,
        for<'a> T: AddAssign<&'a T>,
    {
        let mut tmp = Self::new(self.m, self.prec);
        tmp.pow_mut(self, a);
        tmp
    }

    pub fn is_zero(&self) -> bool {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            if !self.fcvec.fc_ref(v, u, bd).is_zero_g() {
                return false;
            }
        });
        true
    }

    pub fn fourier_coefficient(&self, v: usize, u: i64) -> T {
        let bd = self.u_bds.vec[v] as i64;
        let mut a = T::new_g();
        a.set_g(self.fcvec.fc_ref(v, u, bd));
        a
    }

    pub fn fourier_coefficients(&self, vec: &Vec<(usize, i64)>) -> Vec<T> {
        vec.iter()
            .map(|&(v, u)| self.fourier_coefficient(v, u))
            .collect()
    }

    pub fn fc_vector(&self, len: usize) -> Vec<T> {
        let mut vu_vec = Vec::new();
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            vu_vec.push((v, u));
        });
        let res: Vec<_> = vu_vec
            .iter()
            .take(len)
            .map(|&(v, u)| self.fourier_coefficient(v, u))
            .collect();
        assert_eq!(res.len(), len);
        res
    }

    pub fn fc_vector_all(&self) -> Vec<T>
    where
        T: Clone,
    {
        self.fcvec.vec.clone().into_iter().flat_map(|v| v).collect()
    }

    pub fn fc_vector_u_nonneg(&self) -> Vec<T>
    where
        T: Clone,
    {
        let mut vu_vec = Vec::new();
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            if u >= 0 {
                vu_vec.push((v, u));
            }
        });
        vu_vec
            .iter()
            .map(|&(v, u)| self.fourier_coefficient(v, u))
            .collect()
    }

    pub fn fc_vector_real_quad(&self, len: usize) -> Vec<(Fmpz, Fmpz)>
    where
        T: RealQuadElement<Fmpz>,
    {
        self.fc_vector(len)
            .iter()
            .map(|x| (x.rt_part(), x.ir_part()))
            .collect()
    }

    pub fn set(&mut self, other: &Self) {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).set_g(other.fcvec.fc_ref(v, u, bd));
        })
    }

    pub fn set_zero(&mut self) {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).set_ui_g(0);
        })
    }

    pub fn negate(&mut self) {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).negate_g();
        })
    }

    pub fn diagonal_restriction(&self) -> Vec<T>
    where
        for<'a> T: AddAssign<&'a T>,
    {
        let prec = self.prec;
        let mut vec = Vec::new();
        for _ in 0..(prec + 1) {
            let mut a = T::new_g();
            a.set_ui_g(0);
            vec.push(a);
        }
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            vec[v] += self.fcvec.fc_ref(v, u, bd);
        });
        vec
    }
}

impl HmfGen<Fmpz> {
    pub fn gcd(&self) -> Fmpz {
        let mut res = Fmpz::new();
        let mut tmp = Fmpz::new();
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            tmp.set(&res);
            res.gcd_mut(&tmp, self.fcvec.fc_ref(v, u, bd));
        });
        res
    }
}

impl<T> HmfGen<T>
where
    T: BigNumber + RealQuadElement<Fmpz>,
{
    pub fn gcd(&self) -> Fmpz {
        let mut res = Fmpz::new();
        let mut tmp = Fmpz::new();
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            tmp.set(&res);
            res.gcd_mut(&tmp, &self.fcvec.fc_ref(v, u, bd).ir_part());
            tmp.set(&res);
            res.gcd_mut(&tmp, &self.fcvec.fc_ref(v, u, bd).rt_part());
        });
        let mut a = false;
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            tmp.add_mut(&self.fcvec.fc_ref(v, u, bd).ir_part(),
                        &self.fcvec.fc_ref(v, u, bd).rt_part());
            tmp /= &res;
            if self.m & 0b11 == 1 && !tmp.is_even() {
                a = true;
                break;
            }
        });
        if a {
            &res >> 1
        } else {
            res
        }
    }
}

impl<'a, T> DivAssign<&'a T> for HmfGen<T>
where
    T: BigNumber,
{
    fn div_assign(&mut self, num: &T) {
        let mut tmp = T::R::default();
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).set_divexact_g(num, &mut tmp);
        })
    }
}

impl<'a, T> AddAssign<&'a HmfGen<T>> for HmfGen<T>
where
    for<'b> T: AddAssign<&'b T>,
    T: BigNumber + Clone,
{
    fn add_assign(&mut self, other: &HmfGen<T>) {
        if self.is_zero() {
            self.set(other);
            self.weight = other.weight;
        } else {
            self.weight = weight_add(self.weight, other.weight);
            let prec = min(self.prec, other.prec);
            self.decrease_prec(prec);
            v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
                T::add_assign(
                    self.fcvec.fc_ref_mut(v, u, bd),
                    other.fcvec.fc_ref(v, u, bd),
                );
            })
        }
    }
}

impl<'a, T> SubAssign<&'a HmfGen<T>> for HmfGen<T>
where
    for<'b> T: SubAssign<&'b T>,
    T: BigNumber + Clone,
{
    fn sub_assign(&mut self, other: &HmfGen<T>) {
        self.weight = weight_add(self.weight, other.weight);
        let prec = min(self.prec, other.prec);
        self.decrease_prec(prec);
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            *self.fcvec.fc_ref_mut(v, u, bd) -= other.fcvec.fc_ref(v, u, bd);
        })
    }
}

impl<'a, 'b, T> Add<&'a HmfGen<T>> for &'b HmfGen<T>
where
    T: BigNumber + Clone,
{
    type Output = HmfGen<T>;
    fn add(self, other: &HmfGen<T>) -> HmfGen<T> {
        let prec = min(self.prec, other.prec);
        let mut res = HmfGen::new(self.m, prec);
        res.add_mut(self, other);
        res
    }
}

macro_rules! impl_op_take {
    ($tr: ident, $mth: ident) => {
        impl<T> $tr<HmfGen<T>> for HmfGen<T>
            where T: BigNumber + Clone,
        for<'c> T: AddAssign<&'c T>,
        {
            type Output = HmfGen<T>;
            fn $mth(self, other: Self) -> Self {
                let a = &self;
                a.$mth(&other)
            }
        }

        impl<'a, T> $tr<&'a HmfGen<T>> for HmfGen<T>
            where T: BigNumber + Clone,
        for<'c> T: AddAssign<&'c T>,
        {
            type Output = HmfGen<T>;
            fn $mth(self, other: &Self) -> Self {
                let a = &self;
                a.$mth(other)
            }
        }

        impl<'a, T> $tr<HmfGen<T>> for &'a HmfGen<T>
            where T: BigNumber + Clone,
        for<'c> T: AddAssign<&'c T>,
        {
            type Output = HmfGen<T>;
            fn $mth(self, other: HmfGen<T>) -> HmfGen<T> {
                self.$mth(&other)
            }
        }

    }
}
impl_op_take!(Add, add);
impl_op_take!(Mul, mul);
impl_op_take!(Sub, sub);

impl<'a, 'b, T> Mul<&'a HmfGen<T>> for &'b HmfGen<T>
where
    T: BigNumber + Clone,
    for<'c> T: AddAssign<&'c T>,
{
    type Output = HmfGen<T>;
    fn mul(self, other: &HmfGen<T>) -> HmfGen<T> {
        let prec = min(self.prec, other.prec);
        let mut res = HmfGen::new(self.m, prec);
        res.mul_mut(self, other);
        res
    }
}

impl<'a, 'b, T> Mul<&'a T> for &'b HmfGen<T>
where
    T: BigNumber + Clone,
{
    type Output = HmfGen<T>;
    fn mul(self, other: &T) -> HmfGen<T> {
        let mut res = HmfGen::new(self.m, self.prec);
        res.mul_mut_by_const(self, other);
        res
    }
}

impl<'a, T> Mul<&'a T> for HmfGen<T>
where
    T: BigNumber + Clone,
{
    type Output = HmfGen<T>;
    fn mul(self, other: &T) -> HmfGen<T> {
        let mut res = HmfGen::new(self.m, self.prec);
        res.mul_mut_by_const(&self, other);
        res
    }
}

impl<'a, T> Mul<c_long> for &'a HmfGen<T>
where
    T: BigNumber + Clone,
{
    type Output = HmfGen<T>;
    fn mul(self, other: c_long) -> HmfGen<T> {
        let mut res = HmfGen::new(self.m, self.prec);
        res.mul_mut_by_const(self, &T::from_si_g(other));
        res
    }
}

impl<'b, T> Mul<&'b Fmpz> for HmfGen<T>
where
    for<'c> T: From<&'c Fmpz>,
    T: BigNumber + Clone + RealQuadElement<Fmpz>,
{
    type Output = HmfGen<T>;
    fn mul(self, other: &Fmpz) -> HmfGen<T> {
        let a: T = From::from(other);
        &self * &a
    }
}

impl<T> Mul<c_long> for HmfGen<T>
where
    T: BigNumber + Clone,
{
    type Output = HmfGen<T>;
    fn mul(self, other: c_long) -> HmfGen<T> {
        let mut res = HmfGen::new(self.m, self.prec);
        res.mul_mut_by_const(&self, &T::from_si_g(other));
        res
    }
}

impl<'a, 'b, T> Sub<&'a HmfGen<T>> for &'b HmfGen<T>
where
    T: BigNumber + Clone,
{
    type Output = HmfGen<T>;
    fn sub(self, other: &HmfGen<T>) -> HmfGen<T> {
        let prec = min(self.prec, other.prec);
        let mut res = HmfGen::new(self.m, prec);
        res.sub_mut(self, other);
        res
    }
}

impl<'a, T> Neg for &'a HmfGen<T>
where
    T: BigNumber + Clone,
{
    type Output = HmfGen<T>;
    fn neg(self) -> HmfGen<T> {
        let mut res = HmfGen::new(self.m, self.prec);
        let mone = T::from_si_g(-1);
        res.mul_mut_by_const(self, &mone);
        res
    }
}

impl<T> PartialEq for HmfGen<T>
where
    T: BigNumber + PartialEq,
{
    fn eq(&self, other: &HmfGen<T>) -> bool {
        let prec = min(self.prec, other.prec);
        v_u_bd_iter!((self.m, self.u_bds.take(prec + 1), v, u, bd) {
            if self.fcvec.fc_ref(v, u, bd) != other.fcvec.fc_ref(v, u, bd) {
                return false;
            }
        }
        );
        true
    }
}

/// set (v, u) th F.C. of fc_vec1 * fc_vec2 to a.
/// This function take care the case when fc_vec2 is sparse.
fn _mul_mut_tmp<T>(
    a: &mut T,
    u: i64,
    v: usize,
    fc_vec1: &FcVec<T>,
    fc_vec2: &FcVec<T>,
    u_bds: &UBounds,
    m: u64,
) where
    for<'a> T: AddAssign<&'a T>,
    T: BigNumber,
{
    a.set_ui_g(0);
    let mut tmp = T::new_g();
    for v2 in 0..(v + 1) {
        let bd2 = u_bds.vec[v2] as i64;
        let v2_i64 = v2 as i64;
        for u2 in u_iter!(m, v2_i64, bd2) {
            if !fc_vec2.fc_ref(v2, u2, bd2).is_zero_g() {
                let u1 = u - u2;
                let v1 = v - v2;
                let u1abs = u1.abs() as usize;
                if u1abs * u1abs <= m as usize * v1 * v1 {
                    let bd1 = u_bds.vec[v1] as i64;
                    tmp.mul_mut_g(fc_vec2.fc_ref(v2, u2, bd2), fc_vec1.fc_ref(v1, u1, bd1));
                    T::add_assign(a, &tmp);
                }
            }
        }
    }
}

impl<'a, T> MulAssign<&'a HmfGen<T>> for HmfGen<T>
where
    T: BigNumber + Clone,
    for<'b> T: AddAssign<&'b T>,
{
    fn mul_assign(&mut self, other: &HmfGen<T>) {
        let prec = min(self.prec, other.prec);
        self.weight = weight_mul(self.weight, other.weight);
        self.decrease_prec(prec);
        // We need cloned self.
        let f = self.clone();
        let mut tmp = T::from_ui_g(0);
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            _mul_mut_tmp(&mut tmp, u, v, &f.fcvec, &other.fcvec, &self.u_bds, self.m);
            self.fcvec.fc_ref_mut(v, u, bd).set_g(&tmp);
            })
    }
}

impl<'a, T> MulAssign<&'a Fmpq> for HmfGen<T>
where
    T: BigNumber,
    T: RealQuadElement<Fmpq>,
    for<'b> T: MulAssign<&'b Fmpq>,
{
    fn mul_assign(&mut self, other: &Fmpq) {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            *self.fcvec.fc_ref_mut(v, u, bd) *= other;
        }
        );
    }
}

impl<'a, T> MulAssign<&'a T> for HmfGen<T>
where
    T: BigNumber,
{
    fn mul_assign(&mut self, other: &T) {
        let mut tmp = T::R::default();
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).mul_assign_g(&other, &mut tmp);
        }
        )
    }
}

impl<T> MulAssign<c_ulong> for HmfGen<T>
where
    T: BigNumber + MulAssign<c_ulong>,
{
    fn mul_assign(&mut self, other: c_ulong) {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            *self.fcvec.fc_ref_mut(v, u, bd) *= other;
        }
        );
    }
}

impl<T> MulAssign<c_long> for HmfGen<T>
where
    T: BigNumber + MulAssign<c_long>,
{
    fn mul_assign(&mut self, other: c_long) {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            *self.fcvec.fc_ref_mut(v, u, bd) *= other;
        }
        );
    }
}

impl<T, S> ShlAssign<S> for HmfGen<T>
where
    T: BigNumber + ShlAssign<S>,
    S: Copy,
{
    fn shl_assign(&mut self, other: S) {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            T::shl_assign(self.fcvec.fc_ref_mut(v, u, bd), other);
        }
        );
    }
}

impl<T, S> ShrAssign<S> for HmfGen<T>
where
    T: BigNumber + ShrAssign<S>,
    S: Copy,
{
    fn shr_assign(&mut self, other: S) {
        v_u_bd_iter!((self.m, self.u_bds, v, u, bd) {
            T::shr_assign(self.fcvec.fc_ref_mut(v, u, bd), other);
        }
        );
    }
}

pub fn initial_term<T>(f: &HmfGen<T>) -> Option<(usize, i64, &T)>
where
    T: BigNumber,
{
    for (v, &bd) in f.u_bds.vec.iter().enumerate() {
        let bd = bd as i64;
        let v_i64 = v as i64;
        for u in (-bd..(bd + 1))
            .rev()
            .filter(|&x| !is_1mod4!(f.m) || is_even!(x + v_i64))
        {
            if !f.fcvec.fc_ref(v, u, bd).is_zero_g() {
                return Some((v, u, f.fcvec.fc_ref(v, u, bd)));
            }
        }
    }
    None
}

#[allow(dead_code)]
fn print_vth_cf<T>(f: &HmfGen<T>, v: usize, s: &'static str)
where
    T: BigNumber + std::fmt::Display,
{
    let bd = f.u_bds.vec[v] as i64;
    let mut res = Vec::new();
    for u in -bd..(bd + 1) {
        let a = f.fcvec.fc_ref(v, u, bd);
        if !a.is_zero_g() {
            res.push(format!("({}) * q1**({})", a, u));
        }
    }
    println!("{}: {}", s, res.join(" + "));
}

/// Assuming f is divisible by g in coefficeint T, set res = f/g.
pub fn div_mut<T>(res: &mut HmfGen<T>, f: &HmfGen<T>, g: &HmfGen<T>)
where
    T: BigNumber + Clone,
    for<'a> T: SubAssign<&'a T>,
{
    let (v_init, _, _) = initial_term(&g).unwrap();
    let mut tmp = HmfGen::<T>::new(res.m, g.prec);
    let mut f_cloned = f.clone();
    let prec = f.prec;
    assert!(prec >= v_init);

    res.decrease_prec(prec - v_init);
    res.weight = weight_div(f.weight, g.weight);
    res.set_zero();
    let ref u_bds = f.u_bds;
    for v in v_init..(prec + 1) {
        for i in (1 + v_init)..(v + 1) {
            fcvec::mul_mut(
                &mut tmp.fcvec.vec[v],
                &g.fcvec.vec[i],
                &res.fcvec.vec[v - i],
                i,
                v - i,
                u_bds.vec[i],
                u_bds.vec[v - i],
                u_bds.vec[v],
                u_bds,
                0,
                0,
                0,
                res.m,
            );
            fcvec::sub_assign(
                &mut f_cloned.fcvec.vec[v],
                &tmp.fcvec.vec[v],
                v,
                u_bds,
                res.m,
            );
        }
        fcvec::div_mut(
            &f_cloned.fcvec.vec[v],
            &g.fcvec.vec[v_init],
            &mut res.fcvec.vec[v - v_init],
            v_init,
            v - v_init,
            u_bds.vec[v_init],
            u_bds.vec[v - v_init],
            u_bds.vec[v],
            &u_bds,
        );
    }
}

/// Set `res` to a square root of `f`. `v_init` is the initial exponent of `q = exp(2pi v z)`.
/// And we assume `res[n]` is already set.
pub fn square_root_mut<T>(res: &mut HmfGen<T>, v_init: usize, f: &HmfGen<T>)
where
    T: BigNumber + Clone,
    for<'a> T: SubAssign<&'a T>,
    for<'a> T: std::ops::ShrAssign<u64>,
{
    let mut tmp = HmfGen::<T>::new(res.m, f.prec);
    let mut f_cloned = f.clone();
    let prec = f.prec;
    assert!(prec >= v_init);

    res.decrease_prec(prec - v_init);
    if let Some((k1, k2)) = f.weight {
        res.weight = Some((k1 >> 1, k2 >> 1));
    } else {
        res.weight = None;
    }
    let init_vec = res.fcvec.vec[v_init].clone();
    let ref u_bds = f.u_bds;

    for v in (2 * v_init + 1)..(prec + 1) {
        for i in (1 + v_init)..(v - v_init) {
            fcvec::mul_mut(
                &mut tmp.fcvec.vec[v],
                &res.fcvec.vec[i],
                &res.fcvec.vec[v - i],
                i,
                v - i,
                u_bds.vec[i],
                u_bds.vec[v - i],
                u_bds.vec[v],
                u_bds,
                0,
                0,
                0,
                res.m,
            );
            fcvec::sub_assign(
                &mut f_cloned.fcvec.vec[v],
                &tmp.fcvec.vec[v],
                v,
                u_bds,
                res.m,
            );
        }
        fcvec::div_mut(
            &f_cloned.fcvec.vec[v],
            &init_vec,
            &mut res.fcvec.vec[v - v_init],
            v_init,
            v - v_init,
            u_bds.vec[v_init],
            u_bds.vec[v - v_init],
            u_bds.vec[v],
            &u_bds,
        );
        fcvec::shr_assign(&mut res.fcvec.vec[v - v_init], v - v_init, &u_bds, 1, res.m);
    }
}

// Todo: make the denominator small.
pub fn div_mut_with_denom<T>(res: &mut HmfGen<T>, f: &HmfGen<T>, g: &HmfGen<T>, check: bool) -> T
where
    T: BigNumber + Clone + PartialEq + fmt::Debug,
    for<'a> T: SubAssign<&'a T>,
    for<'a> T: AddAssign<&'a T>,
{
    let mut dnm = T::new_g();
    let mut count = 0;
    let prec = f.prec;
    let (v_init, _, a) = initial_term(&g).unwrap();
    let u_bds = &f.u_bds;
    for v in v_init..(prec + 1) {
        let v_i = v as i64;
        let bd_i = u_bds.vec[v] as i64;
        for _ in u_iter!(res.m, v_i, bd_i) {
            count += 1;
        }
    }
    dnm.pow_mut(a, count);
    let mut f = f.clone();
    f *= &dnm;
    div_mut(res, &f, g);
    if check {
        let f1 = res as &HmfGen<T> * g;
        let mut f = f.clone();
        f.decrease_prec(f1.prec);
        assert_eq!(f, f1);
    }
    dnm
}

pub fn relations_over_z(forms: &[HmfGen<Fmpz>]) -> Vec<Vec<Fmpz>> {
    let vv: Vec<_> = forms.iter().map(|f| f.fc_vector_u_nonneg()).collect();
    let n = forms.len();
    let m = vv[0].len();
    let mut mat = FmpzMat::new(m as i64, n as i64);
    for (i, v) in vv.iter().enumerate() {
        for (j, x) in v.iter().enumerate() {
            mat.set_entry(j as isize, i as isize, x);
        }
    }
    mat.nullspace_basis()
}

pub fn relations_over_q(forms: &[HmfGen<Fmpq>]) -> Vec<Vec<Fmpq>> {
    let vv: Vec<_> = forms.iter().map(|f| f.fc_vector_u_nonneg()).collect();
    let n = forms.len();
    let m = vv[0].len();
    let mut mat = FmpqMat::new(m as i64, n as i64);
    for (i, v) in vv.iter().enumerate() {
        for (j, x) in v.iter().enumerate() {
            mat.set_entry(j as isize, i as isize, x);
        }
    }
    mat.right_kernel_basis()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weight_add_mul() {
        let a: Weight = Some((1, 2));
        let b: Weight = Some((3, 2));
        assert_eq!(weight_mul(a, b), Some((4, 4)));
        assert_eq!(weight_mul(a, None), None);
        assert_eq!(weight_add(a, b), None);
        assert_eq!(weight_add(a, a), a);
    }
}
