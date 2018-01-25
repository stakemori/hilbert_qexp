use elements::UBounds;
use std::ops::SubAssign;
use bignum::BigNumber;

pub fn sub_assign<T>(f_vec: &mut Vec<T>, g_vec: &Vec<T>, v: usize, u_bds: &UBounds)
where
    for<'a> T: SubAssign<&'a T>,
{
    let bd = u_bds.vec[v];
    for i in (0..(bd + 1)).filter(|x| is_even!(x + v)) {
        T::sub_assign(&mut f_vec[bd + i], &g_vec[bd + i]);
    }
    for i in (1..(bd + 1)).filter(|x| is_even!(x + v)) {
        let idx = bd - i;
        T::sub_assign(&mut f_vec[idx], &g_vec[idx])
    }
}

/// set v_g + v_h coefficient (Laurant polynomial of e(u)) to the product of
/// g[v_g], h[v_h].
/// g_vec[gap_g + i] is defined only when i.abs() <= bd_g.
pub fn mul_mut<T>(
    f_vec: &mut Vec<T>,
    g_vec: &Vec<T>,
    h_vec: &Vec<T>,
    v_g: usize,
    v_h: usize,
    gap_g: usize,
    gap_h: usize,
    gap_gh: usize,
    u_bds: &UBounds,
    parity_g: usize,
    parity_h: usize,
    parity_gh: usize,
) where
    T: BigNumber,
{
    let bd_g = u_bds.vec[v_g];
    let bd_h = u_bds.vec[v_h];
    let bd_gh = u_bds.vec[v_g + v_h];

    for i in (0..(bd_gh + 1)).filter(|&x| is_even!(v_g + v_h + x + parity_gh)) {
        f_vec[(gap_gh + i) as usize].set_ui_g(0);
        f_vec[(gap_gh - i) as usize].set_ui_g(0);
    }
    let mut tmp = T::R::default();
    // naive implementation of polynomial multiplication
    // i -> i - bd_g
    for i in (0..(2 * bd_g + 1)).filter(|&x| is_even!(v_g + x + bd_g + parity_g)) {
        for j in (0..(2 * bd_h + 1)).filter(|&x| is_even!(v_h + x + bd_h + parity_h)) {
            f_vec[gap_gh + i + j - bd_g - bd_h].addmul_mut_g(
                &g_vec[gap_g + i - bd_g],
                &h_vec[gap_h + j - bd_h],
                &mut tmp,
            );
        }
    }
}

/// Assuming f_vec is divisible by g_vec in integral coefficients, set h_vec to
/// f_vec/g_vec.
pub fn div_mut<T>(
    f_vec: &Vec<T>,
    g_vec: &Vec<T>,
    h_vec: &mut Vec<T>,
    v_g: usize,
    v_h: usize,
    gap_g: usize,
    gap_h: usize,
    gap_gh: usize,
    u_bds: &UBounds,
) where
    T: BigNumber,
    for<'a> T: SubAssign<&'a T>,
{
    let gap_gh = gap_gh as i64;
    let gap_g = gap_g as i64;
    let gap_h = gap_h as i64;
    let bd_g = u_bds.vec[v_g] as i64;
    let bd_h = u_bds.vec[v_h] as i64;
    for i in -bd_h..(bd_h + 1) {
        h_vec[(gap_h + i) as usize].set_ui_g(0);
    }
    let mut tmp = T::R::default();
    // initial index of g
    let n = ((-bd_g)..(bd_g + 1))
        .rev()
        .filter(|&i| !g_vec[(gap_g + i) as usize].is_zero_g())
        .next()
        .unwrap();

    let mut tmp_elt = T::new_g();
    h_vec[(gap_h + bd_h) as usize].set_g(&f_vec[(gap_gh + n + bd_h) as usize]);
    h_vec[(gap_h + bd_h) as usize].set_divexact_g(&g_vec[(gap_g + n) as usize], &mut tmp);

    for m in ((-bd_h)..bd_h).rev() {
        tmp_elt.set_g(&f_vec[(gap_gh + n + m) as usize]);
        for i in ((n + m - bd_h)..n).filter(|&x| -bd_g <= x) {
            tmp_elt.submul_mut_g(
                &g_vec[(gap_g + i) as usize],
                &h_vec[(gap_h + n + m - i) as usize],
                &mut tmp,
            );

        }
        h_vec[(gap_h + m) as usize].set_g(&tmp_elt);
        h_vec[(gap_h + m) as usize].set_divexact_g(&g_vec[(gap_g + n) as usize], &mut tmp);
    }
}
