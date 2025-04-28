use std::{
    f64::consts::SQRT_2,
    fmt::{Debug, Display},
    fs::read_to_string,
    iter::zip,
    path::Path,
    str::FromStr,
};

use nalgebra::dmatrix;
use tensor::Tensor4;
type Tensor3 = tensor::tensor3::Tensor3<f64>;

use crate::{
    consts::{FACT3, FACT4, FUNIT3, FUNIT4, ICTOP, IPTOC, WAVE},
    f3qcm::F3qcm,
    f4qcm::F4qcm,
    ifrm1::Ifrm1,
    mode::Mode,
    state::State,
    Dmat, Dvec, Spectro,
};

// separate for macro
use crate::f4qcm;

use self::linalg::symm_eigen_decomp;

pub mod linalg;

impl Display for Spectro {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use crate::Curvil::*;
        writeln!(f, "# SPECTRO #############")?;
        for chunk in self.header.chunks(15) {
            for i in chunk {
                write!(f, "{i:5}")?;
            }
            writeln!(f)?;
        }
        writeln!(f, "# GEOM #############")?;
        writeln!(f, "{:5}{:5}", self.geom.atoms.len() + self.dummies.len(), 0)?;
        for atom in &self.geom.atoms {
            writeln!(
                f,
                "{:5.2}{:16.8}{:16.8}{:16.8}",
                atom.atomic_number as f64, atom.x, atom.y, atom.z
            )?;
        }
        for dummy in &self.dummies {
            let atom = dummy.get_vals(&self.geom);
            writeln!(
                f,
                "{:5.2}{:16.8}{:16.8}{:16.8}",
                0.0, atom[0], atom[1], atom[2],
            )?;
        }
        writeln!(f, "# WEIGHT #############")?;
        writeln!(f, "{:5}", self.weights.len())?;
        for weight in &self.weights {
            writeln!(f, "{:5}{:12.6}", weight.0, weight.1)?;
        }
        writeln!(f, "# CURVIL #############")?;
        for curvil in &self.curvils {
            match curvil {
                Bond(i, j) => write!(f, "{i:5}{j:5}")?,
                Bend(i, j, k) => write!(f, "{i:5}{j:5}{k:5}")?,
                Tors(i, j, k, l) => write!(f, "{i:5}{j:5}{k:5}{l:5}")?,
            }
            writeln!(f)?;
        }
        if !self.degmodes.is_empty() {
            writeln!(f, "# DEGMODE #############")?;
            for curvil in &self.degmodes {
                for i in curvil {
                    write!(f, "{i:5}")?;
                }
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

/// parse an entire `line` into a vector of the same type
pub(crate) fn parse_line<T: FromStr>(line: &str) -> Vec<T>
where
    <T as FromStr>::Err: Debug,
{
    line.split_whitespace()
        .map(|s| {
            s.parse::<T>()
                .unwrap_or_else(|_| panic!("failed to parse {line}"))
        })
        .collect::<Vec<_>>()
}

/// helper for sorting the find3r and find4t indices
fn sort_indices<const N: usize>(mut indices: [usize; N]) -> [usize; N] {
    indices.sort();
    indices
}

/// cubic force constant indexing formula. I think it relies on the fortran
/// numbering though, so I need to add one initially and then subtract one at
/// the end. this returns *indices* so it doesn't work directly for computing
/// lengths. For example, `find3r(3, 3, 3) = 19` and `find3r(2, 2, 2) = 9`, so
/// to get the length of the cubic force constant vector for water you have to
/// do `find3r(2, 2, 2) + 1`, where 2 is `nvib-1`
#[inline]
pub(crate) fn find3(i: usize, j: usize, k: usize) -> usize {
    let [i, j, k] = sort_indices([i + 1, j + 1, k + 1]);
    i + (j - 1) * j / 2 + (k - 1) * k * (k + 1) / 6 - 1
}

/// quartic force constant indexing formula. it relies on the fortran numbering,
/// so I need to add one initially and then subtract one at the end
#[inline]
pub(crate) fn find4(i: usize, j: usize, k: usize, l: usize) -> usize {
    let [i, j, k, l] = sort_indices([i + 1, j + 1, k + 1, l + 1]);
    i + (j - 1) * j / 2
        + (k - 1) * k * (k + 1) / 6
        + (l - 1) * l * (l + 1) * (l + 2) / 24
        - 1
}

#[inline]
pub fn ioff(n: usize) -> usize {
    ((n - 1) * n) / 2
}

/// convert harmonic frequencies to wavenumbers. not sure if this works on
/// anything else
pub fn to_wavenumbers(freqs: &Dvec) -> Dvec {
    Dvec::from_iterator(
        freqs.len(),
        freqs.iter().map(|f| {
            if *f < 0.0 {
                -WAVE * f64::sqrt(-f)
            } else {
                WAVE * f64::sqrt(*f)
            }
        }),
    )
}

/// Load a whitespace-separated vector from `infile`.
///
/// Panics if there is an error reading the file or parsing any of the entries
/// as `f64`.
pub fn load_vec<P>(infile: P) -> Vec<f64>
where
    P: AsRef<Path> + std::fmt::Debug,
{
    let data = match read_to_string(&infile) {
        Ok(data) => data,
        Err(e) => panic!("failed to read {infile:?} with {e}"),
    };
    data.split_ascii_whitespace()
        .map(|s| s.parse().unwrap())
        .collect()
}

pub(crate) fn close(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

/// transform the cubic force constants in `f3x` to normal coordinates
pub fn force3(
    n3n: usize,
    mut f3x: Tensor3,
    lx: &Dmat,
    nvib: usize,
    freq: &Dvec,
) -> F3qcm {
    let start = (0, 0);
    let end = (n3n, n3n - 1);
    for kabc in 0..n3n {
        let mut dd =
            Dmat::from_row_slice(n3n, n3n, f3x.submatrix(start, end, kabc));
        dd *= FUNIT3;
        let ee = lx.transpose() * dd * lx;
        f3x.set_submatrix(start, end, kabc, ee.data.as_slice());
    }
    let mut f3qcm =
        F3qcm::with_capacity(find3(nvib - 1, nvib - 1, nvib - 1) + 1);
    for i in 0..nvib {
        let wi = freq[i];
        for j in 0..=i {
            let wj = freq[j];
            for k in 0..=j {
                let wk = freq[k];
                let wijk = wi * wj * wk;
                let fact = FACT3 / wijk.sqrt();
                let mut val = 0.0;
                for l in 0..n3n {
                    val += f3x[(i, j, l)] * lx[(l, k)];
                }
                f3qcm.push(val * fact);
            }
        }
    }
    f3qcm
}

pub fn force4(
    n3n: usize,
    f4x: &Tensor4,
    lx: &Dmat,
    nvib: usize,
    harms: &Dvec,
) -> F4qcm {
    let lxt = lx.transpose();
    let mut f4q = Tensor4::zeros(n3n, n3n, n3n, n3n);
    for kabc in 0..n3n {
        for labc in 0..n3n {
            let mut dd = Dmat::from_row_slice(
                n3n,
                n3n,
                f4x.submatrix((0, 0), (n3n, n3n - 1), kabc, labc),
            );
            dd *= FUNIT4;
            let ee = lxt.clone() * dd * lx.clone();
            // technically this as_slice call is wrong because it comes out in
            // column-major order, but the matrix is symmetric
            f4q.set_submatrix(
                (0, 0),
                (n3n, n3n - 1),
                kabc,
                labc,
                ee.data.as_slice(),
            );
        }
    }
    // now can I include the loop above in here - the holy grail
    let n = nvib - 1;
    let mut f4qcm = f4qcm![0.0; find4(n, n, n, n) + 1];
    for i in 0..nvib {
        let wi = harms[i];
        for j in 0..=i {
            let wj = harms[j];
            let mut dd = Dmat::zeros(n3n, n3n);
            for kabc in 0..n3n {
                for labc in 0..n3n {
                    dd[(kabc, labc)] = f4q[(i, j, kabc, labc)];
                }
            }
            let ee = &lxt * dd * lx;
            for k in 0..nvib {
                let wk = harms[k];
                for l in 0..=k {
                    let wl = harms[l];
                    let wijkl = wi * wj * wk * wl;
                    let sqws = wijkl.sqrt();
                    let fact = FACT4 / sqws;
                    let ijkl = (k, l, i, j);
                    f4qcm[ijkl] = ee[(k, l)] * fact;
                }
            }
        }
    }
    f4qcm
}

/// make E0 for asymmetric tops or the first component of E0 for symmetric tops
pub(crate) fn make_e0(
    modes: &[Mode],
    f4qcm: &F4qcm,
    f3qcm: &F3qcm,
    freq: &Dvec,
    ifrm1: &Ifrm1,
    ifrmchk: &tensor::Tensor3<usize>,
) -> f64 {
    // NOTE: took out some weird IA stuff here and reproduced their results.
    // maybe my signs are actually right and theirs are wrong.
    let (i1mode, _, _) = Mode::partition(modes);
    let f3k: f64 = i1mode
        .iter()
        .map(|&k| -7.0 * f3qcm[(k, k, k)].powi(2) / (576.0 * freq[k]))
        .sum();
    let f4k: f64 = i1mode.iter().map(|&k| f4qcm[(k, k, k, k)] / 64.0).sum();
    let mut f3kkl = 0.0;
    for &k in &i1mode {
        // kkkk and kkk terms
        let wk = freq[k].powi(2);
        // kkl terms
        for &l in i1mode.iter().filter(|&&l| k != l) {
            let wl = freq[l].powi(2);
            let zval1 = f3qcm[(k, k, l)].powi(2);
            let res = ifrm1.get(&k);
            if res.is_some() && *res.unwrap() == l {
                let delta = 2.0 * freq[k] + freq[l];
                f3kkl += 3.0 * zval1 / (128.0 * delta);
            } else {
                let delta = 4.0 * wk - wl;
                f3kkl += 3.0 * zval1 * freq[l] / (64.0 * delta);
            }
        }
    }
    // klm terms
    let mut f3klm = 0.0;
    for &k in &i1mode {
        for &l in i1mode.iter().filter(|&&l| k > l) {
            for &m in i1mode.iter().filter(|&&m| l > m) {
                let zval3 = f3qcm[(k, l, m)].powi(2);
                let xklm = freq[k] * freq[l] * freq[m];
                let d1 = freq[k] + freq[l] + freq[m];
                let d2 = freq[k] - freq[l] + freq[m];
                let d3 = freq[k] + freq[l] - freq[m];
                let d4 = freq[k] - freq[l] - freq[m];
                if ifrmchk[(k, l, m)] != 0 {
                    f3klm -= zval3 * (1.0 / d1 + 1.0 / d2 + 1.0 / d4) / 16.0;
                } else if ifrmchk[(l, m, k)] != 0 {
                    f3klm -= zval3 * (1.0 / d1 + 1.0 / d2 - 1.0 / d3) / 16.0;
                } else if ifrmchk[(k, m, l)] != 0 {
                    f3klm -= zval3 * (1.0 / d1 - 1.0 / d3 + 1.0 / d4) / 16.0;
                } else {
                    f3klm -= zval3 * xklm / (4.0 * d1 * d2 * d3 * d4);
                }
            }
        }
    }
    // biggest differences in f4k and f3klm, but I think it's okay
    f4k + f3k + f3kkl + f3klm
}

/// compute the fundamental frequencies from the harmonic frequencies and the
/// anharmonic constants
pub fn make_funds(freq: &Dvec, nvib: usize, xcnst: &Dmat) -> Vec<f64> {
    let mut fund = Vec::with_capacity(freq.len());
    for i in 0..nvib {
        let mut val = freq[i] + 2.0 * xcnst[(i, i)];
        for j in 0..nvib {
            if j != i {
                val += 0.5 * xcnst[(i, j)];
            }
        }
        fund.push(val);
    }
    fund
}

/// take a vec of energy: state pairs and print them in SPECTRO's format
pub(crate) fn print_vib_states(reng: &[f64], i1sts: &Vec<State>) {
    println!(
        "{:^10}{:^20}{:^20}{:>21}",
        "STATE NO.", "ENERGY (CM-1)", "ABOVE ZPT", "VIBRATIONAL STATE"
    );
    for (i, (energy, state)) in zip(reng, i1sts).enumerate() {
        print!("{:5}{:20.4}{:20.4}", i + 1, energy, *energy - reng[0]);
        print!("{:>21}", "NON-DEG (Vs) :");
        print!("{state:5}");
        println!();
    }
}

pub(crate) fn rsfrm2(
    ijst: usize,
    kst: usize,
    ivib: usize,
    jvib: usize,
    kvib: usize,
    f3qcm: &F3qcm,
    eng: &mut [f64],
) {
    let val = f3qcm[(ivib, jvib, kvib)] / (2.0 * SQRT_2);
    let eres = dmatrix![
    eng[ijst] - eng[0], val;
    val, eng[kst] - eng[0];
    ];
    // TODO left out error measures
    let (eigval, eigvec) = symm_eigen_decomp(eres, false);
    let a = eigvec[(0, 0)];
    let b = eigvec[(1, 0)];
    if a.abs() > b.abs() {
        eng[ijst] = eigval[0] + eng[0];
        eng[kst] = eigval[1] + eng[0];
    } else {
        eng[kst] = eigval[0] + eng[0];
        eng[ijst] = eigval[1] + eng[0];
    }
    // TODO left out properties
}

/// calculate the type-1 fermi resonance contribution to the energy. `deg` is
/// usually false, but true in one case for symmetric tops
pub(crate) fn rsfrm1(
    ist: usize,
    jst: usize,
    ivib: usize,
    jvib: usize,
    f3qcm: &F3qcm,
    eng: &mut [f64],
    deg: bool,
) {
    let val = if deg {
        f3qcm[(ivib, ivib, jvib)] / (2.0 * SQRT_2)
    } else {
        0.25 * f3qcm[(ivib, ivib, jvib)]
    };
    // TODO left out printed error measures
    let eres = dmatrix![
    eng[ist] - eng[0], val;
    val, eng[jst] - eng[0];
    ];
    let (eigval, eigvec) = linalg::symm_eigen_decomp(eres, false);
    let a = eigvec[(0, 0)];
    let b = eigvec[(1, 0)];
    if a.abs() > b.abs() {
        eng[ist] = eigval[0] + eng[0];
        eng[jst] = eigval[1] + eng[0];
    } else {
        eng[ist] = eigval[1] + eng[0];
        eng[jst] = eigval[0] + eng[0];
    }
    // TODO left out the calculation of updated properties
}

/// convert tau to tau prime in wavenumbers
pub(crate) fn tau_prime(maxcor: usize, tau: &Tensor4) -> Dmat {
    let mut taupcm = Dmat::zeros(maxcor, maxcor);
    for i in 0..maxcor {
        for j in 0..maxcor {
            taupcm[(i, j)] = tau[(i, i, j, j)];
            if i != j {
                taupcm[(i, j)] += 2.0 * tau[(i, j, i, j)];
            }
        }
    }
    taupcm
}

pub(crate) fn make_tau(
    maxcor: usize,
    nvib: usize,
    freq: &Dvec,
    primat: &[f64],
    wila: &Dmat,
) -> Tensor4 {
    // convert to cm-1 from the biggest mess you've ever seen
    const CONST1: f64 = 3.833384078e04;
    let mut tau = Tensor4::zeros(maxcor, maxcor, maxcor, maxcor);
    for i in 0..maxcor {
        for j in 0..maxcor {
            for k in 0..maxcor {
                for l in 0..maxcor {
                    let x = ioff(i.max(j) + 1) + i.min(j);
                    let y = ioff(k.max(l) + 1) + k.min(l);
                    let mut sum = 0.0;
                    for v in 0..nvib {
                        let div = freq[v].powi(2)
                            * primat[i]
                            * primat[j]
                            * primat[k]
                            * primat[l];
                        sum += wila[(v, x)] * wila[(v, y)] / div;
                    }
                    tau[(i, j, k, l)] = -0.5 * CONST1 * sum;
                }
            }
        }
    }
    tau
}

/// set up vectors for principal -> cartesian and cartesian -> principal
/// transformations for asymmetric tops
pub(crate) fn princ_cart(irep: usize) -> ([usize; 3], [usize; 3]) {
    let ic = [IPTOC[(0, irep)], IPTOC[(1, irep)], IPTOC[(2, irep)]];
    let id = [ICTOP[(0, irep)], ICTOP[(1, irep)], ICTOP[(2, irep)]];
    (ic, id)
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::assert_abs_diff_eq;
    use nalgebra::dmatrix;

    use crate::{consts::FACT2, load_fc2};

    #[test]
    fn test_find3r() {
        let got = find3(2, 2, 2);
        let want = 9;
        assert_eq!(got, want);
    }

    #[test]
    fn test_taupcm() {
        let s = Spectro::load("testfiles/h2o/spectro.in");
        let fc2 = load_fc2("testfiles/fort.15", s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (harms, lxm) = linalg::symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let (_zmat, wila) = s.zeta(&lxm, &w);
        let tau = make_tau(3, 3, &freq, &s.primat, &wila);
        let got = tau_prime(3, &tau);
        let want = dmatrix![
        -0.08628870,  0.01018052, -0.00283749;
         0.01018052, -0.00839612, -0.00138895;
        -0.00283749, -0.00138895, -0.00093412;
           ];
        assert_abs_diff_eq!(got, want, epsilon = 2e-7);
    }
}
