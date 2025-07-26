use crate::{
    DMat, Siic, Tensor3, foreach,
    geom::Geom,
    hmat::{Hmat, hijs1},
    htens4::{Htens4, h4th1, h4th2},
};

use self::utils::fill3b;

mod display;

pub(crate) mod utils {
    use crate::Tensor3;

    pub(crate) fn fill3a(t: &mut Tensor3, nsx: usize) {
        for p in 0..nsx {
            for n in 0..p {
                for m in 0..n {
                    t[(n, m, p)] = t[(m, n, p)];
                    t[(n, p, m)] = t[(m, n, p)];
                    t[(m, p, n)] = t[(m, n, p)];
                    t[(p, m, n)] = t[(m, n, p)];
                    t[(p, n, m)] = t[(m, n, p)];
                }
                t[(n, p, n)] = t[(n, n, p)];
                t[(p, n, n)] = t[(n, n, p)];
            }
            for m in 0..p {
                t[(p, m, p)] = t[(m, p, p)];
                t[(p, p, m)] = t[(m, p, p)];
            }
        }
    }

    /// copy values across the 3D diagonals
    pub(crate) fn fill3b(t: &mut Tensor3) {
        for m in 0..3 {
            for n in 0..m {
                for p in 0..n {
                    t[(n, m, p)] = t[(m, n, p)];
                    t[(n, p, m)] = t[(m, n, p)];
                    t[(m, p, n)] = t[(m, n, p)];
                    t[(p, m, n)] = t[(m, n, p)];
                    t[(p, n, m)] = t[(m, n, p)];
                }
                t[(n, m, n)] = t[(m, n, n)];
                t[(n, n, m)] = t[(m, n, n)];
            }
            for p in 0..m {
                t[(m, p, m)] = t[(m, m, p)];
                t[(p, m, m)] = t[(m, m, p)];
            }
        }
    }
}

pub struct Htens {
    pub h111: Tensor3,
    pub h112: Tensor3,
    pub h113: Tensor3,
    pub h123: Tensor3,
    pub h221: Tensor3,
    pub h222: Tensor3,
    pub h223: Tensor3,
    pub h331: Tensor3,
    pub h332: Tensor3,
    pub h333: Tensor3,
    pub h411: Tensor3,
    pub h421: Tensor3,
    pub h422: Tensor3,
    pub h431: Tensor3,
    pub h432: Tensor3,
    pub h433: Tensor3,
    pub h441: Tensor3,
    pub h442: Tensor3,
    pub h443: Tensor3,
    pub h444: Tensor3,
}

/// unpack (splat) an s vector into a series of slices
#[macro_export]
macro_rules! splat {
    ($s:expr, $($var:ident => $idx:expr$(,)*)*) => {
	$(
	    let $var = &$s[3*$idx..3*$idx+3];
	)*
    };
}

/// unpack (splat) an s vector into a series of DVecs
#[macro_export]
macro_rules! dsplat {
    ($s:expr, $($var:ident => $idx:expr$(,)*)*) => {
	$(
	    let $var = nalgebra::DVector::from_row_slice(&$s[3*$idx..3*$idx+3]);
	)*
    };
}

/// helper function for calling Htens::new with an [Siic::Stretch]
pub fn hijks1(geom: &Geom, k1: usize, k2: usize) -> Tensor3 {
    Htens::new(geom, &Siic::Stretch(k1, k2)).h111
}

/// helper function for calling Htens::new with an [Siic::Bend]
pub fn hijks2(geom: &Geom, i: usize, j: usize, k: usize) -> Htens {
    Htens::new(geom, &Siic::Bend(i, j, k))
}

impl Htens {
    pub fn zeros() -> Self {
        Self {
            h111: Tensor3::zeros((3, 3, 3)),
            h112: Tensor3::zeros((3, 3, 3)),
            h113: Tensor3::zeros((3, 3, 3)),
            h123: Tensor3::zeros((3, 3, 3)),
            h221: Tensor3::zeros((3, 3, 3)),
            h222: Tensor3::zeros((3, 3, 3)),
            h223: Tensor3::zeros((3, 3, 3)),
            h331: Tensor3::zeros((3, 3, 3)),
            h332: Tensor3::zeros((3, 3, 3)),
            h333: Tensor3::zeros((3, 3, 3)),
            h411: Tensor3::zeros((3, 3, 3)),
            h421: Tensor3::zeros((3, 3, 3)),
            h422: Tensor3::zeros((3, 3, 3)),
            h431: Tensor3::zeros((3, 3, 3)),
            h432: Tensor3::zeros((3, 3, 3)),
            h433: Tensor3::zeros((3, 3, 3)),
            h441: Tensor3::zeros((3, 3, 3)),
            h442: Tensor3::zeros((3, 3, 3)),
            h443: Tensor3::zeros((3, 3, 3)),
            h444: Tensor3::zeros((3, 3, 3)),
        }
    }

    pub fn new(geom: &Geom, siic: &Siic) -> Self {
        use Siic::*;
        let mut h = Htens::zeros();
        match siic {
            // HIJKS1
            Stretch(j, i) => {
                Self::stretch(geom, siic, i, j, &mut h);
            }
            // HIJKS2
            Bend(i, j, k) => {
                Self::bend(geom, siic, i, j, k, &mut h);
            }
            // HIJKS6
            Torsion(i, j, k, l) => {
                Self::torsion(geom, siic, i, j, k, l, &mut h)
            }
            // HIJKS3
            Lin1(i, j, k, _) => Self::lin1(geom, siic, i, j, k, &mut h),
            // hijks7
            Out(i, j, k, l) => Self::out(geom, siic, i, j, k, l, &mut h),
            // hijks8
            Linx(i, j, k, l) => Self::linx(geom, siic, i, j, k, l, &mut h),
            // hijks9
            Liny(i, j, k, l) => Self::liny(geom, i, j, k, l, &mut h),
        }
        h
    }

    /// the name suggests some kind of triple product, but I'm not really sure.
    /// I think it's a Cartesian product of some kind maybe
    fn tripro() -> Tensor3 {
        let mut ret = Tensor3::zeros((3, 3, 3));
        for k in 0..3 {
            let mut vect = nalgebra::vector![0.0, 0.0, 0.0];
            vect[k] = 1.0;
            let rmat = Hmat::mat1(&vect);
            for j in 0..3 {
                for i in 0..3 {
                    ret[(i, j, k)] = rmat[(i, j)];
                }
            }
        }
        ret
    }

    fn stretch(geom: &Geom, siic: &Siic, a: &usize, b: &usize, h: &mut Htens) {
        let v1 = geom.unit(*a, *b);
        let w1 = 1.0 / geom.dist(*a, *b);
        let hm = Hmat::new(geom, siic);
        for k in 0..3 {
            for j in k..3 {
                for i in j..3 {
                    h.h111[(i, j, k)] = -(v1[i] * hm.h11[(k, j)]
                        + v1[j] * hm.h11[(k, i)]
                        + v1[k] * hm.h11[(j, i)])
                        * w1;
                }
            }
        }
        fill3b(&mut h.h111);
    }

    fn bend(
        geom: &Geom,
        siic: &Siic,
        a: &usize,
        b: &usize,
        c: &usize,
        h: &mut Htens,
    ) {
        use Siic::*;
        // copied from h_mat Bend
        let tmp = geom.s_vec(siic);
        splat!(tmp, v1 => a, v3 => c);
        let e21 = geom.unit(*b, *a);
        let e23 = geom.unit(*b, *c);
        let t21 = geom.dist(*b, *a);
        let t23 = geom.dist(*b, *c);
        let h11a = Hmat::new(geom, &Stretch(*a, *b)).h11;
        let h33a = Hmat::new(geom, &Stretch(*c, *b)).h11;
        let phi = geom.angle(*a, *b, *c);
        // end copy
        let hijs2 = Hmat::new(geom, siic);
        let h111a = Self::new(geom, &Stretch(*a, *b)).h111;
        let h333a = Self::new(geom, &Stretch(*c, *b)).h111;
        let sphi = phi.sin();
        let ctphi = phi.cos() / sphi;
        let w1 = 1.0 / t21;
        let w2 = 1.0 / t23;
        let w3 = ctphi * w1;
        let w4 = ctphi * w2;
        for k in 0..3 {
            let w5 = v1[k] * ctphi + e21[k] * w1;
            let w6 = e21[k] * w3;
            let w7 = v1[k] * w1;
            let w8 = v3[k] * ctphi + e23[k] * w2;
            let w9 = e23[k] * w4;
            let w10 = v3[k] * w2;
            for j in 0..3 {
                for i in 0..3 {
                    h.h221[(i, j, k)] = w5 * hijs2.h11[(i, j)]
                        + v1[i] * v1[j] * w6
                        + h11a[(i, j)] * w7;
                    h.h223[(i, j, k)] = w8 * hijs2.h33[(i, j)]
                        + v3[i] * v3[j] * w9
                        + h33a[(i, j)] * w10;
                }
            }
        }

        h.h111 = -(h.h221.clone()
            + h.h221.view().permuted_axes((2, 0, 1))
            + h.h221.view().permuted_axes((0, 2, 1)))
            + w3 * h111a;
        h.h333 = -(h.h223.clone()
            + h.h223.view().permuted_axes((2, 0, 1))
            + h.h223.view().permuted_axes((0, 2, 1)))
            + h333a * w4;

        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    h.h111[(i, j, k)] += v1[i] * v1[j] * v1[k];
                    h.h333[(i, j, k)] += v3[i] * v3[j] * v3[k];
                }
            }
        }

        for i in 0..3 {
            let w3 = v1[i] * ctphi + e21[i] * w1;
            let w4 = v3[i] * ctphi + e23[i] * w2;
            for j in 0..3 {
                for k in 0..3 {
                    h.h221[(i, j, k)] = w3 * hijs2.h31[(k, j)];
                    h.h223[(i, j, k)] = w4 * hijs2.h31[(j, k)];
                }
            }
        }

        h.h113 = -(h.h221.clone() + h.h221.view().permuted_axes((1, 0, 2)));
        h.h331 = -(h.h223.clone() + h.h223.view().permuted_axes((1, 0, 2)));
        let w3 = 1.0 / (sphi * sphi);
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    h.h113[(i, j, k)] +=
                        v3[k] * (v1[i] * v1[j] - h11a[(i, j)] * w1) * w3;
                    h.h331[(i, j, k)] +=
                        v1[k] * (v3[i] * v3[j] - h33a[(i, j)] * w2) * w3;
                }
            }
        }

        h.h123 = -(h.h331.clone().permuted_axes((2, 0, 1)) + &h.h113);
        h.h112 = -(&h.h111 + &h.h113);
        h.h332 = -(&h.h333 + &h.h331);
        h.h221 = -(h.h123.clone().permuted_axes((1, 2, 0))
            + h.h112.view().permuted_axes((2, 1, 0)));
        h.h223 = -(h.h332.clone().permuted_axes((1, 2, 0)) + &h.h123);
        h.h222 = -(h.h223.clone().permuted_axes((2, 0, 1))
            + h.h221.view().permuted_axes((2, 0, 1)));
    }

    fn torsion(
        geom: &Geom,
        siic: &Siic,
        a: &usize,
        b: &usize,
        c: &usize,
        d: &usize,
        h: &mut Htens,
    ) {
        use Siic::*;
        let tmp = geom.s_vec(siic);
        splat!(tmp, v1 => a, v4 => d);

        // unit and non-unit vectors
        let e21 = geom.unit(*b, *a);
        let e23 = geom.unit(*b, *c);
        let e34 = geom.unit(*c, *d);
        let t21 = geom.dist(*b, *a);
        let t23 = geom.dist(*b, *c);
        let t34 = geom.dist(*c, *d);

        // angles
        let p2 = geom.angle(*a, *b, *c);
        let tmp = geom.s_vec(&Bend(*a, *b, *c));
        dsplat!(tmp, bp21 => a, bp22 => b, bp23 => c);

        let p3 = geom.angle(*b, *c, *d);
        let tmp = geom.s_vec(&Bend(*b, *c, *d));
        dsplat!(tmp, bp32 => b, bp33 => c, bp34 => d);

        let c1 = 1.0 / t21;
        let c2 = 1.0 / t34;
        let c3 = 1.0 / t23;
        let c4 = f64::sin(p2);
        let c5 = f64::cos(p2);
        let c6 = c5 / c4;
        let c7 = f64::sin(p3);
        let c8 = f64::cos(p3);
        let c9 = c8 / c7;
        let c10 = 1.0 / (c4 * c4);
        let c11 = 1.0 / (c7 * c7);
        let c12 = c1 * c1;
        let c13 = c2 * c2;
        let c14 = c3 * c3;
        let c15 = t21 * c3;
        let c16 = t34 * c3;
        let w1 = 2.0 * c10 * c12;
        let w2 = 2.0 * c11 * c13;
        let w5 = c10 * c12;
        let w6 = c11 * c13;
        let w7 = w5 * c3;
        let w8 = w6 * c3;
        let w11 = c1 * c3 * c10;
        let w12 = c2 * c3 * c11;
        let w15 = 2.0 * c1;
        let w16 = 2.0 * c12;
        let w17 = 2.0 * c6;
        let w18 = 2.0 * c10;
        let w19 = 2.0 * c3;
        let w20 = c4 * c3;
        let w21 = c4 * c15;
        let w22 = c5 * c15;
        let w23 = w22 * c3;
        let w24 = c3 * c15;
        let w26 = c3 * c4 * c15;
        let w27 = c5 * c14;
        let w28 = 2.0 * c2;
        let w29 = 2.0 * c13;
        let w30 = 2.0 * c9;
        let w31 = 2.0 * c11;
        let w32 = 2.0 * c3;
        let w33 = c7 * c3;
        let w34 = c7 * c16;
        let w35 = c8 * c16;
        let w36 = w35 * c3;
        let w37 = t34 * c14;
        let w38 = c3 * c7 * c16;
        let w39 = c8 * c14;
        let w42 = c5 * c15;
        let w43 = w42 - 1.0;
        let w44 = c8 * c16;
        let w45 = w44 - 1.0;
        let w48 = c5 * c3;
        let w49 = c4 * c15;
        let w50 = c5 * t21 * c14;
        let w51 = c8 * c3;
        let w52 = c7 * c16;
        let w53 = c8 * t34 * c14;
        let w54 = 2.0 * c6;
        let w55 = 2.0 * c9;
        let w56 = 2.0 * c1;
        let w57 = 2.0 * c2;
        let w58 = c5 * c3;
        let w59 = c8 * c3;

        // matrices
        let hm32 = Hmat::mat1(&e23);
        let hm21 = Hmat::mat1(&e21);
        let hm43 = Hmat::mat1(&e34);
        let Hmat { h11, h31, h32, .. } = Hmat::new(geom, &Bend(*a, *b, *c));
        let h44 = Hmat::new(geom, &Stretch(*a, *b)).h11;
        let h42 = Hmat::new(geom, &Stretch(*b, *c)).h11;
        let h43 = 2.0 * (w15 * h44 + c6 * h11 - c10 * &bp21 * bp21.transpose())
            - w16 * DMat::identity(3, 3);
        let h43a = w17 * h31.transpose() - w18 * &bp21 * bp23.transpose();
        let h43b = w19 * &h42 - w17 * h32.transpose()
            + w18 * &bp22 * bp23.transpose()
            - c14 * DMat::identity(3, 3);
        let h43c = -e21 * bp23.transpose() * w20
            + h32.transpose() * w21
            + &bp22 * bp23.transpose() * w22
            - &h42 * w23;
        let h43d = &h43c
            + -w26 * &bp22 * e23.transpose()
            + w27 * (c15 * e23 - e21) * e23.transpose();
        let Hmat {
            h21: h32a,
            h31: h42a,
            h33: h44a,
            ..
        } = Hmat::new(geom, &Bend(*b, *c, *d));
        let h11a = Hmat::new(geom, &Stretch(*d, *c)).h11;
        let h31a = Hmat::new(geom, &Stretch(*c, *b)).h11;
        let h21 = 2.0
            * (w28 * h11a + c9 * h44a - c11 * &bp34 * bp34.transpose())
            - w29 * DMat::identity(3, 3);
        let h21a = w30 * h42a - w31 * &bp34 * bp32.transpose();
        let h21b = w32 * &h31a - w30 * &h32a + w31 * &bp33 * bp32.transpose()
            - c14 * DMat::identity(3, 3);
        let h21c = -e34 * bp32.transpose() * w33
            + &h32a * w34
            + &bp33 * bp32.transpose() * w35
            - &h31a * w36;
        let h21d = &h21c
            + w38 * &bp33 * e23.transpose()
            + w39 * (e34 + c16 * e23) * e23.transpose();
        let Hmat {
            h11: h11b,
            h21: h21e,
            h31: h31b,
            h42: h42b,
            h43: h43e,
            h44: h44b,
            ..
        } = Hmat::new(geom, siic);

        // scratch tensors
        let mut h_alpha = Tensor3::zeros((3, 3, 3));
        for k in 0..3 {
            h_alpha[(0, 0, k)] = e21[k] * c1 + bp21[k] * c6;
            h_alpha[(0, 1, k)] = e34[k] * c2 + bp34[k] * c9;
            h_alpha[(0, 2, k)] = e23[k] * c3 + bp23[k] * w54;
            h_alpha[(1, 0, k)] = -e23[k] * c3 + bp32[k] * w55;
            h_alpha[(1, 1, k)] = e21[k] * w56 + e23[k] * c3 - bp22[k] * w54;
            h_alpha[(1, 2, k)] = e34[k] * w57 - e23[k] * c3 - bp33[k] * w55;
            h_alpha[(2, 0, k)] = e23[k] * w58 + bp23[k] * c4;
            h_alpha[(2, 1, k)] = -e23[k] * w59 + bp32[k] * c7;
        } // end 10

        let mut h_beta = Tensor3::zeros((3, 3, 3));
        for k in 0..3 {
            h_beta[(0, 0, k)] = w52 * bp33[k] + w53 * e23[k] + w51 * e34[k];
            h_beta[(0, 1, k)] = w49 * bp22[k] - w50 * e23[k] + w48 * e21[k];
            h_beta[(0, 2, k)] = -w48 * e21[k] - w49 * bp22[k] + w50 * e23[k];
            h_beta[(1, 0, k)] = -w51 * e34[k] - w52 * bp33[k] - w53 * e23[k];
        }

        let p = Self::tripro();

        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    h.h111[(i, j, k)] = w1 * h_alpha[(0, 0, k)] * hm32[(i, j)]
                        - v1[j] * h43[(i, k)]
                        - 2.0 * h11b[(j, k)] * h_alpha[(0, 0, i)];
                    h.h444[(i, j, k)] = -w2 * h_alpha[(0, 1, k)] * hm32[(i, j)]
                        - v4[j] * h21[(i, k)]
                        - 2.0 * h44b[(j, k)] * h_alpha[(0, 1, i)];
                    h.h113[(i, j, k)] = w5 * h_alpha[(0, 2, k)] * hm32[(i, j)]
                        - w7 * p[(i, j, k)]
                        - v1[j] * h43a[(i, k)]
                        - 2.0 * h_alpha[(0, 0, i)] * h31b[(k, j)];
                    h.h442[(i, j, k)] = -w6 * h_alpha[(1, 0, k)] * hm32[(i, j)]
                        - w8 * p[(i, j, k)]
                        - v4[j] * h21a[(i, k)]
                        - 2.0 * h_alpha[(0, 1, i)] * h42b[(j, k)];
                    h.h123[(i, k, j)] =
                        -w11 * h_alpha[(1, 1, k)] * hm21[(i, j)]
                            + w7 * p[(i, j, k)]
                            + v1[i] * h43b[(k, j)]
                            - h21e[(k, i)] * h_alpha[(0, 2, j)];
                    h.h432[(i, k, j)] =
                        -w12 * h_alpha[(1, 2, k)] * hm43[(i, j)]
                            + w8 * p[(i, j, k)]
                            + v4[i] * h21b[(k, j)]
                            - h43e[(i, k)] * h_alpha[(1, 0, j)];
                    h.h223[(i, j, k)] = -v1[j]
                        * (h43c[(i, k)] + e23[i] * w24 * h_alpha[(2, 0, k)])
                        + v4[j] * h21d[(k, i)]
                        - c15 * h21e[(i, j)] * h_alpha[(2, 0, k)]
                        + c16
                            * (h43e[(j, k)] - c3 * v4[j] * e23[k])
                            * h_alpha[(2, 1, i)]
                        + h42b[(j, i)] * h_beta[(0, 0, k)]
                        + (h31b[(k, j)] - c3 * v1[j] * e23[k])
                            * h_beta[(0, 2, i)];
                    h.h332[(i, j, k)] = v1[j] * h43d[(k, i)]
                        - v4[j]
                            * (h21c[(i, k)]
                                - e23[i] * w37 * h_alpha[(2, 1, k)])
                        - c16 * h43e[(j, i)] * h_alpha[(2, 1, k)]
                        + c15
                            * (h21e[(k, j)] + c3 * v1[j] * e23[k])
                            * h_alpha[(2, 0, i)]
                        + h31b[(i, j)] * h_beta[(0, 1, k)]
                        + (h42b[(j, k)] + c3 * v4[j] * e23[k])
                            * h_beta[(1, 0, i)];
                }
            }
        } // end 102

        // I think this part does need to be done separately because the indices
        // referenced in the other tensors may not be set
        h.h223 += &(w43 * h.h123.clone().permuted_axes((1, 0, 2))
            - w44 * h.h432.clone().permuted_axes((2, 0, 1)));
        h.h332 += &(w45 * h.h432.clone().permuted_axes((1, 0, 2))
            - w42 * h.h123.clone().permuted_axes((2, 0, 1)));

        h.h112 = -(&h.h111 + &h.h113);
        h.h443 = -(h.h444.clone().permuted_axes((1, 2, 0)) + &h.h442);

        h.h221 = -(h.h112.clone().permuted_axes((0, 2, 1))
            + h.h123.view().permuted_axes((2, 1, 0)));
        h.h433 = -(h.h443.clone().permuted_axes((1, 0, 2))
            + h.h432.view().permuted_axes((0, 2, 1)));

        h.h422 = -(h.h432.clone() + h.h442.view().permuted_axes((1, 0, 2)));
        h.h331 = -(h.h123.clone().permuted_axes((1, 2, 0))
            + h.h113.view().permuted_axes((0, 2, 1)));

        h.h222 = -(&h.h221 + &h.h223 + h.h422.view().permuted_axes((1, 2, 0)));
        h.h333 = -(&h.h331 + &h.h332 + h.h433.view().permuted_axes((1, 2, 0)));
    }

    fn lin1(
        geom: &Geom,
        siic: &Siic,
        a: &usize,
        b: &usize,
        c: &usize,
        h: &mut Htens,
    ) {
        use Siic::*;
        let tmp = geom.s_vec(siic);
        splat!(tmp, v1 => a, v3 => c);
        let th = siic.value(geom);
        let e21 = geom.unit(*b, *a);
        let e23 = geom.unit(*b, *c);
        let t21 = geom.dist(*b, *a);
        let t23 = geom.dist(*b, *c);
        let h11a = Hmat::new(geom, &Stretch(*a, *b)).h11;
        let h33a = Hmat::new(geom, &Stretch(*c, *b)).h11;
        let Hmat { h11, h31, h33, .. } = Hmat::new(geom, siic);
        let h111a = Self::new(geom, &Stretch(*a, *b)).h111;
        let h333a = Self::new(geom, &Stretch(*c, *b)).h111;

        let tanth = th.tan();
        let costh = th.cos();

        let w1 = 1.0 / t21;
        let w2 = 1.0 / t23;
        let w3 = tanth * w1;
        let w4 = tanth * w2;
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    h.h221[(i, j, k)] = h11[(i, j)]
                        * (v1[k] * tanth - e21[k] / t21)
                        + v1[k] * v1[j] * e21[i] * tanth / t21
                        - (h11a[(i, j)] * v1[k]) / t21;
                    h.h223[(i, j, k)] = h33[(i, j)]
                        * (v3[k] * tanth - e23[k] / t23)
                        + v3[k] * v3[j] * e23[i] * tanth / t23
                        - (h33a[(i, j)] * v3[k]) / t23;
                }
            }
        } // end 12

        for k in 0..3 {
            for j in k..3 {
                for i in j..3 {
                    h.h111[(i, j, k)] = (h.h221[(i, j, k)]
                        + h.h221[(j, k, i)]
                        + h.h221[(k, i, j)])
                        + v1[i] * v1[j] * v1[k]
                        - h111a[(i, j, k)] * w3;
                    h.h333[(i, j, k)] = (h.h223[(i, j, k)]
                        + h.h223[(j, k, i)]
                        + h.h223[(k, i, j)])
                        + v3[i] * v3[j] * v3[k]
                        - h333a[(i, j, k)] * w4;
                }
            }
        }
        fill3b(&mut h.h111);
        fill3b(&mut h.h333);

        for i in 0..3 {
            let w5 = v1[i] * tanth - e21[i] * w1;
            let w6 = v3[i] * tanth - e23[i] * w2;
            for j in 0..3 {
                for k in 0..3 {
                    h.h221[(i, j, k)] = w5 * h31[(k, j)];
                    h.h223[(i, j, k)] = w6 * h31[(j, k)];
                }
            }
        }

        let w5 = 1.0 / (costh * costh);
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    h.h113[(i, j, k)] =
                        v3[k] * (v1[i] * v1[j] - h11a[(i, j)] * w1) * w5
                            + h.h221[(i, j, k)]
                            + h.h221[(j, i, k)];
                    h.h331[(i, j, k)] =
                        v1[k] * (v3[i] * v3[j] - h33a[(i, j)] * w2) * w5
                            + h.h223[(i, j, k)]
                            + h.h223[(j, i, k)];
                }
            }
        }

        h.h123 = -(&h.h331.clone().permuted_axes((2, 0, 1)) + &h.h113);
        h.h112 = -(&h.h111 + &h.h113);
        h.h332 = -(&h.h333 + &h.h331);

        h.h221 = -(h.h123.clone().permuted_axes((1, 2, 0))
            + h.h112.view().permuted_axes((2, 1, 0)));
        h.h223 = -(h.h332.clone().permuted_axes((1, 2, 0)) + &h.h123);

        h.h222 = -(h.h223.clone().permuted_axes((2, 0, 1))
            + h.h221.view().permuted_axes((2, 0, 1)));
    }

    fn out(
        geom: &Geom,
        siic: &Siic,
        a: &usize,
        b: &usize,
        c: &usize,
        d: &usize,
        h: &mut Htens,
    ) {
        let e21 = geom.unit(*b, *a);
        let e23 = geom.unit(*b, *c);
        let e24 = geom.unit(*b, *d);
        let t21 = geom.dist(*b, *a);
        let t23 = geom.dist(*b, *c);
        let t24 = geom.dist(*b, *d);

        // vect2 call
        let phi = Siic::Bend(*c, *b, *d).value(geom);
        let svec = geom.s_vec(&Siic::Bend(*c, *b, *d));
        splat!(svec, bp3 => c, bp4 => d);

        // vect5 call
        let svec = geom.s_vec(siic);
        let gamma = siic.value(geom);
        splat!(svec, v1 => a, v3 => c, v4 => d);

        // hijs1 call
        let ht11 = Hmat::new(geom, &Siic::Stretch(*a, *b)).h11;

        // hijs2 call
        let Hmat {
            h11: hp33,
            h31: hp43,
            h33: hp44,
            ..
        } = Hmat::new(geom, &Siic::Bend(*c, *b, *d));

        // hijs7 call
        let Hmat {
            h11,
            h31,
            h33,
            h41,
            h43,
            h44,
            ..
        } = Hmat::new(geom, siic);

        // hijks1 call
        let ht111 = Htens::new(geom, &Siic::Stretch(*a, *b)).h111;

        // hijks2 call
        let Htens {
            h113: hp334,
            h331: hp443,
            ..
        } = Htens::new(geom, &Siic::Bend(*c, *b, *d));

        let cp21 = Hmat::mat1(&e21);
        let cp24 = Hmat::mat1(&e24);
        let cp2124 = e21.cross(&e24);
        let cg = gamma.cos();
        let tg = gamma.tan();
        let sp = phi.sin();
        let ctp = phi.cos() / sp;
        let c1 = 1.0 / t21;
        let s2g = 1.0 / (cg * cg);
        let c3 = 1.0 / t23;
        let c4 = 1.0 / t24;
        let c2p = 1.0 / (sp * sp);
        let c1111 = tg * c1;
        let c1112 = s2g * c1;
        let c3331 = t24 * c3 / sp;
        let c3332 = c3331 * tg;
        let c3333 = c3331 * s2g;
        let c3335 = c3 * c3;
        let c3334 = 2.0 * c3335;
        let c4411 = t23 * c4 / sp;
        let c4412 = c4411 * s2g;
        let c431 = c3 * c4 / (cg * sp);
        let c4311 = c431 * c1;
        let c4312 = c431 * tg;
        let c4313 = c3333 * c4;
        let c4431 = c4411 * tg;
        let c4442 = c4 * c4;
        let c4441 = 2.0 * c4442;
        for i in 0..3 {
            for j in 0..=i {
                for k in 0..=j {
                    h.h111[(i, j, k)] = s2g * v1[i] * v1[j] * v1[k]
                        + tg * h11[(i, j)] * v1[k]
                        + tg * h11[(i, k)] * v1[j]
                        + c1111 * (e21[i] * v1[j] * v1[k] - ht111[(i, j, k)])
                        - c1112 * ht11[(j, k)] * v1[i]
                        - c1 * (e21[i] * h11[(j, k)]
                            + e21[j] * h11[(i, k)]
                            + e21[k] * h11[(i, j)]
                            + ht11[(i, j)] * v1[k]
                            + ht11[(i, k)] * v1[j]);
                    h.h333[(i, j, k)] = c3331
                        * (h33[(i, j)] * bp4[k] + hp43[(k, i)] * v3[j]
                            - (v3[j] * bp4[k] + tg * hp43[(k, j)])
                                * (c3 * e23[i] + ctp * bp3[i]))
                        + c3332 * hp334[(i, j, k)]
                        + c3333 * hp43[(k, j)] * v3[i]
                        + h33[(i, k)]
                            * (tg * v3[j] - ctp * bp3[j] - c3 * e23[j])
                        + v3[k]
                            * (v3[i] * v3[j] * s2g
                                + h33[(i, j)] * tg
                                + bp3[i] * bp3[j] * c2p
                                - hp33[(i, j)] * ctp
                                + e23[i] * e23[j] * c3334);
                    h.h444[(i, j, k)] = c4411
                        * (h44[(i, j)] * bp3[k] + hp43[(i, k)] * v4[j]
                            - (v4[j] * bp3[k] + hp43[(j, k)] * tg)
                                * (e24[i] * c4 + bp4[i] * ctp))
                        + c4431 * hp443[(i, j, k)]
                        + c4412 * hp43[(j, k)] * v4[i]
                        + h44[(i, k)]
                            * (tg * v4[j] - ctp * bp4[j] - c4 * e24[j])
                        + v4[k]
                            * (v4[i] * v4[j] * s2g
                                + h44[(i, j)] * tg
                                + bp4[i] * bp4[j] * c2p
                                - hp44[(i, j)] * ctp
                                + e24[i] * e24[j] * c4441);
                    if i == j {
                        h.h333[(i, j, k)] -= v3[k] * c3335;
                        h.h444[(i, j, k)] -= v4[k] * c4442;
                    }
                }
            }
        } // end 12 loop
        fill3b(&mut h.h111);
        fill3b(&mut h.h333);
        fill3b(&mut h.h444);

        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    h.h113[(j, k, i)] = s2g * v3[i] * v1[j] * v1[k]
                        + tg * h31[(i, j)] * v1[k]
                        + tg * h31[(i, k)] * v1[j]
                        - c1 * (e21[j] * h31[(i, k)] + e21[k] * h31[(i, j)])
                        - c1112 * ht11[(j, k)] * v3[i];
                    h.h411[(i, j, k)] = s2g * v4[i] * v1[j] * v1[k]
                        + tg * h41[(i, j)] * v1[k]
                        + tg * h41[(i, k)] * v1[j]
                        - c1 * (e21[j] * h41[(i, k)] + e21[k] * h41[(i, j)])
                        - c1112 * ht11[(j, k)] * v4[i];
                    h.h433[(i, j, k)] = c3331
                        * (h43[(i, j)] * bp4[k]
                            + hp44[(i, k)] * v3[j]
                            + (c4 * e24[i] - ctp * bp4[i])
                                * (v3[j] * bp4[k] + tg * hp43[(k, j)]))
                        + c3332 * hp443[(i, k, j)]
                        + c3333 * hp43[(k, j)] * v4[i]
                        + h43[(i, k)]
                            * (tg * v3[j] - ctp * bp3[j] - c3 * e23[j])
                        + v3[k]
                            * (tg * h43[(i, j)] + s2g * v4[i] * v3[j]
                                - ctp * hp43[(i, j)]
                                + c2p * bp4[i] * bp3[j]);
                }
            }
        } // end 22 loop

        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    h.h331[(i, j, k)] = c3331 * h31[(i, k)] * bp4[j]
                        + c3333 * hp43[(j, i)] * v1[k]
                        + (tg * v3[i] - ctp * bp3[i] - c3 * e23[i])
                            * h31[(j, k)]
                        + tg * h31[(i, k)] * v3[j]
                        + s2g * v3[i] * v3[j] * v1[k];
                    h.h441[(i, j, k)] = c4411 * h41[(i, k)] * bp3[j]
                        + c4412 * hp43[(i, j)] * v1[k]
                        + (tg * v4[i] - ctp * bp4[i] - c4 * e24[i])
                            * h41[(j, k)]
                        + tg * h41[(i, k)] * v4[j]
                        + s2g * v4[i] * v4[j] * v1[k];
                    h.h443[(i, j, k)] = c4411
                        * (h43[(j, k)] * bp3[i]
                            + hp33[(i, k)] * v4[j]
                            + (c3 * e23[k] - ctp * bp3[k])
                                * (bp3[i] * v4[j] + tg * hp43[(j, i)]))
                        + c4431 * hp334[(i, k, j)]
                        + c4412 * hp43[(j, i)] * v3[k]
                        + h43[(i, k)]
                            * (tg * v4[j] - ctp * bp4[j] - c4 * e24[j])
                        + v4[i]
                            * (tg * h43[(j, k)] + s2g * v4[j] * v3[k]
                                - ctp * hp43[(j, k)]
                                + c2p * bp4[j] * bp3[k]);
                }
            }
        } // end 32 loop

        let prod = Self::tripro();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    h.h431[(i, j, k)] = c4311
                        * (prod[(k, j, i)]
                            + e21[k] * cp21[(i, j)]
                            + e24[i] * cp24[(j, k)]
                            - e24[i] * e21[k] * cp2124[j])
                        + c4312 * v1[k] * (e24[i] * cp2124[j] - cp21[(i, j)])
                        + tg * h41[(i, k)] * v3[j]
                        + s2g * v4[i] * v3[j] * v1[k]
                        + h31[(j, k)] * (tg * v4[i] - ctp * bp4[i])
                        + c4313 * e24[i] * bp4[j] * v1[k]
                        + c3331 * h41[(i, k)] * bp4[j]
                        + c3333 * hp44[(i, j)] * v1[k];
                }
            }
        } // end 42

        h.h112 = -(&h.h111 + &h.h113 + h.h411.view().permuted_axes((1, 2, 0)));
        h.h123 = -(h.h113.clone()
            + h.h331.view().permuted_axes((2, 0, 1))
            + h.h431.view().permuted_axes((2, 0, 1)));
        h.h332 = -(&h.h331 + &h.h333 + h.h433.view().permuted_axes((1, 2, 0)));
        h.h421 = -(&h.h411 + &h.h431 + &h.h441);
        h.h432 = -(&h.h431 + &h.h433 + h.h443.view().permuted_axes((0, 2, 1)));
        h.h442 = -(&h.h441 + &h.h443 + &h.h444);

        h.h221 = -(h.h421.clone()
            + h.h112.view().permuted_axes((0, 2, 1))
            + h.h123.view().permuted_axes((2, 1, 0)));
        h.h223 = -(h.h123.clone()
            + h.h332.view().permuted_axes((0, 2, 1))
            + h.h432.view().permuted_axes((0, 2, 1)));
        h.h422 = -(h.h421.clone()
            + h.h432.view().permuted_axes((0, 2, 1))
            + h.h442.view().permuted_axes((0, 2, 1)));

        h.h222 = -(&h.h221 + &h.h223 + h.h422.view().permuted_axes((1, 2, 0)));
    }

    fn linx(
        geom: &Geom,
        siic: &Siic,
        a: &usize,
        b: &usize,
        c: &usize,
        d: &usize,
        h: &mut Htens,
    ) {
        use Siic::*;
        let qb = geom.unit(*c, *b);
        let qc = geom.unit(*c, *d);
        let r23 = geom.dist(*b, *c);
        let bend = Bend(*a, *b, *c);
        let s = geom.s_vec(&bend);
        splat!(s, q3 => c);
        // hijs2
        let Hmat {
            h31: q31, h32: q32, ..
        } = Hmat::new(geom, &bend);
        // hijks1 call
        let Htens { h111: q444, .. } = Htens::new(geom, &Stretch(*d, *c));
        let q4444 = h4th1(geom, *c, *d);

        let q44 = hijs1(geom, *d, *c);
        let Htens {
            h113: q113,
            h123: q123,
            h223: q223,
            ..
        } = hijks2(geom, *a, *b, *c);

        let Htens4 {
            h1113: q1113,
            h1123: q1123,
            h1223: q1223,
            h2223: q2223,
            ..
        } = h4th2(geom, *a, *b, *c);

        // hijs8 call
        let Hmat {
            h11: q11,
            h21: q21,
            h41: q41,
            h22: q22,
            h42: q42,
            h44: q44a,
            ..
        } = Hmat::new(geom, siic);

        // vect8 call
        let s = geom.s_vec(siic);
        splat!(s, q1 => a, q2 => b, q4 => d);
        let w = siic.value(geom);
        let q22a = hijs1(geom, *b, *c);
        let q222 = hijks1(geom, *b, *c);

        // 4
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        h.h444[(i, j, k)] -= r23 * q4444[(l, i, j, k)] * q3[l];
                        h.h441[(i, j, k)] -=
                            r23 * q444[(l, i, j)] * q31[(l, k)];
                        h.h442[(i, j, k)] -=
                            r23 * q444[(l, i, j)] * q32[(l, k)];
                        h.h111[(i, j, k)] -= r23 * q1113[(i, j, k, l)] * qc[l];
                        h.h112[(i, j, k)] -= r23 * q1123[(i, j, k, l)] * qc[l];
                        h.h221[(i, j, k)] -= r23 * q1223[(k, j, i, l)] * qc[l];
                        h.h222[(i, j, k)] -= r23 * q2223[(i, j, k, l)] * qc[l];
                        h.h421[(i, j, k)] -=
                            r23 * q123[(k, j, l)] * q44[(l, i)];
                        h.h422[(i, j, k)] -=
                            r23 * q223[(j, k, l)] * q44[(l, i)];
                        h.h411[(i, j, k)] -=
                            r23 * q113[(j, k, l)] * q44[(l, i)];
                    }
                    h.h442[(i, j, k)] += q44a[(i, j)] * qb[k] / r23;
                    h.h421[(i, j, k)] += q41[(i, k)] * qb[j] / r23;
                    h.h112[(i, j, k)] += q11[(i, j)] * qb[k] / r23;
                    h.h222[(i, j, k)] += q22[(j, k)] * qb[i] / r23
                        + (q22[(i, k)] * qb[j] + q22[(i, j)] * qb[k]) / r23
                        + (q22a[(i, j)] * q2[k] + q22a[(j, k)] * q2[i]) / r23
                        + q22a[(i, k)] * q2[j] / r23
                        + 6.0 * w * qb[i] * qb[j] * qb[k] / r23.powi(3)
                        - (qb[i] * qb[j] * q2[k]
                            + qb[j] * qb[k] * q2[i]
                            + qb[i] * qb[k] * q2[j])
                            * 2.0
                            / r23.powi(2)
                        + w * q222[(i, j, k)] / r23
                        - 2.0
                            * w
                            * (q22a[(i, j)] * qb[k]
                                + q22a[(i, k)] * qb[j]
                                + q22a[(j, k)] * qb[i])
                            / r23.powi(2);
                    h.h221[(i, j, k)] +=
                        (q21[(i, k)] * qb[j] + q21[(j, k)] * qb[i]) / r23
                            + q22a[(i, j)] * q1[k] / r23
                            - 2.0 * qb[i] * qb[j] * q1[k] / r23.powi(2);
                    h.h422[(i, j, k)] +=
                        (q42[(i, k)] * qb[j] + q42[(i, j)] * qb[k]) / r23
                            + q22a[(j, k)] * q4[i] / r23
                            - 2.0 * qb[j] * qb[k] * q4[i] / r23.powi(2);
                }
            }
        }

        h.h223 = -&h.h222 - &h.h221 - h.h422.view().permuted_axes((1, 2, 0));
        h.h113 = -&h.h112 - &h.h111 - h.h411.view().permuted_axes((1, 2, 0));
        h.h123 = -h.h112.clone().permuted_axes((0, 2, 1))
            - h.h221.view().permuted_axes((2, 1, 0))
            - h.h421.view().permuted_axes((2, 1, 0));
        h.h443 = -&h.h442 - &h.h441 - &h.h444;
        h.h431 = -&h.h421 - &h.h411 - &h.h441;
        h.h432 = -&h.h422 - h.h421.view().permuted_axes((0, 2, 1)) - &h.h442;

        h.h331 = -h.h431.clone()
            - h.h123.view().permuted_axes((1, 2, 0))
            - h.h113.view().permuted_axes((0, 2, 1));
        h.h332 = -h.h432.clone()
            - h.h223.view().permuted_axes((0, 2, 1))
            - h.h123.view().permuted_axes((0, 2, 1));
        h.h433 =
            -h.h431.clone() - &h.h432 - h.h443.view().permuted_axes((1, 2, 0));

        h.h333 = -h.h433.clone()
            - h.h331.view().permuted_axes((2, 0, 1))
            - h.h332.view().permuted_axes((2, 0, 1));
    }

    fn liny(
        geom: &Geom,
        a: &usize,
        b: &usize,
        c: &usize,
        d: &usize,
        h: &mut Htens,
    ) {
        let out = Siic::Out(*d, *c, *b, *a);
        let sout = geom.s_vec(&out);
        splat!(sout, e4 => d, e3 => c, e2 => b);
        let tout = out.value(geom);
        let w = -tout.sin();
        let cosy = tout.cos();

        // hijs7 call
        let Hmat {
            h11: q44,
            h21: q34,
            h31: q24,
            h22: q33,
            h32: q23,
            h33: q22,
            ..
        } = Hmat::new(geom, &out);

        // hijks7 call
        let Htens {
            h111: q444,
            h112: q443,
            h221: q334,
            h222: q333,
            h113: q442,
            h123: q432,
            h223: q332,
            h331: q224,
            h332: q223,
            h333: q222,
            ..
        } = Htens::new(geom, &out);

        // 1
        foreach!(k, i, j,
                h.h222[(i, j, k)] = cosy * e2[i] * e2[j] * e2[k]
                    - cosy * q222[(i, j, k)]
                    - w * (e2[i] * q22[(j, k)]
                        + e2[j] * q22[(i, k)]
                        + e2[k] * q22[(i, j)]);
                h.h223[(i, j, k)] = cosy * e2[i] * e2[j] * e3[k]
                    - cosy * q223[(i, j, k)]
                    - w * (e2[i] * q23[(j, k)]
                        + e2[j] * q23[(i, k)]
                        + e3[k] * q22[(i, j)]);
                h.h422[(i, j, k)] = cosy * e4[i] * e2[j] * e2[k]
                    - cosy * q224[(j, k, i)]
                    - w * (e2[k] * q24[(j, i)]
                        + e2[j] * q24[(k, i)]
                        + e4[i] * q22[(j, k)]);
                h.h333[(i, j, k)] = cosy * e3[i] * e3[j] * e3[k]
                    - cosy * q333[(i, j, k)]
                    - w * (e3[k] * q33[(j, i)]
                        + e3[j] * q33[(k, i)]
                        + e3[i] * q33[(j, k)]);
                h.h433[(i, j, k)] = cosy * e4[i] * e3[j] * e3[k]
                    - cosy * q334[(j, k, i)]
                    - w * (e3[k] * q34[(j, i)]
                        + e3[j] * q34[(k, i)]
                        + e4[i] * q33[(j, k)]);
                h.h332[(i, j, k)] = cosy * e3[i] * e3[j] * e2[k]
                    - cosy * q332[(i, j, k)]
                    - w * (e3[i] * q23[(k, j)]
                        + e3[j] * q23[(k, i)]
                        + e2[k] * q33[(i, j)]);
                h.h432[(i, j, k)] = cosy * e4[i] * e3[j] * e2[k]
                    - cosy * q432[(i, j, k)]
                    - w * (e4[i] * q23[(k, j)]
                        + e3[j] * q24[(k, i)]
                        + e2[k] * q34[(j, i)]);
                h.h444[(i, j, k)] = cosy * e4[i] * e4[j] * e4[k]
                    - cosy * q444[(i, j, k)]
                    - w * (e4[i] * q44[(k, j)]
                        + e4[j] * q44[(k, i)]
                        + e4[k] * q44[(i, j)]);
                h.h443[(i, j, k)] = cosy * e4[i] * e4[j] * e3[k]
                    - cosy * q443[(i, j, k)]
                    - w * (e4[i] * q34[(k, j)]
                        + e4[j] * q34[(k, i)]
                        + e3[k] * q44[(i, j)]);
                h.h442[(i, j, k)] = cosy * e4[i] * e4[j] * e2[k]
                    - cosy * q442[(i, j, k)]
                    - w * (e4[i] * q24[(k, j)]
                        + e4[j] * q24[(k, i)]
                        + e2[k] * q44[(i, j)]);
        );

        h.h221 = -&h.h222 - &h.h223 - h.h422.view().permuted_axes((1, 2, 0));
        h.h331 = -&h.h332 - &h.h333 - h.h433.view().permuted_axes((1, 2, 0));
        h.h123 = -h.h332.clone().permuted_axes((0, 2, 1))
            - &h.h223
            - h.h432.view().permuted_axes((0, 2, 1));
        h.h441 = -&h.h442 - &h.h443 - &h.h444;
        h.h431 = -&h.h432
            - h.h433.view().permuted_axes((0, 2, 1))
            - h.h443.view().permuted_axes((0, 2, 1));
        h.h421 = -&h.h422
            - h.h432.view().permuted_axes((0, 2, 1))
            - h.h442.view().permuted_axes((0, 2, 1));

        h.h112 = -h.h421.clone().permuted_axes((0, 2, 1))
            - h.h123.view().permuted_axes((2, 0, 1))
            - h.h221.view().permuted_axes((0, 2, 1));
        h.h113 = -h.h431.clone().permuted_axes((0, 2, 1))
            - h.h331.view().permuted_axes((0, 2, 1))
            - h.h123.view().permuted_axes((1, 0, 2));
        h.h411 = -&h.h441 - &h.h431 - &h.h421;

        h.h111 = -h.h411.clone().permuted_axes((1, 2, 0)) - &h.h113 - &h.h112;
    }
}
