use std::fmt::Display;

use crate::{DMat, Siic, Vec3, dsplat, geom::Geom, htens::Htens};

#[derive(Debug)]
pub struct Hmat {
    pub h11: DMat,
    pub h21: DMat,
    pub h31: DMat,
    pub h22: DMat,
    pub h32: DMat,
    pub h33: DMat,
    pub h41: DMat,
    pub h42: DMat,
    pub h43: DMat,
    pub h44: DMat,
}

/// helper function for calling Hmat::new on `geom` with an [Siic::Stretch] made
/// from atoms `i` and `j`
pub fn hijs1(geom: &Geom, i: usize, j: usize) -> DMat {
    Hmat::new(geom, &Siic::Stretch(i, j)).h11
}

pub fn hijs2(geom: &Geom, i: usize, j: usize, k: usize) -> Hmat {
    Hmat::new(geom, &Siic::Bend(i, j, k))
}

impl Hmat {
    pub fn zeros() -> Self {
        Self {
            h11: DMat::zeros(3, 3),
            h21: DMat::zeros(3, 3),
            h31: DMat::zeros(3, 3),
            h22: DMat::zeros(3, 3),
            h32: DMat::zeros(3, 3),
            h33: DMat::zeros(3, 3),
            h41: DMat::zeros(3, 3),
            h42: DMat::zeros(3, 3),
            h43: DMat::zeros(3, 3),
            h44: DMat::zeros(3, 3),
        }
    }

    // making block matrices to pack into sr in machx
    pub fn new(geom: &Geom, s: &Siic) -> Self {
        use crate::Siic::*;
        let mut h = Self::zeros();
        match s {
            // from HIJS1
            Stretch(a, b) => {
                let v1 = geom.unit(*a, *b);
                let t21 = geom.dist(*a, *b);
                h.h11 = DMat::identity(3, 3);
                h.h11 -= v1 * v1.transpose();
                h.h11 /= t21;
            }
            // from HIJS2
            Bend(i, j, k) => {
                let s = geom.s_vec(s);
                dsplat!(s, v1 => i, v3 => k);
                let e21 = geom.unit(*j, *i);
                let e23 = geom.unit(*j, *k);
                let t21 = geom.dist(*j, *i);
                let t23 = geom.dist(*j, *k);
                let h11a = Self::new(geom, &Stretch(*i, *j)).h11;
                let h33a = Self::new(geom, &Stretch(*k, *j)).h11;
                let phi = geom.angle(*i, *j, *k);
                let sphi = phi.sin();
                let w1 = phi.cos() / sphi;
                let w2 = 1.0 / t21;
                let w3 = w1 * w2;
                let w4 = 1.0 / t23;
                let w5 = w1 * w4;
                let w6 = 1.0 / (t21 * sphi);
                h.h11 = w3 * &h11a;
                h.h11 += -w1 * &v1 * v1.transpose()
                    - w2 * (e21 * v1.transpose() + &v1 * e21.transpose());
                //
                h.h33 = w5 * &h33a;
                h.h33 += -w1 * &v3 * v3.transpose()
                    - w4 * (e23 * v3.transpose() + &v3 * e23.transpose());
                //
                h.h31 = -w6 * h33a;
                h.h31 -= v3 * (w1 * v1 + w2 * e21).transpose();
                //
                h.h21 = -(&h.h11 + &h.h31);
                h.h32 = -(&h.h31 + &h.h33);
                h.h22 = -(h.h21.transpose() + &h.h32);
            }
            // from HIJS6
            Torsion(i, j, k, l) => {
                Self::torsion(geom, s, i, j, k, l, &mut h);
            }
            // from HIJS3
            Lin1(i, j, k, l) => {
                let tmp = geom.s_vec(s);
                dsplat!(tmp, v1 => i, v3 => k);
                let th = s.value(geom);
                let e21 = geom.unit(*j, *i);
                let e23 = geom.unit(*j, *k);
                let t21 = geom.dist(*j, *i);
                let t23 = geom.dist(*j, *k);
                let h11a = Self::new(geom, &Stretch(*i, *j)).h11;
                let h33a = Self::new(geom, &Stretch(*k, *j)).h11;
                let ea = geom[*l] / geom[*l].magnitude();
                let tanth = th.tan();
                let costh = th.cos();
                let em = Self::mat1(&-ea);
                let w1 = 1.0 / t21;
                let w2 = 1.0 / t23;
                h.h11 = tanth * (-h11a * w1 + &v1 * v1.transpose());
                h.h11 += -w1 * (e21 * v1.transpose() + &v1 * e21.transpose());

                h.h31 = (em * &h33a).transpose() / costh;
                h.h31 += -&v3 * e21.transpose();
                h.h31 /= t21;
                h.h31 += &v3 * v1.transpose() * tanth;

                h.h33 = tanth * (-h33a * w2 + &v3 * v3.transpose());
                h.h33 -= w2 * (e23 * v3.transpose() + v3 * e23.transpose());

                h.h21 = -(&h.h11 + &h.h31);
                h.h32 = -(&h.h31 + &h.h33);
                h.h22 = -(h.h21.transpose() + &h.h32);
            }
            // from HIJS7
            Out(i, j, k, l) => Self::out(geom, i, j, k, l, s, &mut h),
            // hijs8
            &Linx(i, j, k, l) => {
                Self::linx(geom, i, j, k, l, &mut h);
            }
            // hijs9
            &Liny(k1, k2, k3, k4) => {
                // vect5 call
                let out = Siic::Out(k4, k3, k2, k1);
                let tout = out.value(geom);
                let s = geom.s_vec(&out);
                dsplat!(s, e1 => k1, e2 => k2, e3 => k3, e4 => k4);
                let w = tout.sin();
                let cosy = tout.cos();
                let Hmat {
                    h11: q44,
                    h21: q34,
                    h31: q24,
                    h41: q14,
                    h22: q33,
                    h32: q23,
                    h42: q13,
                    h33: q22,
                    h43: q12,
                    h44: q11,
                } = Hmat::new(geom, &out);
                h.h11 = w * &e1 * e1.transpose() - cosy * q11.transpose();
                h.h33 = w * &e3 * e3.transpose() - cosy * q33.transpose();
                h.h44 = w * &e4 * e4.transpose() - cosy * q44.transpose();
                h.h31 = w * &e3 * e1.transpose() - cosy * q13.transpose();
                h.h41 = w * &e4 * e1.transpose() - cosy * q14.transpose();
                h.h43 = w * &e4 * e3.transpose() - cosy * q34.transpose();
                h.h21 = w * &e2 * e1.transpose() - cosy * q12.transpose();
                h.h32 = w * &e3 * e2.transpose() - cosy * q23.transpose();
                h.h42 = w * &e4 * e2.transpose() - cosy * q24.transpose();
                h.h22 = w * &e2 * e2.transpose() - cosy * q22.transpose();
            }
        }
        h
    }

    /// helper function for constructing some kind of matrix
    pub fn mat1(v: &Vec3) -> DMat {
        let mut em = DMat::zeros(3, 3);
        em[(1, 0)] = -v[2];
        em[(2, 0)] = v[1];
        em[(2, 1)] = -v[0];
        em[(0, 1)] = -em[(1, 0)];
        em[(0, 2)] = -em[(2, 0)];
        em[(1, 2)] = -em[(2, 1)];
        em
    }

    fn out(
        geom: &Geom,
        a: &usize,
        b: &usize,
        c: &usize,
        d: &usize,
        s: &Siic,
        h: &mut Self,
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
        dsplat!(svec, bp3 => c, bp4 => d);

        // vect5 call
        let svec = geom.s_vec(s);
        dsplat!(svec, v1 => a, v3 => c, v4 => d);
        let gamma = s.value(geom);

        // hijs2 call
        let Hmat {
            h31: hp43,
            h33: hp44,
            ..
        } = Hmat::new(geom, &Siic::Bend(*c, *b, *d));

        let v5 = e23.cross(&e24);
        let v6 = e24.cross(&e21);
        let cp31 = Self::mat1(&e24);
        let cp41 = Self::mat1(&e23);
        let cp43 = Self::mat1(&e21);
        let sp = phi.sin();
        let cg = gamma.cos();
        let tg = gamma.tan();
        let c21 = 1.0 / t21;
        let c23 = 1.0 / t23;
        let c24 = 1.0 / t24;
        let ctp = 1.0 / (sp / phi.cos());
        let c11 = tg * c21 * c21;
        let c312 = c21 / (cg * sp);
        let c311 = c312 * c23;
        let c313 = c312 * ctp;
        let c411 = c312 * c24;
        let c3 = c23 / sp;
        let c331 = t24 * c3;
        let c332 = c331 * tg;
        let c441 = t23 * (c24 / sp);
        let c442 = c441 * tg;
        let c431 = c3 * c24 / cg;
        let c432 = tg;
        let c434 = tg * c3;
        let c435 = t24 * c3;
        let c436 = c435 * tg;

        h.h11 = -c11 * DMat::identity(3, 3);
        h.h11 += (tg * &v1 - c21 * e21) * v1.transpose()
            + c11 * e21 * e21.transpose()
            - c21 * &v1 * e21.transpose();

        h.h33 = c331 * &v3 * bp4.transpose() + c332 * hp43.transpose();
        h.h33 += (tg * &v3 - e23 * c23 - &bp3 * ctp) * v3.transpose();

        h.h44 = &v4 * bp3.transpose() * c441 + hp43 * c442;
        h.h44 += (tg * &v4 - e24 * c24 - &bp4 * ctp) * v4.transpose();

        let xj = (tg * &v1 - e21 * c21).transpose();
        h.h31 = -cp31 * c311;
        h.h31 += &v3 * xj - (c311 * e23 + c313 * bp3) * v5.transpose();

        h.h41 = cp41 * c411;
        h.h41 += &v4 * xj - (c411 * e24 + c313 * &bp4) * v5.transpose();

        h.h43 = c431 * cp43.transpose();
        h.h43 += -c431 * e24 * v6.transpose()
            + (c432 * &v4 - ctp * &bp4) * v3.transpose()
            + (c434 * e24 + c435 * v4) * bp4.transpose()
            + c436 * hp44;

        h.h21 = -(&h.h11 + &h.h31 + &h.h41);
        h.h32 = -(&h.h31 + &h.h33 + h.h43.transpose());
        h.h42 = -(&h.h41 + &h.h43 + h.h44.transpose());
        h.h22 = -(h.h21.transpose() + &h.h32 + &h.h42);
    }

    fn torsion(
        geom: &Geom,
        s: &Siic,
        a: &usize,
        b: &usize,
        c: &usize,
        d: &usize,
        h: &mut Hmat,
    ) {
        // unpack the s vector. mine are in the opposite order of the
        // fortran
        let tmp = geom.s_vec(s);
        dsplat!(tmp, v1 => a, v4 => d);
        // unit and non-unit vectors
        let e21 = geom.unit(*b, *a);
        let e23 = geom.unit(*b, *c);
        let e34 = geom.unit(*c, *d);
        let t21 = geom.dist(*b, *a);
        let t23 = geom.dist(*b, *c);
        let t34 = geom.dist(*c, *d);
        // angles
        let p2 = geom.angle(*a, *b, *c);
        let tmp = geom.s_vec(&Siic::Bend(*a, *b, *c));
        dsplat!(tmp, bp21 => a, bp22 => b, bp23 => c);

        let p3 = geom.angle(*b, *c, *d);
        let tmp = geom.s_vec(&Siic::Bend(*b, *c, *d));
        dsplat!(tmp, bp32 => b, bp34 => d);

        let x = t21 * p2.sin().powi(2);
        let y = t34 * p3.sin().powi(2);
        let w1 = 1.0 / (t21 * x);
        let w2 = 1.0 / (t23 * x);
        let w3 = 1.0 / (t34 * y);
        let w4 = 1.0 / (t23 * y);

        h.h11 = -w1 * Self::mat1(&e23);
        h.h31 = w2 * Self::mat1(&e21);
        h.h44 = w3 * Self::mat1(&e23);
        h.h42 = -w4 * Self::mat1(&e34);

        // these are cotans
        let xx = p2.cos() / p2.sin();
        let xy = p3.cos() / p3.sin();
        h.h11 -= 2.0 * (e21 / t21 + bp21 * xx) * v1.transpose();
        h.h31 -= (e23 / t23 + 2.0 * bp23 * xx) * v1.transpose();
        h.h44 -= 2.0 * (e34 / t34 + &bp34 * xy) * v4.transpose();
        h.h42 += &v4 * (e23 / t23 - 2.0 * &bp32 * xy).transpose();

        h.h21 = -(&h.h11 + &h.h31);
        h.h43 = -(&h.h44 + &h.h42);

        let x1 = t21 / t23;
        let y1 = t34 / t23;
        let x2 = p2.cos();
        let y2 = p2.sin();
        let x3 = p3.cos();
        let y3 = p3.sin();
        let c1 = x1 * x2 - 1.0;
        let c2 = -x3 * y1;
        let c3 = -x2 / t23;
        let c4 = -x1 * y2;
        let c5 = x1 * x2 / t23;
        let c6 = y1 * y3;
        let c7 = -y1 * x3 / t23;

        h.h22 = c1 * &h.h21 + c2 * h.h42.transpose();
        h.h22 += (c3 * e21 + c4 * bp22 + c5 * e23) * v1.transpose()
            + (c6 * bp32 + c7 * e23) * v4.transpose();

        h.h32 = -(h.h21.transpose() + &h.h22 + &h.h42);
        h.h33 = -(&h.h31 + &h.h32 + h.h43.transpose());
    }

    /// compute the h matrix for a LINX
    fn linx(geom: &Geom, a: usize, b: usize, c: usize, d: usize, h: &mut Hmat) {
        let e2 = geom.unit(c, b);
        let e4 = geom.unit(c, d);
        let t32 = geom.dist(b, c);
        // vect2 call
        let bend = Siic::Bend(a, b, c);
        dsplat!(geom.s_vec(&bend), q3 => c);
        let t = e4.dot(&q3);
        let w = -t32 * t;
        let stre = Siic::Stretch(d, c);
        let Hmat { h11: e44, .. } = Self::new(geom, &stre);
        let Hmat {
            h31: q31, h32: q32, ..
        } = Self::new(geom, &bend);
        let Hmat { h11: e22, .. } = Self::new(geom, &Siic::Stretch(b, c));
        let Htens { h111: q444, .. } = Htens::new(geom, &stre);
        for j in 0..3 {
            for k in 0..3 {
                for i in 0..3 {
                    h.h44[(j, k)] -= t32 * q444[(i, j, k)] * q3[i];
                }
            }
        }
        let Htens {
            h113: q113,
            h123: q123,
            h223: q223,
            ..
        } = Htens::new(geom, &bend);

        h.h22 = w * e22 / t32;
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    h.h11[(j, k)] -= t32 * e4[i] * q113[(j, k, i)];
                    h.h21[(j, k)] -=
                        e4[i] * (e2[j] * q31[(i, k)] + t32 * q123[(k, j, i)]);
                    h.h22[(j, k)] -= e4[i]
                        * (e2[j] * q32[(i, k)]
                            + e2[k] * q32[(i, j)]
                            + t32 * q223[(j, k, i)]);
                }
            }
        }

        h.h41 -= t32 * e44.transpose() * q31;
        h.h42 -= e44.transpose() * (t32 * q32 + q3 * e2.transpose());
        h.h31 = -&h.h11 - &h.h21 - &h.h41;
        h.h32 = -h.h21.transpose() - &h.h22 - &h.h42;
        h.h43 = -&h.h41 - &h.h42 - &h.h44;
        h.h33 = -&h.h31 - &h.h32 - h.h43.transpose();
    }
}

impl Display for Hmat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "H11 =\n{:.6}", self.h11)?;
        writeln!(f, "H21 =\n{:.6}", self.h21)?;
        writeln!(f, "H31 =\n{:.6}", self.h31)?;
        writeln!(f, "H22 =\n{:.6}", self.h22)?;
        writeln!(f, "H32 =\n{:.6}", self.h32)?;
        writeln!(f, "H33 =\n{:.6}", self.h33)?;
        writeln!(f, "H41 =\n{:.6}", self.h41)?;
        writeln!(f, "H42 =\n{:.6}", self.h42)?;
        writeln!(f, "H43 =\n{:.6}", self.h43)?;
        writeln!(f, "H44 =\n{:.6}", self.h44)?;
        Ok(())
    }
}
