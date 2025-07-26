use crate::{ANGBOHR, DVec, Siic, Vec3, hmat::Hmat};
use nalgebra as na;
use serde::{Deserialize, Serialize};
use std::{fmt::Display, ops::Index};

#[derive(Debug, Default, PartialEq, Clone, Serialize, Deserialize)]
pub struct Geom(pub Vec<Vec3>);

impl Geom {
    pub fn new() -> Self {
        Geom(Vec::new())
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn push(&mut self, it: Vec3) {
        self.0.push(it)
    }

    /// convert the geometry from angstroms to bohr
    pub fn to_bohr(&mut self) {
        for a in self.0.iter_mut() {
            *a /= ANGBOHR;
        }
    }

    /// return the unit vector from atom i to atom j
    pub fn unit(&self, i: usize, j: usize) -> Vec3 {
        let diff = self[j] - self[i];
        diff / diff.magnitude()
    }

    /// distance between atoms i and j in angstroms, assuming `self` in bohr
    pub fn dist(&self, i: usize, j: usize) -> f64 {
        ANGBOHR * (self[j] - self[i]).magnitude()
    }

    /// return the unit vector from atom i to atom j and the distance between
    /// the atoms
    pub fn vect1(&self, i: usize, j: usize) -> (Vec3, f64) {
        (self.unit(i, j), self.dist(i, j))
    }

    /// angle in radians between atoms i, j, and k, where j is the central atom
    pub fn angle(&self, i: usize, j: usize, k: usize) -> f64 {
        let e_ji = Self::unit(self, j, i);
        let e_jk = Self::unit(self, j, k);
        (e_ji.dot(&e_jk)).acos()
    }

    pub fn s_vec(&self, ic: &Siic) -> Vec<f64> {
        let mut tmp = vec![0.0; 3 * self.len()];
        match ic {
            Siic::Stretch(a, b) => {
                let e_12 = self.unit(*a, *b);
                for i in 0..3 {
                    tmp[3 * a + i] = -e_12[i];
                    tmp[3 * b + i] = e_12[i];
                }
            }
            Siic::Bend(a, b, c) => {
                let e_21 = self.unit(*b, *a);
                let e_23 = self.unit(*b, *c);
                let t_12 = self.dist(*b, *a);
                let t_32 = self.dist(*b, *c);
                let w = e_21.dot(&e_23);
                let sp = (1.0 - w * w).sqrt();
                let c1 = 1.0 / (t_12 * sp);
                let c2 = 1.0 / (t_32 * sp);
                for i in 0..3 {
                    tmp[3 * a + i] = (w * e_21[i] - e_23[i]) * c1;
                    tmp[3 * c + i] = (w * e_23[i] - e_21[i]) * c2;
                    tmp[3 * b + i] = -tmp[3 * a + i] - tmp[3 * c + i];
                }
            }
            Siic::Torsion(a, b, c, d) => {
                let e_21 = self.unit(*b, *a);
                let e_32 = self.unit(*c, *b);
                let e_43 = self.unit(*d, *c);
                let t_21 = self.dist(*b, *a);
                let t_32 = self.dist(*c, *b);
                let t_43 = self.dist(*d, *c);
                let v5 = e_21.cross(&e_32);
                let v6 = e_43.cross(&e_32);
                let cp2 = -e_21.dot(&e_32);
                let cp3 = -e_43.dot(&e_32);
                let sp2 = 1.0 - cp2 * cp2;
                let sp3 = 1.0 - cp3 * cp3;
                // terminal atoms
                let w1 = 1.0 / (t_21 * sp2);
                let w2 = 1.0 / (t_43 * sp3);
                for i in 0..3 {
                    tmp[3 * a + i] = -w1 * v5[i];
                    tmp[3 * d + i] = -w2 * v6[i];
                }
                let w3 = (t_32 - t_21 * cp2) * w1 / t_32;
                let w4 = cp3 / (t_32 * sp3);
                let w5 = (t_32 - t_43 * cp3) * w2 / t_32;
                let w6 = cp2 / (t_32 * sp2);
                for i in 0..3 {
                    tmp[3 * b + i] = w3 * v5[i] + w4 * v6[i];
                    tmp[3 * c + i] = w5 * v6[i] + w6 * v5[i];
                }
            }
            Siic::Lin1(a, b, c, d) => {
                let e21 = self.unit(*b, *a);
                let e23 = self.unit(*c, *b);
                let t21 = self.dist(*b, *a);
                let t23 = self.dist(*c, *b);
                let ea = self[*d];
                let d = 1.0 / ea.dot(&ea).sqrt();
                let ea = d * ea;
                let e2m = e23.cross(&e21);
                let stheta = ea.dot(&e2m);
                let w = f64::asin(stheta);
                let ctheta = w.cos();
                let ttheta = stheta / ctheta;
                let v4 = ea.cross(&e23);
                let v5 = ea.cross(&e21);
                let c1 = 1.0 / (ctheta * t21);
                let c2 = ttheta / t21;
                let c3 = 1.0 / (ctheta * t23);
                let c4 = ttheta / t23;
                for i in 0..3 {
                    // V1, V3, V2 in VECT3
                    let v1i = c1 * v4[i] - c2 * e21[i];
                    let v3i = -(c3 * v5[i] + c4 * e23[i]);
                    tmp[3 * a + i] = -v1i;
                    tmp[3 * c + i] = v3i;
                    tmp[3 * b + i] = -(-v1i + v3i);
                }
            }
            Siic::Out(a, b, c, d) => {
                let e21 = self.unit(*a, *b);
                let e23 = self.unit(*c, *b);
                let e24 = self.unit(*d, *b);
                let t21 = self.dist(*a, *b);
                let t23 = self.dist(*c, *b);
                let t24 = self.dist(*d, *b);
                let v5 = e23.cross(&e24);
                let w1 = e21.dot(&e23);
                let w2 = e21.dot(&e24);
                let phi = self.angle(*c, *b, *d);
                let sphi = phi.sin();
                let w = e21.dot(&v5);
                let w = (w / sphi).asin();
                let w = if w1 + w2 > 0.0 {
                    std::f64::consts::PI.copysign(w) - w
                } else {
                    w
                };
                let cg = w.cos();
                let sg = w.sin();
                let tg = sg / cg;
                let w1 = cg * sphi;
                let w2 = 1.0 / (t21 * w1);
                let w3 = tg / t21;
                let w4 = 1.0 / (t23 * w1);
                let w5 = t24 * sg * w4;
                let w6 = 1.0 / (t24 * w1);
                let w7 = t23 * sg * w6;
                let v6 = e24.cross(&e21);
                let v7 = e21.cross(&e23);
                let svec = self.s_vec(&Siic::Bend(*c, *b, *d));
                let b3p = &svec[3 * c..3 * c + 3];
                let b4p = &svec[3 * d..3 * d + 3];
                for i in 0..3 {
                    // V1, V3, V4, V2 in VECT5
                    let v1i = v5[i] * w2 - e21[i] * w3;
                    let v3i = v6[i] * w4 - b4p[i] * w5;
                    let v4i = v7[i] * w6 - b3p[i] * w7;
                    tmp[3 * a + i] = v1i;
                    tmp[3 * c + i] = v3i;
                    tmp[3 * d + i] = v4i;
                    tmp[3 * b + i] = -v1i - v3i - v4i;
                }
            }
            &Siic::Linx(a, b, c, d) => {
                let e32 = self.unit(c, b);
                let e34 = self.unit(c, d);
                let t32 = self.dist(c, b);
                let s = self.s_vec(&Siic::Bend(a, b, c));
                let s3 = &s[3 * c..3 * c + 3];
                let e3 = na::vector![s3[0], s3[1], s3[2]];
                let w = -t32 * e34.dot(&e3);
                let h1 = Hmat::new(self, &Siic::Stretch(c, d));
                let h2 = Hmat::new(self, &Siic::Bend(a, b, c));
                for i in 0..3 {
                    tmp[3 * a + i] = 0.0;
                    tmp[3 * d + i] = 0.0;
                    tmp[3 * b + i] = w * e32[i] / t32;
                    for j in 0..3 {
                        tmp[3 * a + i] -= t32 * e34[j] * h2.h31[(j, i)];
                        tmp[3 * b + i] -= t32 * e34[j] * h2.h32[(j, i)];
                        tmp[3 * d + i] -= t32 * e3[j] * h1.h11[(j, i)];
                    }
                    tmp[3 * c + i] =
                        -tmp[3 * a + i] - tmp[3 * b + i] - tmp[3 * d + i];
                }
            }
            &Siic::Liny(a, b, c, d) => {
                let out = &Siic::Out(d, c, b, a);
                let tout = out.value(self);
                let cosy = tout.cos();
                let s = self.s_vec(out);
                let e1 = &s[3 * a..3 * a + 3];
                let e2 = &s[3 * b..3 * b + 3];
                let e3 = &s[3 * c..3 * c + 3];
                let e4 = &s[3 * d..3 * d + 3];
                for i in 0..3 {
                    tmp[3 * a + i] = -cosy * e1[i];
                    tmp[3 * b + i] = -cosy * e2[i];
                    tmp[3 * c + i] = -cosy * e3[i];
                    tmp[3 * d + i] = -cosy * e4[i];
                }
            }
        }
        tmp
    }
}

impl From<psqs::geom::Geom> for Geom {
    /// panics if any element of `geom` has a length other than 3 or if `geom`
    /// is not Cartesian
    fn from(geom: psqs::geom::Geom) -> Self {
        let geom = geom.xyz().unwrap();
        Self(
            geom.iter()
                .map(|atom| Vec3::from_row_slice(&atom.coord()))
                .collect(),
        )
    }
}

impl From<symm::Molecule> for Geom {
    fn from(geom: symm::Molecule) -> Self {
        Self(
            geom.atoms
                .iter()
                .map(|atom| Vec3::from_row_slice(&atom.coord()))
                .collect(),
        )
    }
}

impl From<&DVec> for Geom {
    fn from(dvec: &DVec) -> Self {
        Self(
            dvec.as_slice()
                .chunks(3)
                .map(Vec3::from_row_slice)
                .collect(),
        )
    }
}

impl From<Geom> for DVec {
    fn from(val: Geom) -> Self {
        let mut geom = Vec::with_capacity(val.len());
        for c in &val {
            geom.extend(&c);
        }
        DVec::from(geom)
    }
}

impl IntoIterator for &Geom {
    type Item = Vec3;

    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.clone().into_iter()
    }
}

impl Index<usize> for Geom {
    type Output = Vec3;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Display for Geom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for atom in &self.0 {
            writeln!(f, "{:20.10}{:20.10}{:20.10}", atom[0], atom[1], atom[2])?;
        }
        Ok(())
    }
}
