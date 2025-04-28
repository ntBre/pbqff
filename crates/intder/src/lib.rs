use htens4::fill4a;
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    error::Error,
    fmt::{Display, Formatter},
    fs::File,
    io::{BufRead, BufReader, Read, Write},
    sync::LazyLock,
};

pub mod geom;
pub mod hmat;
pub mod htens;
pub mod htens4;

use geom::Geom;
use hmat::Hmat;
use htens::Htens;
use nalgebra as na;
use regex::Regex;
use symm::{Axis, Irrep};

use crate::htens::utils::fill3a;
pub type Tensor3 = ndarray::Array3<f64>;
pub type Tensor4 = ndarray::Array4<f64>;

/// from <https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0>
pub const ANGBOHR: f64 = 0.529_177_210_9;
/// constants from the fortran version
pub const HART: f64 = 4.3597482;
// const DEBYE: f64 = 2.54176548;

// flags
pub static VERBOSE: LazyLock<bool> =
    LazyLock::new(|| std::env::var("INTDER_DEBUG").is_ok());

pub fn is_verbose() -> bool {
    *VERBOSE
}

// TODO make these input or flag params
const TOLDISP: f64 = 1e-14;
const MAX_ITER: usize = 20;

type Vec3 = na::Vector3<f64>;
pub type DMat = na::DMatrix<f64>;
pub type DVec = na::DVector<f64>;

#[derive(Debug, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub enum Siic {
    /// bond stretch between two atoms
    Stretch(usize, usize),

    /// central atom is second like normal people would expect
    Bend(usize, usize, usize),

    /// angle between planes formed by i, j, k and j, k, l
    Torsion(usize, usize, usize, usize),

    /// linear bend of atoms `i`, `j`, and `k`, about `d`, a dummy atom
    /// perpendicular to the line `i`-`j`-`k` and extending from atom `j`
    Lin1(usize, usize, usize, usize),

    /// bend of atom `i` out of the plane formed by atoms `j`, `k`, and `l`
    Out(usize, usize, usize, usize),

    /// the x component of the c → d unit vector in the local coordinate system
    /// in which the b → c vector defines the +z axis and the a atom lies in the
    /// xz plane in the +x direction
    Linx(usize, usize, usize, usize),

    /// the y component of the c → d unit vector in the local coordinate system
    /// in which the b → c vector defines the +z axis and the a atom lies in the
    /// xz plane in the +x direction
    Liny(usize, usize, usize, usize),
}

impl Display for Siic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut lin_helper = |i, j, k, l, x| {
            write!(f, "LIN{}({}-{}-{}-{})", x, i + 1, j + 1, k + 1, l + 1)
        };
        match &self {
            Siic::Stretch(i, j) => write!(f, "r({}-{})", i + 1, j + 1),
            Siic::Bend(i, j, k) => {
                write!(f, "∠({}-{}-{})", i + 1, j + 1, k + 1)
            }
            Siic::Torsion(i, j, k, l) => {
                write!(f, "τ({}-{}-{}-{})", i + 1, j + 1, k + 1, l + 1)
            }
            Siic::Lin1(i, j, k, l) => lin_helper(i, j, k, l, "1"),
            Siic::Out(i, j, k, l) => {
                write!(f, "OUT({}-{}-{}-{})", i + 1, j + 1, k + 1, l + 1)
            }
            Siic::Linx(i, j, k, l) => lin_helper(i, j, k, l, "X"),
            Siic::Liny(i, j, k, l) => lin_helper(i, j, k, l, "Y"),
        }
    }
}

impl Siic {
    pub fn value(&self, geom: &Geom) -> f64 {
        use Siic::*;
        match self {
            Stretch(a, b) => geom.dist(*a, *b),
            Bend(a, b, c) => geom.angle(*a, *b, *c),
            // vect6
            Torsion(a, b, c, d) => {
                let e_21 = geom.unit(*b, *a);
                let e_32 = geom.unit(*c, *b);
                let e_43 = geom.unit(*d, *c);
                let v5 = e_21.cross(&e_32);
                let v6 = e_43.cross(&e_32);
                let w2 = e_21.dot(&e_32);
                let w3 = e_43.dot(&e_32);
                let sp2 = (1.0 - w2 * w2).sqrt();
                let sp3 = (1.0 - w3 * w3).sqrt();
                let w2 = e_21.dot(&v6);
                let w3 = -v5.dot(&v6);
                let w = w2 / (sp2 * sp3);
                let w_size = w.abs() - 1.0;
                let w = w.clamp(-1.0, 1.0).asin();
                if w.is_nan() {
                    panic!("nan calling sin on {w_size}");
                }
                if w3 < 0.0 {
                    std::f64::consts::PI - w
                } else {
                    w
                }
            }
            // vect3
            Lin1(a, b, c, d) => {
                let e21 = geom.unit(*b, *a);
                let e23 = geom.unit(*b, *c);
                let ea = geom[*d] / geom[*d].magnitude();
                let e2m = e23.cross(&e21);
                let stheta = ea.dot(&e2m);
                f64::asin(stheta)
            }
            // vect4
            Out(a, b, c, d) => {
                let e21 = geom.unit(*b, *a);
                let e23 = geom.unit(*b, *c);
                let e24 = geom.unit(*b, *d);
                let v5 = e23.cross(&e24);
                let w1 = e21.dot(&e23);
                let w2 = e21.dot(&e24);
                let phi = geom.angle(*c, *b, *d);
                let sphi = phi.sin();
                let w = e21.dot(&v5);
                let w = (w / sphi).asin();
                if w1 + w2 > 0.0 {
                    std::f64::consts::PI.copysign(w) - w
                } else {
                    w
                }
            }
            // vect8
            Linx(a, b, c, d) => {
                let e34 = geom.unit(*c, *d);
                let t32 = geom.dist(*c, *b);
                let s = geom.s_vec(&Self::Bend(*a, *b, *c));
                let s3 = &s[3 * c..3 * c + 3];
                let e3 = na::vector![s3[0], s3[1], s3[2]];
                -t32 * e34.dot(&e3)
            }
            // vect9
            Liny(a, b, c, d) => -Out(*d, *c, *b, *a).value(geom).sin(),
        }
    }
}

static DEFAULT_WEIGHTS: LazyLock<HashMap<&'static str, usize>> =
    LazyLock::new(|| {
        HashMap::from([
            ("X", 0),
            ("H", 1),
            ("He", 4),
            ("Li", 7),
            ("Be", 9),
            ("B", 11),
            ("C", 12),
            ("N", 14),
            ("O", 16),
            ("F", 19),
            ("Ne", 20),
            ("Na", 23),
            ("Mg", 24),
            ("Al", 27),
            ("Si", 28),
            ("P", 31),
            ("S", 32),
            ("Cl", 35),
            ("Ar", 40),
        ])
    });

#[derive(Debug, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub label: String,
    pub weight: usize,
}

impl Atom {
    pub fn new(label: impl Into<String>, weight: usize) -> Self {
        Self {
            label: label.into(),
            weight,
        }
    }
}

impl From<&symm::Atom> for Atom {
    fn from(value: &symm::Atom) -> Self {
        let label = value.label();
        Self {
            label: label.to_string(),
            weight: DEFAULT_WEIGHTS[label],
        }
    }
}

impl Display for Intder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use Siic::*;
        writeln!(f, "# INTDER ###############")?;
        for (i, op) in self.input_options.iter().enumerate() {
            if i == 16 {
                writeln!(f)?;
            }
            write!(f, "{op:5}")?;
        }
        writeln!(f)?;
        let lin_helper = |i, j, k, l, x, f: &mut Formatter<'_>| {
            writeln!(f, " LIN{x}{:5}{:5}{:5}{:5}", i + 1, j + 1, k + 1, l + 1)
        };
        for siic in &self.simple_internals {
            match siic {
                Stretch(i, j) => {
                    writeln!(f, "{:<5}{:5}{:5}", "STRE", i + 1, j + 1)?
                }
                Bend(i, j, k) => writeln!(
                    f,
                    "{:<5}{:5}{:5}{:5}",
                    "BEND",
                    i + 1,
                    j + 1,
                    k + 1
                )?,
                Torsion(i, j, k, l) => writeln!(
                    f,
                    "{:<5}{:5}{:5}{:5}{:5}",
                    "TORS",
                    i + 1,
                    j + 1,
                    k + 1,
                    l + 1
                )?,
                Lin1(i, j, k, l) => lin_helper(i, j, k, l, "1", f)?,
                Out(i, j, k, l) => writeln!(
                    f,
                    "{:<5}{:5}{:5}{:5}{:5}",
                    "OUT",
                    i + 1,
                    j + 1,
                    k + 1,
                    l + 1
                )?,
                Linx(i, j, k, l) => lin_helper(i, j, k, l, "X", f)?,
                Liny(i, j, k, l) => lin_helper(i, j, k, l, "Y", f)?,
            }
        }
        for (i, sic) in self.symmetry_internals.iter().enumerate() {
            write!(f, "{:5}", i + 1)?;
            for (j, s) in sic.iter().enumerate() {
                if *s != 0.0 {
                    let sign = s.signum();
                    write!(f, "{:4}{:14.9}", j + 1, sign)?;
                }
            }
            writeln!(f)?;
        }
        writeln!(f, "{:5}", 0)?;
        write!(f, "{}", self.geom)?;
        if !self.disps.is_empty() {
            writeln!(f, "DISP{:4}", self.disps.len())?;
            for disp in &self.disps {
                for (i, d) in disp.iter().enumerate() {
                    if *d != 0.0 {
                        writeln!(f, "{:5}{:20.10}", i + 1, d)?;
                    }
                }
                writeln!(f, "{:5}", 0)?;
            }
        } else {
            // assume freqs
            for (i, Atom { label, weight }) in self.atoms.iter().enumerate() {
                let width = match i {
                    0 => 11,
                    1 => 13,
                    _ => 12,
                };
                write!(
                    f,
                    "{:>width$}",
                    format!("{}{}", label.to_uppercase(), weight)
                )?;
            }
            writeln!(f)?;
            let nsic = self.symmetry_internals.len();
            for i in 1..=nsic {
                for j in 1..=i {
                    if let Some(v) = self.fc2.get(fc2_index(nsic, i, j)) {
                        writeln!(f, "{:5}{:5}{:5}{:5}{:20.12}", i, j, 0, 0, v)?;
                    }
                }
            }
            writeln!(f, "{:5}", 0)?;
            for i in 1..=nsic {
                for j in 1..=i {
                    for k in 1..=j {
                        if let Some(v) = self.fc3.get(fc3_index(i, j, k)) {
                            writeln!(
                                f,
                                "{:5}{:5}{:5}{:5}{:20.12}",
                                i, j, k, 0, v,
                            )?;
                        }
                    }
                }
            }
            writeln!(f, "{:5}", 0)?;
            for i in 1..=nsic {
                for j in 1..=i {
                    for k in 1..=j {
                        for l in 1..=k {
                            if let Some(v) = self.fc4.get(fc4_index(i, j, k, l))
                            {
                                writeln!(f, "{i:5}{j:5}{k:5}{l:5}{v:20.12}")?;
                            }
                        }
                    }
                }
            }
            writeln!(f, "{:5}", 0)?;
        }
        Ok(())
    }
}

#[derive(Debug, Default, PartialEq, Clone, Serialize, Deserialize)]
pub struct Intder {
    pub input_options: Vec<usize>,
    pub simple_internals: Vec<Siic>,
    pub symmetry_internals: Vec<Vec<f64>>,
    /// cartesian geometry in bohr
    pub geom: Geom,
    /// SIC displacements to be converted to Cartesian coordinates
    pub disps: Vec<Vec<f64>>,
    /// Atom labels and weights, for use in force constant conversions
    pub atoms: Vec<Atom>,
    /// second order SIC force constants to be converted to cartesian
    /// coordinates
    pub fc2: Vec<f64>,
    /// third order SIC force constants to be converted to cartesian
    /// coordinates
    pub fc3: Vec<f64>,
    /// fourth order SIC force constants to be converted to cartesian
    /// coordinates
    pub fc4: Vec<f64>,
}

/// compute the index in the second-order force constant array, assuming `i`
/// and `j` have minimum values of 1
pub fn fc2_index(ncoords: usize, i: usize, j: usize) -> usize {
    let mut sp = [i, j];
    sp.sort();
    ncoords * (sp[0] - 1) + sp[1] - 1
}

pub fn fc3_index(i: usize, j: usize, k: usize) -> usize {
    let mut sp = [i, j, k];
    sp.sort();
    sp[0] + (sp[1] - 1) * sp[1] / 2 + (sp[2] - 1) * sp[2] * (sp[2] + 1) / 6 - 1
}

pub fn fc4_index(i: usize, j: usize, k: usize, l: usize) -> usize {
    let mut sp = [i, j, k, l];
    sp.sort();
    sp[0]
        + (sp[1] - 1) * sp[1] / 2
        + (sp[2] - 1) * sp[2] * (sp[2] + 1) / 6
        + (sp[3] - 1) * sp[3] * (sp[3] + 1) * (sp[3] + 2) / 24
        - 1
}

fn get_fc3(v: &[f64], i: usize, j: usize, k: usize) -> f64 {
    *v.get(fc3_index(i + 1, j + 1, k + 1)).unwrap_or(&0.0)
}

fn get_fc4(v: &[f64], i: usize, j: usize, k: usize, l: usize) -> f64 {
    *v.get(fc4_index(i + 1, j + 1, k + 1, l + 1)).unwrap_or(&0.0)
}

/// the opposite of the Kronecker delta, 0 when equal, 1 otherwise
#[inline]
fn delta(i: usize, j: usize) -> f64 {
    if i == j {
        0.0
    } else {
        1.0
    }
}

#[derive(Debug)]
pub enum IntderError {
    DispError(String),
    FreqError,
}

impl Display for IntderError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl Error for IntderError {}

impl Intder {
    pub fn new() -> Self {
        Intder {
            input_options: Vec::new(),
            simple_internals: Vec::new(),
            symmetry_internals: Vec::new(),
            geom: Geom::new(),
            disps: Vec::new(),
            atoms: Vec::new(),
            fc2: Vec::new(),
            fc3: Vec::new(),
            fc4: Vec::new(),
        }
    }

    pub fn load_file(infile: &str) -> Self {
        let f = match File::open(infile) {
            Ok(f) => f,
            Err(_) => {
                eprintln!("failed to open infile '{infile}'");
                std::process::exit(1);
            }
        };
        Self::load(f)
    }

    /// helper function for parsing a single simple internal coordinate line
    fn parse_simple_internal(sp: Vec<&str>) -> Siic {
        let vals: Vec<_> = sp
            .iter()
            .skip(1)
            .map(|s| s.parse::<usize>().unwrap() - 1)
            .collect();
        match sp[0] {
            "STRE" => Siic::Stretch(vals[0], vals[1]),
            "BEND" => Siic::Bend(vals[0], vals[1], vals[2]),
            "TORS" => Siic::Torsion(vals[0], vals[1], vals[2], vals[3]),
            "LIN1" => Siic::Lin1(vals[0], vals[1], vals[2], vals[3]),
            "OUT" => Siic::Out(vals[0], vals[1], vals[2], vals[3]),
            "LINX" => Siic::Linx(vals[0], vals[1], vals[2], vals[3]),
            "LINY" => Siic::Liny(vals[0], vals[1], vals[2], vals[3]),
            e => {
                panic!("unknown coordinate type '{e}'");
            }
        }
    }

    pub fn load<R: Read>(r: R) -> Self {
        let siic = Regex::new(r"STRE|BEND|TORS|OUT|LIN|SPF|RCOM").unwrap();
        // something like "    1   1   1.000000000   2   1.000000000"
        let syic =
            Regex::new(r"^\s*(\d+\s+)(\d+\s+[0-9-]\.[0-9]+(\s+|$))+").unwrap();
        let zero = Regex::new(r"^\s*0\s*$").unwrap();
        // just a bunch of integers like 3    3    3    0    0    3 ..."
        let iops = Regex::new(r"^\s*(\d+(\s+|$))+$").unwrap();
        let geom = Regex::new(r"\s*([0-9-]+\.[0-9]+(\s+|$)){3}").unwrap();
        let atoms = Regex::new(r"[A-Za-z]+\d+").unwrap();
        let atom = Regex::new(r"([A-Za-z]+)(\d+)").unwrap();
        let fcs = Regex::new(r"^\s*(\d+\s+){4}[0-9-]+\.\d+\s*$").unwrap();

        let mut intder = Intder::new();
        let reader = BufReader::new(r);
        let mut in_disps = false;
        let mut disp_tmp = vec![];
        for line in reader.lines().map(Result::unwrap) {
            if line.contains("# INTDER #") || line.is_empty() {
                continue;
            } else if in_disps {
                // build up tmp until we hit a zero just like with syics
                if zero.is_match(&line) {
                    intder.disps.push(disp_tmp.clone());
                    disp_tmp = vec![0.0; intder.symmetry_internals.len()];
                    continue;
                }
                let sp: Vec<&str> = line.split_whitespace().collect();
                disp_tmp[sp[0].parse::<usize>().unwrap() - 1] = match sp.get(1)
                {
                    Some(s) => s,
                    None => panic!("line '{line}' too short"),
                }
                .parse::<f64>()
                .unwrap();
            } else if siic.is_match(&line) {
                intder.simple_internals.push(Self::parse_simple_internal(
                    line.split_whitespace().collect(),
                ));
            } else if syic.is_match(&line) {
                // this has to come after the simple internals
                assert!(!intder.simple_internals.is_empty());
                let mut tmp = vec![0.0; intder.simple_internals.len()];
                let mut sp = line.split_whitespace();
                sp.next();
                let mut idx = usize::default();
                let mut i_max = 0;
                for (i, c) in sp.enumerate() {
                    if i % 2 == 0 {
                        idx = c.parse::<usize>().unwrap() - 1;
                    } else {
                        tmp[idx] = c.parse().unwrap();
                    }
                    i_max += 1;
                }
                // i_max is 2 times the number of components in the SIC, if the
                // SIC is well-formed. for a single component (simple internal =
                // symmetry internal), i_max is 2, and so on. the prefactor for
                // the SICs is √n/n where n is the number of components
                assert!(i_max % 2 == 0);
                let nparts = (i_max / 2) as f64;
                if i_max > 2 {
                    let fac = nparts.sqrt() / nparts;
                    for t in &mut tmp {
                        *t *= fac;
                    }
                }
                intder.symmetry_internals.push(tmp);
            } else if zero.is_match(&line) {
                continue;
            } else if iops.is_match(&line) {
                intder.input_options.extend(
                    line.split_whitespace()
                        .map(|x| x.parse::<usize>().unwrap()),
                );
            } else if geom.is_match(&line) {
                let mut sp = line.split_ascii_whitespace().peekable();
                if let Some(s) = sp.peek() {
                    if let Some(c) = s.chars().next() {
                        if c.is_alphabetic() {
                            let atom = sp.next().unwrap();
                            intder
                                .atoms
                                .push(Atom::new(atom, DEFAULT_WEIGHTS[atom]));
                        }
                    }
                }
                if let [x, y, z] = sp
                    .map(|x| x.parse::<f64>().unwrap())
                    .collect::<Vec<f64>>()[..]
                {
                    intder.geom.push(na::Vector3::new(x, y, z));
                }
            } else if line.contains("DISP") {
                in_disps = true;
                disp_tmp = vec![0.0; intder.symmetry_internals.len()];
            } else if atoms.is_match(&line) {
                for cap in atom.captures_iter(&line) {
                    intder.atoms.push(Atom {
                        label: String::from(&cap[1]),
                        weight: cap[2].parse().unwrap(),
                    });
                }
            } else if fcs.is_match(&line) {
                let mut sp = line.split_whitespace().collect::<Vec<_>>();
                let val = sp.pop().unwrap().parse::<f64>().unwrap();
                let sp = sp
                    .iter()
                    .map(|s| s.parse::<usize>().unwrap())
                    .collect::<Vec<_>>();
                // here
                intder.add_fc(sp, val);
            }
        }
        intder
    }

    pub fn add_fc(&mut self, sp: Vec<usize>, val: f64) {
        let (idx, target) = match (sp[2], sp[3]) {
            (0, 0) => (
                fc2_index(self.symmetry_internals.len(), sp[0], sp[1]),
                &mut self.fc2,
            ),
            (_, 0) => (fc3_index(sp[0], sp[1], sp[2]), &mut self.fc3),
            (_, _) => (fc4_index(sp[0], sp[1], sp[2], sp[3]), &mut self.fc4),
        };
        if target.len() <= idx {
            target.resize(idx + 1, 0.0);
        }
        target[idx] = val;
    }

    /// return the number of symmetry internal coordinates
    pub fn nsym(&self) -> usize {
        self.symmetry_internals.len()
    }

    /// return the number of cartesian coordinates
    pub fn ncart(&self) -> usize {
        3 * self.geom.len()
    }

    /// return the number of dummy atoms
    pub fn ndum(&self) -> usize {
        self.input_options[7]
    }

    fn get_fc3(&self, i: usize, j: usize, k: usize) -> f64 {
        get_fc3(&self.fc3, i, j, k)
    }

    fn get_fc4(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        get_fc4(&self.fc4, i, j, k, l)
    }

    pub fn print_geom(&self) {
        for atom in &self.geom {
            for c in &atom {
                print!("{c:20.10}");
            }
            println!();
        }
    }

    /// print the simple internal coordinate values, assuming that vals are in
    /// the same order as self.simple_internals for unit purposes
    pub fn print_simple_values(&self, vals: &[f64]) {
        for (i, v) in vals.iter().enumerate() {
            if let Siic::Bend(_, _, _) = self.simple_internals[i] {
                println!("{:5}{:>18.10}", i, v.to_degrees());
            } else {
                println!("{i:5}{v:>18.10}");
            }
        }
    }

    /// print the symmetry internal coordinate values
    pub fn print_symmetry_values(&self, vals: &[f64]) {
        for (i, v) in vals.iter().enumerate() {
            println!("{i:5}{v:>18.10}");
        }
    }

    pub fn print_sics<W: std::io::Write>(&self, w: &mut W, irreps: &[Irrep]) {
        assert_eq!(self.symmetry_internals.len(), irreps.len());
        for (i, sic) in self.symmetry_internals.iter().enumerate() {
            write!(w, "S{i:<2}({}) = ", irreps[i]).unwrap();
            // number of siics printed so far
            let mut nprt = 0;
            for (j, s) in sic.iter().enumerate() {
                if *s != 0.0 {
                    if nprt > 0 {
                        let sign = match s.signum() as isize {
                            -1 => "-",
                            1 => "+",
                            _ => panic!("it's NaN"),
                        };
                        write!(w, " {sign} ").unwrap();
                    }
                    write!(w, "{}", &self.simple_internals[j]).unwrap();
                    nprt += 1;
                }
            }
            writeln!(w).unwrap();
        }
    }

    /// currently returns a vector of simple internal values in Ångstroms or
    /// radians
    pub fn simple_values(&self, geom: &Geom) -> Vec<f64> {
        self.simple_internals
            .iter()
            .map(|s| s.value(geom))
            .collect()
    }

    /// currently returns a vector of symmetry internal values in Ångstroms or
    /// radians
    pub fn symmetry_values(&self, geom: &Geom) -> Vec<f64> {
        let mut ret = Vec::with_capacity(self.symmetry_internals.len());
        let siics = self.simple_values(geom);
        for sic in &self.symmetry_internals {
            let mut sum = f64::default();
            for (i, s) in sic.iter().enumerate() {
                sum += s * siics[i];
            }
            ret.push(sum);
        }
        ret
    }

    /// return the B matrix in simple internal coordinates
    pub fn b_matrix(&self, geom: &Geom) -> DMat {
        let rows = self.simple_internals.len();
        let cols = 3 * geom.len();
        let mut b_mat = Vec::with_capacity(rows * cols);
        for ic in &self.simple_internals {
            b_mat.extend(geom.s_vec(ic));
        }
        DMat::from_row_slice(rows, cols, &b_mat)
    }

    /// return the U matrix, used for converting from simple internals to
    /// symmetry internals. dimensions are (number of symmetry internals) x
    /// (number of simple internals) since each symm. int. is a vector simp.
    /// int. long
    pub fn u_mat(&self) -> DMat {
        let r = self.symmetry_internals.len();
        let c = self.simple_internals.len();
        DMat::from_row_iterator(
            r,
            c,
            self.symmetry_internals.iter().flatten().cloned(),
        )
    }

    /// return the symmetry internal coordinate B matrix by computing the simple
    /// internal B and converting it
    pub fn sym_b_matrix(&self, geom: &Geom) -> DMat {
        let b = self.b_matrix(geom);
        self.u_mat() * b
    }

    /// Let D = BBᵀ and return A = BᵀD⁻¹
    pub fn a_matrix(b: &DMat) -> DMat {
        let d = b * b.transpose();
        match na::Cholesky::new(d.clone()) {
            Some(c) => b.transpose() * c.inverse(),
            None => {
                // compute it again to avoid cloning on the happy path
                eprintln!("{:.8}", b * b.transpose());
                eprintln!("cholesky decomposition failed");
                let c = na::LU::new(d);
                b.transpose()
                    * c.try_inverse().expect("lu decomposition failed")
            }
        }
    }

    /// print the initial geometry stuff
    fn print_init(&self) {
        let simple_vals = self.simple_values(&self.geom);
        let sic_vals = self.symmetry_values(&self.geom);
        println!();
        println!("NUCLEAR CARTESIAN COORDINATES (BOHR)\n");
        self.print_geom();
        println!();
        println!(
        "VALUES OF SIMPLE INTERNAL COORDINATES (ANG. or DEG.) FOR REFERENCE \
	     GEOMETRY\n"
    );
        self.print_simple_values(&simple_vals);
        println!();
        println!(
        "VALUES OF SYMMETRY INTERNAL COORDINATES (ANG. or RAD.) FOR REFERENCE \
	     GEOMETRY\n"
    );
        self.print_symmetry_values(&sic_vals);
        println!();
        println!();
    }

    pub fn print_cart<W: Write>(w: &mut W, cart: &DVec) {
        for i in 0..cart.len() / 3 {
            for j in 0..3 {
                write!(w, "{:20.10}", cart[3 * i + j]).unwrap();
            }
            writeln!(w).unwrap();
        }
    }

    /// convert the displacements in `self.disps` from (symmetry) internal
    /// coordinates to Cartesian coordinates
    pub fn convert_disps(&self) -> Result<Vec<DVec>, IntderError> {
        if is_verbose() {
            self.print_init();
        }
        let sic0 = DVec::from(self.symmetry_values(&self.geom));
        let cart0: DVec = self.geom.clone().into();
        let mut ret = Vec::new();
        for (i, disp) in self.disps.iter().enumerate() {
            // initialize sics and carts to those from the input file
            let mut sic_current = sic0.clone();
            let mut cart_current: DVec = cart0.clone();

            // get a vector from the displacement from the input file
            let disp = DVec::from(disp.clone());
            let sic_desired = &sic_current + &disp;

            if is_verbose() {
                println!("DISPLACEMENT{i:5}\n");
                println!("INTERNAL DISPLACEMENTS\n");
                for (i, d) in disp.iter().enumerate() {
                    if *d != 0.0 {
                        println!("{i:5}{d:20.10}");
                    }
                }
                println!();
                println!("SYMMETRY INTERNAL COORDINATE FINAL VALUES\n");
                self.print_symmetry_values(sic_desired.as_slice());
                println!();
            }

            // measure convergence by max internal deviation between the current
            // SICs and desired SICs
            let mut iter = 1;
            while (&sic_current - &sic_desired).abs().max() > TOLDISP {
                let b_sym = self.sym_b_matrix(&Geom::from(&cart_current));
                let d = &b_sym * b_sym.transpose();
                let a = Intder::a_matrix(&b_sym);

                if is_verbose() {
                    println!(
                        "ITER={:5} MAX INTERNAL DEVIATION = {:.4e}",
                        iter,
                        (&sic_current - &sic_desired).abs().max()
                    );
                    println!("B*BT MATRIX FOR (SYMMETRY) INTERNAL COORDINATES");
                    println!("{d:.6}");

                    println!(
                        "DETERMINANT OF B*BT MATRIX={:8.4}",
                        d.determinant()
                    );

                    println!();
                    println!("A MATRIX FOR (SYMMETRY) INTERNAL COORDINATES");
                    println!("{a:.8}");
                    println!();
                }

                let step = a * (&sic_desired - &sic_current);
                cart_current += step / ANGBOHR;

                sic_current = DVec::from(
                    self.symmetry_values(&Geom::from(&cart_current)),
                );

                iter += 1;

                if MAX_ITER > 0 && iter > MAX_ITER {
                    return Err(IntderError::DispError(format!(
                        "max iterations exceeded on disp {i}"
                    )));
                }
            }

            if is_verbose() {
                println!(
                    "ITER={:5} MAX INTERNAL DEVIATION = {:.4e}\n",
                    iter,
                    (&sic_current - &sic_desired).abs().max()
                );
                println!("NEW CARTESIAN GEOMETRY (BOHR)\n");
                Self::print_cart(&mut std::io::stdout(), &cart_current);
                println!();
            }
            ret.push(cart_current);
        }
        Ok(ret)
    }

    /// return the SR matrix in symmetry internal coordinates
    pub fn machx(&self) -> Vec<DMat> {
        use Siic::*;
        let nc = self.ncart();
        let nsym = self.nsym();
        let u = self.u_mat();
        if nsym == 0 {
            todo!("using only simple internals is unimplemented");
        }
        // simple internal X and SR matrices
        let mut srs_sim = Vec::new();
        // let nsim = self.simple_internals.len();
        for s in &self.simple_internals {
            let mut sr = DMat::zeros(nc, nc);
            let h = Hmat::new(&self.geom, s);
            let size = (3, 3);
            match s {
                Stretch(a, b) => {
                    let l1 = 3 * a;
                    let l2 = 3 * b;
                    sr.view_mut((l1, l1), size).copy_from(&h.h11);
                    sr.view_mut((l2, l2), size).copy_from(&h.h11);
                    sr.view_mut((l1, l2), size).copy_from(&-&h.h11);
                    sr.view_mut((l2, l1), size).copy_from(&-&h.h11);
                }
                Bend(a, b, c) | Lin1(a, b, c, _) => {
                    let l1 = 3 * a;
                    let l2 = 3 * b;
                    let l3 = 3 * c;
                    sr.view_mut((l1, l1), size).copy_from(&h.h11);
                    sr.view_mut((l2, l1), size).copy_from(&h.h21);
                    sr.view_mut((l3, l1), size).copy_from(&h.h31);
                    sr.view_mut((l1, l2), size).copy_from(&h.h21.transpose());
                    sr.view_mut((l2, l2), size).copy_from(&h.h22);
                    sr.view_mut((l3, l2), size).copy_from(&h.h32);
                    sr.view_mut((l1, l3), size).copy_from(&h.h31.transpose());
                    sr.view_mut((l2, l3), size).copy_from(&h.h32.transpose());
                    sr.view_mut((l3, l3), size).copy_from(&h.h33);
                }
                Torsion(a, b, c, d)
                | Out(a, b, c, d)
                | Linx(a, b, c, d)
                | Liny(a, b, c, d) => {
                    let l1 = 3 * a;
                    let l2 = 3 * b;
                    let l3 = 3 * c;
                    let l4 = 3 * d;
                    sr.view_mut((l1, l1), size).copy_from(&h.h11);
                    sr.view_mut((l2, l1), size).copy_from(&h.h21);
                    sr.view_mut((l3, l1), size).copy_from(&h.h31);
                    sr.view_mut((l4, l1), size).copy_from(&h.h41);
                    sr.view_mut((l1, l2), size).copy_from(&h.h21.transpose());
                    sr.view_mut((l2, l2), size).copy_from(&h.h22);
                    sr.view_mut((l3, l2), size).copy_from(&h.h32);
                    sr.view_mut((l4, l2), size).copy_from(&h.h42);
                    sr.view_mut((l1, l3), size).copy_from(&h.h31.transpose());
                    sr.view_mut((l2, l3), size).copy_from(&h.h32.transpose());
                    sr.view_mut((l3, l3), size).copy_from(&h.h33);
                    sr.view_mut((l4, l3), size).copy_from(&h.h43);
                    sr.view_mut((l1, l4), size).copy_from(&h.h41.transpose());
                    sr.view_mut((l2, l4), size).copy_from(&h.h42.transpose());
                    sr.view_mut((l3, l4), size).copy_from(&h.h43.transpose());
                    sr.view_mut((l4, l4), size).copy_from(&h.h44);
                }
            }
            srs_sim.push(sr);
        }
        // TODO if nsym = 0, just return the sim versions
        let mut srs_sym = Vec::with_capacity(nsym);
        for r in 0..nsym {
            let mut sr_sic = DMat::zeros(nc, nc);
            for (i, sr) in srs_sim.iter().enumerate() {
                for n in 0..nc {
                    for m in 0..nc {
                        sr_sic[(m, n)] += u[(r, i)] * sr[(m, n)];
                    }
                }
            }
            srs_sym.push(sr_sic);
        }
        srs_sym
    }

    /// returns the SR matrix in symmetry internal coordinates
    pub fn machy(&self) -> Vec<Tensor3> {
        use Siic::*;
        let nc = self.ncart();
        let nsx = self.nsym();
        if nsx == 0 {
            todo!("using only simple internals is unimplemented");
        }
        let u = self.u_mat();
        let mut srs_sim = Vec::new();
        for s in &self.simple_internals {
            let mut sr = Tensor3::zeros((nc, nc, nc));
            let h = Htens::new(&self.geom, s);
            match s {
                Stretch(a, b) => {
                    hsry2(&mut sr, 3 * a, 3 * b, &h);
                }
                Bend(a, b, c) | Lin1(a, b, c, _) => {
                    hsry3(&mut sr, 3 * a, 3 * b, 3 * c, &h);
                }
                Torsion(a, b, c, d)
                | Out(a, b, c, d)
                | Linx(a, b, c, d)
                | Liny(a, b, c, d) => {
                    hsry4(&mut sr, 3 * a, 3 * b, 3 * c, 3 * d, &h);
                }
            }
            srs_sim.push(sr);
        }

        // convert SR to symmetry internals
        let mut ret = Vec::with_capacity(nsx);
        for r in 0..nsx {
            let mut sr_sic = Tensor3::zeros((nc, nc, nc));
            for (i, sr) in srs_sim.iter().enumerate() {
                let w1 = u[(r, i)];
                for p in 0..nc {
                    for n in 0..=p {
                        for m in 0..=n {
                            sr_sic[(m, n, p)] += w1 * sr[(m, n, p)];
                        }
                    }
                }
            }
            fill3a(&mut sr_sic, nc);
            ret.push(sr_sic);
        }
        ret
    }

    /// represent fc2 as a symmetric matrix
    fn mat_fc2(&self, nsy: usize) -> DMat {
        let mut v = DMat::from_row_slice(nsy, nsy, &self.fc2);
        for row in 0..nsy {
            for col in row..nsy {
                v[(col, row)] = v[(row, col)];
            }
        }
        v
    }

    fn lintr_fc2(&self, bs: &DMat) -> DMat {
        let nsx = self.ncart() - 3 * self.ndum();
        let nsy = self.symmetry_internals.len();
        let f2 = bs.transpose() * self.mat_fc2(nsy) * bs;
        f2.resize(nsx, nsx, 0.0) * ANGBOHR * ANGBOHR / HART
    }

    pub fn lintr_fc3(&self, bs: &DMat) -> Tensor3 {
        let nsx = self.ncart() - 3 * self.ndum();
        let nsy = self.nsym();

        let mut f3 = Tensor3::zeros((nsx, nsx, nsx));
        for i in 0..nsy {
            for j in 0..=i {
                let dij = delta(i, j);
                for k in 0..=j {
                    let djk = delta(j, k);
                    let vik = self.get_fc3(i, j, k);
                    for p in 0..nsx {
                        f3[(i, j, p)] += vik * bs[(k, p)];
                        f3[(i, k, p)] += vik * bs[(j, p)] * djk;
                        f3[(j, k, p)] += vik * bs[(i, p)] * dij;
                    }
                }
            }
        }
        // end of 1138 loop, looking good so far

        let f3_disk = f3;
        let mut f3 = Tensor3::zeros((nsx, nsx, nsx));
        for i in 0..nsy {
            for j in 0..=i {
                let dij = delta(i, j);
                for p in 0..nsx {
                    let vik = f3_disk[(i, j, p)];
                    for n in 0..=p {
                        f3[(i, n, p)] += vik * bs[(j, n)];
                        f3[(j, n, p)] += vik * bs[(i, n)] * dij;
                    }
                }
            }
        }
        // end of 1146 loop, looking good again

        let f3_disk2 = f3;
        let mut f3 = Tensor3::zeros((nsx, nsx, nsx));
        for i in 0..nsy {
            for p in 0..nsx {
                for n in 0..=p {
                    let vik = f3_disk2[(i, n, p)];
                    for m in 0..=n {
                        f3[(m, n, p)] += vik * bs[(i, m)];
                    }
                }
            }
        }
        // end of 1152 loop, still looking good

        fill3a(&mut f3, nsx);
        f3
    }

    pub fn lintr_fc4(&self, bs: &DMat) -> Tensor4 {
        let nsx = self.ncart() - 3 * self.ndum();
        let nsy = self.nsym();
        let mut f4 = Tensor4::zeros((nsx, nsx, nsx, nsx));
        for i in 0..nsy {
            for j in 0..=i {
                let dij = delta(i, j);
                for k in 0..=j {
                    let djk = delta(j, k);
                    for l in 0..=k {
                        let dkl = delta(k, l);
                        let vik = self.get_fc4(i, j, k, l);
                        for q in 0..nsx {
                            f4[(i, j, k, q)] += vik * bs[(l, q)];
                            f4[(i, j, l, q)] += vik * bs[(k, q)] * dkl;
                            f4[(i, k, l, q)] += vik * bs[(j, q)] * djk;
                            f4[(j, k, l, q)] += vik * bs[(i, q)] * dij;
                        }
                    }
                }
            }
        }
        // end 179 loop, looking good so far

        let f4_disk = f4;
        let mut f4 = Tensor4::zeros((nsx, nsx, nsx, nsx));

        // loop starting at line 444
        for i in 0..nsy {
            for j in 0..=i {
                let dij = delta(i, j);
                for k in 0..=j {
                    let djk = delta(j, k);
                    for q in 0..nsx {
                        let vik = f4_disk[(i, j, k, q)];
                        for p in 0..=q {
                            f4[(i, j, p, q)] += vik * bs[(k, p)];
                            f4[(i, k, p, q)] += vik * bs[(j, p)] * djk;
                            f4[(j, k, p, q)] += vik * bs[(i, p)] * dij;
                        }
                    }
                }
            }
        }
        // end 200 loop

        let f4_disk2 = f4;
        let mut f4 = Tensor4::zeros((nsx, nsx, nsx, nsx));

        // start of loop at 514
        for i in 0..nsy {
            for j in 0..=i {
                let dij = delta(i, j);
                for q in 0..nsx {
                    for p in 0..=q {
                        let vik = f4_disk2[(i, j, p, q)];
                        for n in 0..=p {
                            f4[(i, n, p, q)] += vik * bs[(j, n)];
                            f4[(j, n, p, q)] += vik * bs[(i, n)] * dij;
                        }
                    }
                }
            }
        }
        // end of 214 loop, looking good

        let f4_disk2 = f4;
        let mut f4 = Tensor4::zeros((nsx, nsx, nsx, nsx));

        // begin 224 loop
        for i in 0..nsy {
            for q in 0..nsx {
                for p in 0..=q {
                    for n in 0..=p {
                        let vik = f4_disk2[(i, n, p, q)];
                        for m in 0..=n {
                            f4[(m, n, p, q)] += vik * bs[(i, m)];
                        }
                    }
                }
            }
        }
        // end 224 loop

        fill4a(&mut f4, nsx);
        f4
    }

    /// harmonic contributions of S_x to the cubic and quartic force constants
    fn xf2(
        &self,
        mut f3: Tensor3,
        mut f4: Tensor4,
        bs: &DMat,
        xrs: &[DMat],
    ) -> (Tensor3, Tensor4) {
        let ns = self.symmetry_internals.len();
        let nc = self.ncart() - 3 * self.ndum();
        let v = self.mat_fc2(ns);
        let xs = v * bs;
        for (r, xr) in xrs.iter().enumerate() {
            for k in 0..nc {
                for j in 0..=k {
                    for i in 0..=j {
                        let w = xr[(i, j)] * xs[(r, k)]
                            + xr[(i, k)] * xs[(r, j)]
                            + xr[(j, k)] * xs[(r, i)];
                        f3[(i, j, k)] += w;
                    }
                }
            }
        }
        fill3a(&mut f3, nc);

        // begin F4 part

        // this should never be used as zeros. it will be set on the first
        // iteration when r == s == 0
        let v = self.mat_fc2(ns);
        for (r, xt) in xrs.iter().enumerate() {
            let mut xr = DMat::zeros(nc, nc);
            for (s, xs) in xrs.iter().enumerate() {
                for j in 0..nc {
                    for i in 0..nc {
                        xr[(i, j)] += v[(r, s)] * xs[(i, j)];
                    }
                }
            }
            for l in 0..nc {
                for k in 0..=l {
                    for j in 0..=k {
                        for i in 0..=j {
                            let w = xt[(i, j)] * xr[(k, l)]
                                + xt[(i, k)] * xr[(j, l)]
                                + xt[(i, l)] * xr[(j, k)];
                            f4[(i, j, k, l)] += w;
                        }
                    }
                }
            }
        }
        fill4a(&mut f4, nc);
        (f3, f4)
    }

    /// cubic contributions of S_x to the quartic force constants
    fn xf3(&self, mut f4: Tensor4, bs: &DMat, xrs: &[DMat]) -> Tensor4 {
        let ns = self.nsym();
        let nc = self.ncart() - 3 * self.ndum();
        for (r, xr) in xrs.iter().enumerate() {
            let mut xt = DMat::zeros(ns, nc);
            for l in 0..nc {
                for n in 0..ns {
                    for p in 0..ns {
                        let x = self.get_fc3(r, n, p);
                        xt[(n, l)] += x * bs[(p, l)]
                    }
                }
            }
            let xs = bs.transpose() * xt;
            for l in 0..nc {
                for k in 0..=l {
                    for j in 0..=k {
                        for i in 0..=j {
                            let w = xr[(i, j)] * xs[(k, l)]
                                + xr[(i, k)] * xs[(j, l)]
                                + xr[(j, k)] * xs[(i, l)]
                                + xr[(k, l)] * xs[(i, j)]
                                + xr[(j, l)] * xs[(i, k)]
                                + xr[(i, l)] * xs[(j, k)];
                            f4[(i, j, k, l)] += w;
                        }
                    }
                }
            }
        }
        fill4a(&mut f4, nc);
        f4
    }

    /// harmonic contributions of S_y to the quartic force constants
    fn yf2(&self, mut f4: Tensor4, bs: &DMat, yrs: &[Tensor3]) -> Tensor4 {
        let nc = self.ncart() - 3 * self.ndum();
        let ns = self.nsym();
        let xs = self.mat_fc2(ns);
        let xr = xs * bs;
        for (r, yr) in yrs.iter().enumerate() {
            for l in 0..nc {
                for k in 0..=l {
                    for j in 0..=k {
                        for i in 0..=j {
                            let w = yr[(i, j, k)] * xr[(r, l)]
                                + yr[(i, j, l)] * xr[(r, k)]
                                + yr[(i, k, l)] * xr[(r, j)]
                                + yr[(j, k, l)] * xr[(r, i)];
                            f4[(i, j, k, l)] += w;
                        }
                    }
                }
            }
        }
        fill4a(&mut f4, nc);
        f4
    }

    /// Perform the linear transformation of the force constants and convert the
    /// units to those desired by SPECTRO.
    pub fn lintr(
        &self,
        bs: &DMat,
        xr: &[DMat],
        yr: &[Tensor3],
    ) -> (DMat, Vec<f64>, Vec<f64>) {
        let f2 = self.lintr_fc2(bs);
        let (f3, f4) = self.xf2(self.lintr_fc3(bs), self.lintr_fc4(bs), bs, xr);
        let f4 = self.xf3(f4, bs, xr);
        let f4 = self.yf2(f4, bs, yr);

        // convert f3 to the proper units and return it as a Vec in the order
        // desired by spectro
        let nsx = self.ncart() - 3 * self.ndum();
        const CF3: f64 = ANGBOHR * ANGBOHR * ANGBOHR / HART;
        let mut f3_out = Vec::new();
        for m in 0..nsx {
            for n in 0..=m {
                for p in 0..=n {
                    f3_out.push(CF3 * f3[(m, n, p)]);
                }
            }
        }

        // convert f4 to the proper units and return it as a Vec in the order
        // desired by spectro
        const CF4: f64 = ANGBOHR * ANGBOHR * ANGBOHR * ANGBOHR / HART;
        let mut f4_out = Vec::new();
        for m in 0..nsx {
            for n in 0..=m {
                for p in 0..=n {
                    for q in 0..=p {
                        f4_out.push(CF4 * f4[(m, n, p, q)]);
                    }
                }
            }
        }
        (f2, f3_out, f4_out)
    }

    /// convert the force constants in `self.fc[234]` from (symmetry) internal
    /// coordinates to Cartesian coordinates. returns (fc2, fc3, fc4) in the
    /// order printed in the fort.{15,30,40} files for spectro.
    pub fn convert_fcs(&self) -> (DMat, Vec<f64>, Vec<f64>) {
        if is_verbose() {
            self.print_init();
        }
        // let sics = DVec::from(self.symmetry_values(&self.geom));
        let b_sym = self.sym_b_matrix(&self.geom);
        let srs = self.machx();
        let srsy = self.machy();
        let (f2, f3, f4) = self.lintr(&b_sym, &srs, &srsy);

        (f2, f3, f4)
    }

    pub fn dump_fcs(dir: &str, f2: &DMat, f3: &[f64], f4: &[f64]) {
        let f2 = f2.as_slice();
        let pairs = [(f2, "fort.15"), (f3, "fort.30"), (f4, "fort.40")];
        for (fcs, file) in pairs {
            let filename = format!("{dir}/{file}");
            let mut f = match File::create(&filename) {
                Ok(f) => f,
                Err(e) => panic!("failed to create `{filename}` for `{e}`"),
            };
            for chunk in fcs.chunks(3) {
                for c in chunk {
                    write!(f, " {c:>19.10}").unwrap();
                }
                writeln!(f).unwrap();
            }
        }
    }

    /// detect the dummy atoms needed in `self` as the atoms extending from the
    /// [LIN1]s and add two dummy atoms per LIN1, in the two directions
    /// perpendicular to axis. I think this still assumes the molecule is linear
    /// because it takes the `axis` coordinate of the real geometry and combines
    /// that with 1.111111111 and 0.0 as the coordinate of the dummy atom. This
    /// might work because this section of the molecule probably has to be on an
    /// axis to use a LIN1, but it could be wrong. We could take the non-1.11111
    /// coordinate from the real geometry as well
    pub fn add_dummies(&mut self, axis: Axis) -> usize {
        let dummies: HashSet<usize> = self
            .simple_internals
            .iter()
            .filter_map(|sic| {
                if let Siic::Lin1(_, b, _, _) = sic {
                    Some(*b)
                } else {
                    None
                }
            })
            .collect();

        if dummies.is_empty() {
            return 0;
        }

        // add the dummy atoms
        let mut ndum = 0;
        for dummy in dummies {
            // atom the dummy attaches to
            let real_coord = self.geom[dummy];
            let nonzero = axis as usize;
            let (a, b) = axis.not();
            let mut coord = [0.0; 3];
            // match the nonzero field in the real geometry
            coord[nonzero] = real_coord[nonzero];
            coord[a as usize] = 1.1111111111;
            coord[b as usize] = 0.0;
            self.geom.push(na::Vector3::from(coord));

            let mut coord = [0.0; 3];
            // match the nonzero field in the real geometry
            coord[nonzero] = real_coord[nonzero];
            coord[b as usize] = 1.1111111111;
            coord[a as usize] = 0.0;
            self.geom.push(na::Vector3::from(coord));

            // push dummy atoms perpendicular in both directions

            ndum += 2;
        }
        self.input_options[7] = ndum;
        ndum
    }
}

fn hsry2(sr: &mut Tensor3, l1: usize, l2: usize, h: &Htens) {
    for k in 0..3 {
        for j in 0..3 {
            for i in 0..3 {
                let z = h.h111[(i, j, k)];
                sr[(l1 + i, l1 + j, l1 + k)] = z;
                sr[(l1 + i, l2 + j, l2 + k)] = z;
                sr[(l2 + i, l1 + j, l2 + k)] = z;
                sr[(l2 + i, l2 + j, l1 + k)] = z;
                sr[(l1 + i, l1 + j, l2 + k)] = -z;
                sr[(l1 + i, l2 + j, l1 + k)] = -z;
                sr[(l2 + i, l2 + j, l2 + k)] = -z;
                sr[(l2 + i, l1 + j, l1 + k)] = -z;
            }
        }
    }
}

fn hsry3(sr: &mut Tensor3, l1: usize, l2: usize, l3: usize, h: &Htens) {
    for k in 0..3 {
        for j in 0..3 {
            for i in 0..3 {
                sr[(l1 + i, l1 + j, l1 + k)] = h.h111[(i, j, k)];
                sr[(l1 + i, l1 + j, l2 + k)] = h.h112[(i, j, k)];
                sr[(l1 + i, l1 + j, l3 + k)] = h.h113[(i, j, k)];
                sr[(l1 + i, l2 + j, l1 + k)] = h.h112[(i, k, j)];
                sr[(l1 + i, l2 + j, l2 + k)] = h.h221[(j, k, i)];
                sr[(l1 + i, l2 + j, l3 + k)] = h.h123[(i, j, k)];
                sr[(l2 + i, l2 + j, l1 + k)] = h.h221[(i, j, k)];
                sr[(l2 + i, l2 + j, l2 + k)] = h.h222[(i, j, k)];
                sr[(l2 + i, l2 + j, l3 + k)] = h.h223[(i, j, k)];
                sr[(l2 + i, l1 + j, l1 + k)] = h.h112[(j, k, i)];
                sr[(l2 + i, l1 + j, l2 + k)] = h.h221[(i, k, j)];
                sr[(l2 + i, l1 + j, l3 + k)] = h.h123[(j, i, k)];
                sr[(l1 + i, l3 + j, l1 + k)] = h.h113[(i, k, j)];
                sr[(l1 + i, l3 + j, l2 + k)] = h.h123[(i, k, j)];
                sr[(l1 + i, l3 + j, l3 + k)] = h.h331[(j, k, i)];
                sr[(l2 + i, l3 + j, l1 + k)] = h.h123[(k, i, j)];
                sr[(l2 + i, l3 + j, l2 + k)] = h.h223[(i, k, j)];
                sr[(l2 + i, l3 + j, l3 + k)] = h.h332[(j, k, i)];
                sr[(l3 + i, l1 + j, l1 + k)] = h.h113[(j, k, i)];
                sr[(l3 + i, l1 + j, l2 + k)] = h.h123[(j, k, i)];
                sr[(l3 + i, l1 + j, l3 + k)] = h.h331[(i, k, j)];
                sr[(l3 + i, l2 + j, l1 + k)] = h.h123[(k, j, i)];
                sr[(l3 + i, l2 + j, l2 + k)] = h.h223[(j, k, i)];
                sr[(l3 + i, l2 + j, l3 + k)] = h.h332[(i, k, j)];
                sr[(l3 + i, l3 + j, l1 + k)] = h.h331[(i, j, k)];
                sr[(l3 + i, l3 + j, l2 + k)] = h.h332[(i, j, k)];
                sr[(l3 + i, l3 + j, l3 + k)] = h.h333[(i, j, k)];
            }
        }
    }
}

fn hsry4(
    sr: &mut Tensor3,
    l1: usize,
    l2: usize,
    l3: usize,
    l4: usize,
    h: &Htens,
) {
    for k in 0..3 {
        for j in 0..3 {
            for i in 0..3 {
                sr[(l1 + i, l1 + j, l1 + k)] = h.h111[(i, j, k)];
                sr[(l1 + i, l1 + j, l2 + k)] = h.h112[(i, j, k)];
                sr[(l1 + i, l1 + j, l3 + k)] = h.h113[(i, j, k)];
                sr[(l1 + i, l1 + j, l4 + k)] = h.h411[(k, j, i)];
                sr[(l1 + i, l2 + j, l1 + k)] = h.h112[(i, k, j)];
                sr[(l1 + i, l2 + j, l2 + k)] = h.h221[(j, k, i)];
                sr[(l1 + i, l2 + j, l3 + k)] = h.h123[(i, j, k)];
                sr[(l1 + i, l2 + j, l4 + k)] = h.h421[(k, j, i)];
                sr[(l1 + i, l3 + j, l1 + k)] = h.h113[(i, k, j)];
                sr[(l1 + i, l3 + j, l2 + k)] = h.h123[(i, k, j)];
                sr[(l1 + i, l3 + j, l3 + k)] = h.h331[(j, k, i)];
                sr[(l1 + i, l3 + j, l4 + k)] = h.h431[(k, j, i)];
                sr[(l1 + i, l4 + j, l1 + k)] = h.h411[(j, k, i)];
                sr[(l1 + i, l4 + j, l2 + k)] = h.h421[(j, k, i)];
                sr[(l1 + i, l4 + j, l3 + k)] = h.h431[(j, k, i)];
                sr[(l1 + i, l4 + j, l4 + k)] = h.h441[(j, k, i)];
                sr[(l2 + i, l1 + j, l1 + k)] = h.h112[(j, k, i)];
                sr[(l2 + i, l1 + j, l2 + k)] = h.h221[(i, k, j)];
                sr[(l2 + i, l1 + j, l3 + k)] = h.h123[(j, i, k)];
                sr[(l2 + i, l1 + j, l4 + k)] = h.h421[(k, i, j)];
                sr[(l2 + i, l2 + j, l1 + k)] = h.h221[(i, j, k)];
                sr[(l2 + i, l2 + j, l2 + k)] = h.h222[(i, j, k)];
                sr[(l2 + i, l2 + j, l3 + k)] = h.h223[(i, j, k)];
                sr[(l2 + i, l2 + j, l4 + k)] = h.h422[(k, j, i)];
                sr[(l2 + i, l3 + j, l1 + k)] = h.h123[(k, i, j)];
                sr[(l2 + i, l3 + j, l2 + k)] = h.h223[(i, k, j)];
                sr[(l2 + i, l3 + j, l3 + k)] = h.h332[(j, k, i)];
                sr[(l2 + i, l3 + j, l4 + k)] = h.h432[(k, j, i)];
                sr[(l2 + i, l4 + j, l1 + k)] = h.h421[(j, i, k)];
                sr[(l2 + i, l4 + j, l2 + k)] = h.h422[(j, i, k)];
                sr[(l2 + i, l4 + j, l3 + k)] = h.h432[(j, k, i)];
                sr[(l2 + i, l4 + j, l4 + k)] = h.h442[(j, k, i)];
                sr[(l3 + i, l1 + j, l1 + k)] = h.h113[(j, k, i)];
                sr[(l3 + i, l1 + j, l2 + k)] = h.h123[(j, k, i)];
                sr[(l3 + i, l1 + j, l3 + k)] = h.h331[(i, k, j)];
                sr[(l3 + i, l1 + j, l4 + k)] = h.h431[(k, i, j)];
                sr[(l3 + i, l2 + j, l1 + k)] = h.h123[(k, j, i)];
                sr[(l3 + i, l2 + j, l2 + k)] = h.h223[(j, k, i)];
                sr[(l3 + i, l2 + j, l3 + k)] = h.h332[(i, k, j)];
                sr[(l3 + i, l2 + j, l4 + k)] = h.h432[(k, i, j)];
                sr[(l3 + i, l3 + j, l1 + k)] = h.h331[(i, j, k)];
                sr[(l3 + i, l3 + j, l2 + k)] = h.h332[(i, j, k)];
                sr[(l3 + i, l3 + j, l3 + k)] = h.h333[(i, j, k)];
                sr[(l3 + i, l3 + j, l4 + k)] = h.h433[(k, i, j)];
                sr[(l3 + i, l4 + j, l1 + k)] = h.h431[(j, i, k)];
                sr[(l3 + i, l4 + j, l2 + k)] = h.h432[(j, i, k)];
                sr[(l3 + i, l4 + j, l3 + k)] = h.h433[(j, i, k)];
                sr[(l3 + i, l4 + j, l4 + k)] = h.h443[(j, k, i)];
                sr[(l4 + i, l1 + j, l1 + k)] = h.h411[(i, j, k)];
                sr[(l4 + i, l1 + j, l2 + k)] = h.h421[(i, k, j)];
                sr[(l4 + i, l1 + j, l3 + k)] = h.h431[(i, k, j)];
                sr[(l4 + i, l1 + j, l4 + k)] = h.h441[(i, k, j)];
                sr[(l4 + i, l2 + j, l1 + k)] = h.h421[(i, j, k)];
                sr[(l4 + i, l2 + j, l2 + k)] = h.h422[(i, j, k)];
                sr[(l4 + i, l2 + j, l3 + k)] = h.h432[(i, k, j)];
                sr[(l4 + i, l2 + j, l4 + k)] = h.h442[(i, k, j)];
                sr[(l4 + i, l3 + j, l1 + k)] = h.h431[(i, j, k)];
                sr[(l4 + i, l3 + j, l2 + k)] = h.h432[(i, j, k)];
                sr[(l4 + i, l3 + j, l3 + k)] = h.h433[(i, j, k)];
                sr[(l4 + i, l3 + j, l4 + k)] = h.h443[(i, k, j)];
                sr[(l4 + i, l4 + j, l1 + k)] = h.h441[(i, j, k)];
                sr[(l4 + i, l4 + j, l2 + k)] = h.h442[(i, j, k)];
                sr[(l4 + i, l4 + j, l3 + k)] = h.h443[(i, j, k)];
                sr[(l4 + i, l4 + j, l4 + k)] = h.h444[(i, j, k)];
            }
        }
    }
}

// NOTE: these functions were extracted automatically by the LSP, so their
// interface could probably be cleaned up a bit if desired
