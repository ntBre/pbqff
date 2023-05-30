[![check](https://github.com/ntBre/pbqff/actions/workflows/check.yml/badge.svg)](https://github.com/ntBre/pbqff/actions/workflows/check.yml)
[![test](https://github.com/ntBre/pbqff/actions/workflows/test.yml/badge.svg)](https://github.com/ntBre/pbqff/actions/workflows/test.yml)

# pbqff
pbqff automates the construction of quartic force fields (QFFs) and the process
of generating [spectral data](https://github.com/ntBre/spectro) from them

# Installation

Assuming you have the [Rust toolchain](https://www.rust-lang.org/tools/install)
installed, run

```bash
make install
```

As you can see in the Makefile, this simply runs

```bash
cargo build --features vers --bin rust-pbqff --release
sudo ln -s target/release/rust-pbqff /usr/bin/pbqff
```

to build the binary in release mode, and link it into your `$PATH` under the
name `pbqff`.

## Dependencies

If you're installing `pbqff` on a "normal" machine, you're very likely to have
most of these programs already. But if you install on a fresh, minimal Ubuntu
installation (like I do in [this
video](https://www.youtube.com/watch?v=y-FH-LBqzXM)), you might need to install
some or all of these:

- curl (for rustup install)
- rust nightly toolchain
- git
- make
- openssl (on arch) or libssl-dev (on ubuntu) (for vers feature)
- pkg-config (for locating openssl)
- gcc (for linking)
- python3, python3-tk, idle3 (for qffbuddy)
- cmake, gfortran, libblas-dev, liblapack-dev (for MOPAC)

For installing or building `pbqff` itself, you can skip this last set, which are
required for building MOPAC from source. However, if you want to run the tests
for `pbqff`, you will need MOPAC installed at `/opt/mopac/mopac`, so these
dependencies are necessary in that case.

# Coordinate Types
pbqff supports running QFFs in the following coordinate systems:
- Symmetry-internal coordinates (SICs) via
  [intder](https://github.com/ntBre/intder)
- Cartesian coordinates
- Normal coordinates

The normal coordinates are determined automatically by running a Cartesian
harmonic force field.

# Programs and Queues
pbqff supports Molpro and Mopac for computing single-point energies and the PBS
and Slurm queuing systems via [psqs](https://github.com/ntBre/psqs)

# Input
An example input file for a Mopac QFF on *c*-C<sub>3</sub>H<sub>2</sub> looks like:
```
geometry = """
C
C 1 CC
C 1 CC 2 CCC
H 2 CH 1 HCC 3 180.0
H 3 CH 1 HCC 2 180.0

CC =                  1.42101898
CCC =                55.60133141
CH =                  1.07692776
HCC =               147.81488230
"""
optimize = true
charge = 0
step_size = 0.005
coord_type = "sic"
program = "mopac"
queue = "slurm"
sleep_int = 2
job_limit = 2048
chunk_size = 1
template = """scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1 \
external=testfiles/params.dat"""
check_int = 100%
```

# qffbuddy
An optional GUI for preparing input files and running pbqff is also included in
the qffbuddy directory

# Citations
For `pbqff` itself, please cite B. R. Westbrook and R. C. Fortenberry. "pbqff:
Push-Button Quartic Force Fields." *J. Chem. Theory Comput.*, 2023. DOI:
[10.1021/acs.jctc.3c00129](https://doi.org/10.1021/acs.jctc.3c00129)

```
@article{Westbrook23_pbqff,
  author = {Brent R. Westbrook and Ryan C. Fortenberry},
  title = {pbqff: Push-Button Quartic Force Fields},
  journal = {J. Chem. Theory Comput.},
  volume = 19,
  number = 9,
  pages = {2606-2615},
  year = 2023,
}
```

If you use symmetry-internal coordinates, you may also want to cite the original
`INTDER` code by Wesley Allen:

```
@misc{intder,
 author = {W. D. Allen and coworkers},
 note = {$INTDER\ 2005$ is a General Program Written by W. D. Allen and Coworkers, which Performs Vibrational Analysis and Higher-Order Non-Linear Transformations.},
 year = {2005}
}
```

And for the original VPT2 code in `SPECTRO`, you can cite

```
@incollection{spectro91,
 address = {Greenwich, Connecticut},
 author = {J. F. Gaw and A. Willets and W. H. Green and N. C. Handy},
 booktitle = {Advances in Molecular Vibrations and Collision Dynamics},
 editor = {Joel M. Bowman and Mark A. Ratner},
 pages = {170-185},
 publisher = {JAI Press, Inc.},
 title = {{SPECTRO: A} Program for the Derivation of Spectroscopic Constants From Provided Quartic Force Fields and Cubic Dipole Fields},
 year = {1991}
}
```

# TODOs
- [ ] factor out commonality in first part of `CoordType::run`
- [ ] use or delete `Fitted` trait
  - `Findiff` actually has some default methods that make it useful, but
    `Fitted` doesn't provide any default methods and isn't used as a bound
