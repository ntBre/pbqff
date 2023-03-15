# pbqff
pbqff automates the construction of quartic force fields (QFFs) and the process
of generating [spectral data](https://github.com/ntBre/spectro) from them

# Coordinate Types
pbqff supports running QFFs in the following coordinate systems:
- Symmetry-internal coordinates (SICs) via [intder](https://github.com/ntBre/intder)
- Cartesian coordinates
- Normal coordinates

The normal coordinates are determined automatically by running a Cartesian harmonic force field.

# Programs and Queues
pbqff supports Molpro and Mopac for computing single-point energies and the
PBS and Slurm queuing systems via [psqs](https://github.com/ntBre/psqs)

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
An optional GUI for preparing input files and running pbqff is also included in the qffbuddy directory
