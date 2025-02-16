# INTDER
intder performs simple-/symmetry-internal coordinate transformations for geometry
displacements and derivatives. This is a port of the original Fortran INTDER program
by Wesley D. Allen and coworkers.

# Input
The input format is compatible with the original Fortran version,
but whitespace is much less significant. Fields just need to be separated by
some number of tabs or spaces.

## Displacements
A simple input file for converting 4 SIC displacements for water looks like:
```text
# INTDER ##########################
    3    3    3    0    0    3    0    0    0    1    0    0    0    1    1    0
   14
STRE     1    2
STRE     2    3
BEND     1    2    3
    1   1   1.000000000   2   1.000000000
    2   3   1.000000000
    3   1   1.000000000   2  -1.000000000
    0
      0.000000000000     1.431390244079     0.986041163966
      0.000000000000     0.000000000000    -0.124238450265
      0.000000000000    -1.431390244079     0.986041163966
DISP   4
    1       0.0050000000
    0
    2       0.0050000000
    0
    3       0.0050000000
    0
    1      -0.0050000000
    2      -0.0050000000
    3      -0.0100000000
    0
```

## Force Constants
An input file for converting force constants from SICs to Cartesian coordinates,
again for water, looks like:
```text
# INTDER ##########################
    3    3    3    4    0    3    2    0    0    1    3    0    0    0    0    0
   14
STRE     1    2
STRE     2    3
BEND     1    2    3
    1   1   1.000000000   2   1.000000000
    2   3   1.000000000
    3   1   1.000000000   2  -1.000000000
    0
        0.0000000000        1.4313273344        0.9860352735
        0.0000000000        0.0000000000       -0.1242266321
        0.0000000000       -1.4313273344        0.9860352735
         H1          O16          H1
    1    1    0    0      8.360863692425
    2    1    0    0      0.364250381719
    2    2    0    0      0.705590041836
    3    3    0    0      8.562725561924
    0
    1    1    1    0    -41.638868371821
    2    1    1    0     -0.611029974345
    2    2    1    0     -0.447356783193
    2    2    2    0     -0.701565547377
    3    3    1    0    -41.484904978222
    3    3    2    0      0.392943158463
    0
    1    1    1    1    181.917347386091
    2    1    1    1     -0.292134838642
    2    2    1    1      0.372522604963
    2    2    2    1      1.034069182728
    2    2    2    2     -0.655831802100
    3    3    1    1    182.206178030730
    3    3    2    1     -1.233550191665
    3    3    2    2     -0.820459303061
    3    3    3    3    183.621273961541
    0
```

# References
  - D. F. McIntosh, K. H. Michaelian, and M. R. Peterson. Can. J. Chem. Vol. 56,
    1978
  - E. B. Wilson, Jr., J. C. Decius, and P. C. Cross. Molecular
    Vibrations, 1955.
  - W. D. Allen and A. G. Csaszar. J. Chem. Phys. 98, 1993.
  - A. L. L. East, W. D. Allen, and S. J. Klippenstein. J. Chem. Phys. 102, 1995

