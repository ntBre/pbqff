# SPECTRO
spectro computes spectroscopic data for asymmetric and symmetric tops using
rotational and vibrational perturbation theory. These data include
- Anharmonic vibrational frequencies
- Vibrationally-averaged principal rotational constants
- Quartic distortion constants
- Sextic distortion constants

Additionally, spectro corrects these values for the effects of type-1 and type-2
Fermi resonances and Coriolis resonances.

# Input
The main input to the spectro binary (`spectro-bin`) is the molecular geometry:
```text
# SPECTRO ##########################################
    1    1    2    1    0    0    0    4    0    0   00    0    0    0    0
    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0
# GEOM #######################################
   3   1
 1.00      0.0000000000      1.4313901416      0.9860410955
 8.00      0.0000000000      0.0000000000     -0.1242384417
 1.00      0.0000000000     -1.4313901416      0.9860410955
```
The header isn't really needed in this version, but is still required for compatibility.
spectro also requires the harmonic, cubic, and quartic force constants to be in files named
`fort.15`, `fort.30`, and `fort.40`, respectively, in the same directory as the input file.
These should be of the form
```text
        0.0000324161        0.0000000000        0.0000000000
        0.0000000000        0.0000000000        0.0000000000
        0.0000000000        0.0000000000        0.0000000000
        0.0000000000        0.3745242435        0.2350135637
...
```
See [the original manual](https://r410berry.com/static/media/spectro_manual.8f984ff8.pdf)
for details of the ordering therein.

# References
This version of spectro is a translation of the original Fortran code by J. F.
Gaw and A. Willetts, with additional modifications by J. M. L. Martin and T. J.
Lee. In addition to that source code, the additional references were of use for
this implementation:
  - Nielsen41 - The Vibration-Rotation Energies of Polyatomic Molecules
  - Nielsen51 - The Vibration-Rotation Energies of Molecules
  - Watson68 - Simplification of the Molecular Vibration-Rotation Hamiltonian
  - Hoy72 - Anharmonic Force Constant Calculations
  - Mills76 - Vibrationally Averaged Interatomic Distances
  - Hecht61 - The Vibration-Rotation Energies of Tetrahedral XY_4 Molecules Part
    I. Theory of Spherical Top Molecules
  - Herranz61 - The Rotational Structure of the Fundamental Infrared Bands of
    Methane-Type Molecules
  - Green87 - Kinetic Anharmonic Coupling in the Trihalomethanes: A Mechanism
    for Rapid Intramolecular Redistribution of CH Stretch Vibrational Energy
  - Willetts90 - Anharmonic Corrections to Vibrational Transition Intenslties
  - Pliva90 - Anharmonic Constants for Degenerate Modes of Symmetric Top
    Molecules
  - Martin95 - The Anharmonic Force Field of Ethylene, C<sub>2</sub>H<sub>4</sub>,
    by Means of Accurate *ab initio* Calculations
  - Martin97 - Accurate *ab initio* Quartic Force Field for *trans*-HNNH and Treatment
    of Resonance Polyads
  - Lehmann88 - Beyond the x-K relations: Calculations of 1–1 and 2–2 resonance
    constants with application to HCN and DCN
