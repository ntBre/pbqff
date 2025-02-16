# psqs

quantum chemistry ProgramS and Queuing systemS

# Description

This is a crate for Rust interfaces to quantum chemistry programs (like Molpro,
Gaussian, and MOPAC) and queuing systems like Slurm and PBS. It contains traits
defining the high-level behavior of these `Program`s and `Queue`s, as well as
concrete implementations for the aforementioned programs. The list of currently
supported programs is as follows:

## Programs

The following quantum chemistry programs are supported:

- [x] [MOPAC][mopac]
- [x] [Molpro][molpro]
- [x] [CFOUR][cfour]
- [x] [DFTB+][dftb+]
- [ ] Gaussian

## Queuing systems

The following queuing systems are supported:

- [x] Slurm
- [x] PBS

[mopac]: http://openmopac.net/
[molpro]: https://www.molpro.net/
[cfour]: https://cfour.uni-mainz.de/cfour/
[dftb+]: https://dftbplus.org/
