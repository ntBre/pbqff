//! Bundled submission-script templates used by PBQFF.

pub const PBS_MOPAC: &str = r#"#!/bin/sh
#PBS -N {{.basename}}
#PBS -S /bin/bash
#PBS -j oe
#PBS -o {{.filename}}.out
#PBS -W umask=022
#PBS -l walltime=1000:00:00
#PBS -l ncpus=1
#PBS -l mem=1gb
#PBS -q workq

module load openpbs

export WORKDIR=$PBS_O_WORKDIR
cd $WORKDIR

export LD_LIBRARY_PATH=/ddnlus/r2518/Packages/mopac/build
export MOPAC_CMD=/ddnlus/r2518/Packages/mopac/build/mopac
"#;

pub const PBS_MOLPRO: &str = r#"#!/bin/sh
#PBS -N {{.basename}}
#PBS -S /bin/bash
#PBS -j oe
#PBS -o {{.basename}}.out
#PBS -W umask=022
#PBS -l walltime=1000:00:00
#PBS -l ncpus=1
#PBS -l mem=8gb
#PBS -q workq

module load openpbs molpro

export WORKDIR=$PBS_O_WORKDIR
export TMPDIR=/tmp/$USER/$PBS_JOBID
cd $WORKDIR
mkdir -p $TMPDIR
trap 'rm -rf $TMPDIR' EXIT

export MOLPRO_CMD="molpro -t $NCPUS --no-xml-output"
"#;

pub const PBS_CFOUR: &str = r#"#!/bin/sh
#PBS -N {{.basename}}
#PBS -S /bin/bash
#PBS -j oe
#PBS -o {{.filename}}.out
#PBS -W umask=022
#PBS -l walltime=1000:00:00
#PBS -l ncpus=1
#PBS -l mem=8gb
#PBS -q workq

module load openpbs

export WORKDIR=$PBS_O_WORKDIR
cd $WORKDIR

CFOUR_CMD="/ddnlus/r2518/bin/c4ext_new.sh $NCPUS"
"#;

pub const PBS_DFTBPLUS: &str = r#"#!/bin/sh
#PBS -N {{.basename}}
#PBS -S /bin/bash
#PBS -j oe
#PBS -o {{.filename}}.out
#PBS -W umask=022
#PBS -l walltime=1000:00:00
#PBS -l ncpus=1
#PBS -l mem=8gb
#PBS -q workq

module load openpbs

export WORKDIR=$PBS_O_WORKDIR
cd $WORKDIR

export DFTB_CMD=/ddnlus/r2518/.conda/envs/dftb/bin/dftb+
"#;

pub const PBS_ORCA: &str = r#"#!/bin/bash
#PBS -N {{.basename}}
#PBS -S /bin/bash
#PBS -j oe
#PBS -o {{.filename}}.out
#PBS -l mem=16gb
#PBS -l nodes=1:ppn=8
#PBS -l walltime=2150:00:00
#PBS -q workq

module load openpbs

export OMP_NUM_THREADS=8
export ORCA_DIR=/ddnlus/r3750/orca/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export XTBEXE=$ORCA_DIR/xtb/xtb-dist/bin/xtb
export PATH=$ORCA_DIR:$ORCA_DIR/xtb/xtb-dist/bin:$PATH
export WORKDIR=$PBS_O_WORKDIR
export TMPDIR=/tmp/$USER.$PBS_JOBID

mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd "$WORKDIR"

run_orca() {
    input=$1
    name=$(basename "${input%.inp}")
    scratch="$TMPDIR/$name"
    mkdir -p "$scratch"
    cp "$input" "$scratch/$name.inp" || return
    (
        cd "$scratch" || exit 1
        "$ORCA_DIR/orca" "$name.inp" '--bind-to core:overload-allowed --use-hwthread-cpus'
    )
    status=$?
    rm -rf "$scratch"
    return $status
}
ORCA_CMD=run_orca
"#;

pub const SLURM_MOPAC: &str = include_str!("../../templates/slurm/mopac");
pub const SLURM_MOLPRO: &str = include_str!("../../templates/slurm/molpro");
pub const SLURM_CFOUR: &str = "";
pub const SLURM_DFTBPLUS: &str = "";
pub const SLURM_ORCA: &str = "ORCA_CMD=orca\n";

pub const LOCAL_MOPAC: &str = "export MOPAC_CMD=/opt/mopac/mopac
export LD_LIBRARY_PATH=/opt/mopac/\n";
pub const LOCAL_MOLPRO: &str = "";
pub const LOCAL_CFOUR: &str = "CFOUR_CMD=/opt/cfour/cfour\n";
pub const LOCAL_DFTBPLUS: &str = "DFTB_CMD=/opt/dftb+/dftb+\n";
pub const LOCAL_ORCA: &str = "ORCA_CMD=orca\n";
