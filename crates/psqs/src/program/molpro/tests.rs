use std::{fs::read_to_string, path::Path, str::FromStr};

use crate::{
    geom::Geom,
    program::{molpro::Molpro, Procedure, Program, Template},
};

fn opt_templ() -> Template {
    Template::from(
        "
memory,1,g

gthresh,energy=1.d-12,zero=1.d-22,oneint=1.d-22,twoint=1.d-22;
gthresh,optgrad=1.d-8,optstep=1.d-8;
nocompress;

geometry={
{{.geom}}
basis={
default,cc-pVTZ-f12
}
set,charge={{.charge}}
set,spin=0
hf,accuracy=16,energy=1.0d-10
{CCSD(T)-F12,thrden=1.0d-8,thrvar=1.0d-10}
{optg,grms=1.d-8,srms=1.d-8}
",
    )
}

fn single_templ() -> Template {
    Template::from(
        "
memory,1,g

gthresh,energy=1.d-12,zero=1.d-22,oneint=1.d-22,twoint=1.d-22;
gthresh,optgrad=1.d-8,optstep=1.d-8;
nocompress;

geometry={
{{.geom}}
basis={
default,cc-pVTZ-f12
}
set,charge={{.charge}}
set,spin=0
hf,accuracy=16,energy=1.0d-10
{CCSD(T)-F12,thrden=1.0d-8,thrvar=1.0d-10}
",
    )
}

enum Type {
    Opt,
    Spt,
}

fn test_molpro(t: Type, path: impl AsRef<Path>) -> Molpro {
    let path = path.as_ref();
    Molpro::new(
        path.to_string_lossy().to_string(),
        match t {
            Type::Opt => opt_templ(),
            Type::Spt => single_templ(),
        },
        0,
        Geom::from_str(
            "C
C 1 CC
C 1 CC 2 CCC
H 2 CH 1 HCC 3 180.0
H 3 CH 1 HCC 2 180.0

CC =                  1.42101898
CCC =                55.60133141
CH =                  1.07692776
HCC =               147.81488230
",
        )
        .unwrap(),
    )
}

/// in these names, the first word is the template type (opt => optg line
/// included in template, for example), and the second word is the Procedure
pub(crate) mod write_input {
    use super::*;

    use tempfile::NamedTempFile;
    use test_case::test_case;

    #[macro_export]
    macro_rules! check {
        ($want_file: expr) => {
            check!($want_file, "/tmp/opt.inp");
        };
        ($want_file: expr, $got_file: expr) => {
            let want_file = $want_file;
            let got = read_to_string($got_file).expect("file not found");
            let want = read_to_string(want_file).unwrap();
            if got != want {
                panic!("\ngot:\n{}\nwant:\n{}", got, want);
            }
        };
    }

    #[test_case(Type::Opt, Procedure::Opt, "testfiles/molpro/opt.want")]
    #[test_case(Type::Opt, Procedure::SinglePt, "testfiles/molpro/single.want")]
    #[test_case(Type::Spt, Procedure::Opt, "testfiles/molpro/opt.want")]
    #[test_case(Type::Spt, Procedure::SinglePt, "testfiles/molpro/single.want")]
    fn opt_opt(ty: Type, proc: Procedure, want: &str) {
        let dir = NamedTempFile::new().unwrap();
        let mut m = test_molpro(ty, &dir);
        m.write_input(proc);
        check!(want, &*dir.path().with_extension("inp").to_string_lossy());
    }
}

mod read_output {
    use crate::program::ProgramResult;
    use symm::Atom;

    use super::*;

    #[test]
    fn opt() {
        let got = Molpro::read_output("testfiles/molpro/opt").unwrap();
        let want = ProgramResult {
            energy: -76.369839620286,
            cart_geom: Some(vec![
                //
                Atom::new_from_label(
                    "O",
                    0.0000000000,
                    0.0000000000,
                    -0.0657441581,
                ),
                Atom::new_from_label(
                    "H",
                    0.0000000000,
                    0.7574590773,
                    0.5217905246,
                ),
                Atom::new_from_label(
                    "H",
                    0.0000000000,
                    -0.7574590773,
                    0.5217905246,
                ),
            ]),
            time: 27.13,
        };

        assert_eq!(got, want);
    }

    #[test]
    fn dzccr() {
        let got = Molpro::read_output("testfiles/molpro/dzccr");
        let got = got.unwrap_or_else(|e| panic!("{e:#?}"));
        let want = ProgramResult {
            energy: -76.470698498340,
            cart_geom: None,
            time: 4.73,
        };

        assert_eq!(got, want);
    }

    #[test]
    fn error() {
        let got = Molpro::read_output("testfiles/molpro/error");
        let Err(e) = got else {
            panic!("expected error got {got:?}");
        };
        assert!(e.is_error_in_output());
    }

    #[test]
    fn ignore_error() {
        let got = Molpro::read_output("testfiles/molpro/ignore_error");
        assert!(got.is_ok());
    }
}
