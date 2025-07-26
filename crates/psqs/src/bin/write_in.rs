use psqs::{
    geom::Geom,
    program::{
        Procedure, Program, Template,
        mopac::{Mopac, Params},
    },
};

fn test_mopac() -> Mopac {
    let names = vec![
        "USS", "ZS", "BETAS", "GSS", "USS", "UPP", "ZS", "ZP", "BETAS",
        "BETAP", "GSS", "GPP", "GSP", "GP2", "HSP",
    ];
    let atoms = vec![
        "H", "H", "H", "H", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
        "C",
    ];
    #[rustfmt::skip]
        let values = vec![
            -11.246958000000, 1.268641000000, -8.352984000000,
            14.448686000000, -51.089653000000, -39.937920000000,
            2.047558000000, 1.702841000000, -15.385236000000,
            -7.471929000000, 13.335519000000, 10.778326000000,
            11.528134000000, 9.486212000000, 0.717322000000,
        ];
    Mopac::new_full(
        String::from("/tmp/test"),
        Some(Params::from(
            names.iter().map(|s| s.to_string()).collect(),
            atoms.iter().map(|s| s.to_string()).collect(),
            values,
        )),
        Geom::Xyz(Vec::new()),
        0,
        Template::from("scfcrt=1.D-21 aux(precision=14) PM6 A0"),
    )
}

fn main() {
    let mut tm = test_mopac();
    tm.param_dir = Some("/tmp".to_string());
    let mut res = Vec::new();
    for _ in 0..1000 {
        tm.write_input(Procedure::SinglePt);
        res.push(());
    }
}
