use std::str::FromStr;

use insta::{assert_debug_snapshot, assert_snapshot};

use super::*;

#[test]
fn write_input() {
    let template = Template::from(
        "
Geometry = xyzFormat {
{{.geom}}
}

Hamiltonian = DFTB {
  Scc = Yes
  SlaterKosterFiles = Type2FileNames {
    Prefix = \"/opt/dftb+/slako/mio/mio-1-1/\"
    Separator = \"-\"
    Suffix = \".skf\"
  }
  MaxAngularMomentum {
    O = \"p\"
    H = \"s\"
  }
  Charge = {{.charge}}
}

Options {
}

Analysis {
  CalculateForces = Yes
}

ParserOptions {
  ParserVersion = 12
}
",
    );

    let mut d = DFTBPlus {
        filename: "/tmp".into(),
        template,
        charge: 0,
        geom: Geom::from_str(
            "    3
Geometry Step: 9
    O      0.00000000     -0.71603315      0.00000000
    H      0.00000000     -0.14200298      0.77844804
    H     -0.00000000     -0.14200298     -0.77844804
",
        )
        .unwrap(),
    };

    d.write_input(Procedure::Opt);

    assert_snapshot!(read_to_string("/tmp/dftb_in.hsd").unwrap(), @r#"
    Geometry = xyzFormat {
    3

    O       0.0000000000  -0.7160331500   0.0000000000
    H       0.0000000000  -0.1420029800   0.7784480400
    H      -0.0000000000  -0.1420029800  -0.7784480400


    }

    Hamiltonian = DFTB {
      Scc = Yes
      SlaterKosterFiles = Type2FileNames {
        Prefix = "/opt/dftb+/slako/mio/mio-1-1/"
        Separator = "-"
        Suffix = ".skf"
      }
      MaxAngularMomentum {
        O = "p"
        H = "s"
      }
      Charge = 0
    }

    Options {
    }

    Analysis {
      CalculateForces = Yes
    }

    ParserOptions {
      ParserVersion = 12
    }
    Driver = GeometryOptimization {
      Optimizer = Rational {}
      MovedAtoms = 1:-1
      MaxSteps = 100
      OutputPrefix = "geom.out"
      Convergence {
            Energy = 1e-8
            GradElem = 1e-7
            GradNorm = 1e-7
            DispElem = 1e-7
            DispNorm = 1e-7
    }
    }
    "#);

    d.write_input(Procedure::SinglePt);

    assert_snapshot!(read_to_string("/tmp/dftb_in.hsd").unwrap(), @r#"
    Geometry = xyzFormat {
    3

    O       0.0000000000  -0.7160331500   0.0000000000
    H       0.0000000000  -0.1420029800   0.7784480400
    H      -0.0000000000  -0.1420029800  -0.7784480400


    }

    Hamiltonian = DFTB {
      Scc = Yes
      SlaterKosterFiles = Type2FileNames {
        Prefix = "/opt/dftb+/slako/mio/mio-1-1/"
        Separator = "-"
        Suffix = ".skf"
      }
      MaxAngularMomentum {
        O = "p"
        H = "s"
      }
      Charge = 0
    }

    Options {
    }

    Analysis {
      CalculateForces = Yes
    }

    ParserOptions {
      ParserVersion = 12
    }
    "#);

    // test that we can also handle a provided geometry driver
    let template = Template::from(
        "
Geometry = xyzFormat {
{{.geom}}
}

Hamiltonian = DFTB {
  Scc = Yes
  SlaterKosterFiles = Type2FileNames {
    Prefix = \"/opt/dftb+/slako/mio/mio-1-1/\"
    Separator = \"-\"
    Suffix = \".skf\"
  }
  MaxAngularMomentum {
    O = \"p\"
    H = \"s\"
  }
  Charge = {{.charge}}
}

Options {
}

Analysis {
  CalculateForces = Yes
}

ParserOptions {
  ParserVersion = 12
}
Driver = GeometryOptimization {
  Optimizer = Rational {}
  MovedAtoms = 1:-1
  MaxSteps = 100
  OutputPrefix = \"geom.out\"
  Convergence {
        Energy = 1e-8
        GradElem = 1e-8
        GradNorm = 1e-7
        DispElem = 1e-7
        DispNorm = 1e-7
}
}
",
    );

    let mut d = DFTBPlus {
        filename: "/tmp".into(),
        template,
        charge: 0,
        geom: Geom::from_str(
            "    3
Geometry Step: 9
    O      0.00000000     -0.71603315      0.00000000
    H      0.00000000     -0.14200298      0.77844804
    H     -0.00000000     -0.14200298     -0.77844804
",
        )
        .unwrap(),
    };

    d.write_input(Procedure::Opt);
    assert_snapshot!(read_to_string("/tmp/dftb_in.hsd").unwrap(), @r#"
    Geometry = xyzFormat {
    3

    O       0.0000000000  -0.7160331500   0.0000000000
    H       0.0000000000  -0.1420029800   0.7784480400
    H      -0.0000000000  -0.1420029800  -0.7784480400


    }

    Hamiltonian = DFTB {
      Scc = Yes
      SlaterKosterFiles = Type2FileNames {
        Prefix = "/opt/dftb+/slako/mio/mio-1-1/"
        Separator = "-"
        Suffix = ".skf"
      }
      MaxAngularMomentum {
        O = "p"
        H = "s"
      }
      Charge = 0
    }

    Options {
    }

    Analysis {
      CalculateForces = Yes
    }

    ParserOptions {
      ParserVersion = 12
    }
    Driver = GeometryOptimization {
      Optimizer = Rational {}
      MovedAtoms = 1:-1
      MaxSteps = 100
      OutputPrefix = "geom.out"
      Convergence {
            Energy = 1e-8
            GradElem = 1e-8
            GradNorm = 1e-7
            DispElem = 1e-7
            DispNorm = 1e-7
    }
    }
    "#);

    d.write_input(Procedure::SinglePt);
    assert_snapshot!(read_to_string("/tmp/dftb_in.hsd").unwrap(), @r#"
    Geometry = xyzFormat {
    3

    O       0.0000000000  -0.7160331500   0.0000000000
    H       0.0000000000  -0.1420029800   0.7784480400
    H      -0.0000000000  -0.1420029800  -0.7784480400


    }

    Hamiltonian = DFTB {
      Scc = Yes
      SlaterKosterFiles = Type2FileNames {
        Prefix = "/opt/dftb+/slako/mio/mio-1-1/"
        Separator = "-"
        Suffix = ".skf"
      }
      MaxAngularMomentum {
        O = "p"
        H = "s"
      }
      Charge = 0
    }

    Options {
    }

    Analysis {
      CalculateForces = Yes
    }

    ParserOptions {
      ParserVersion = 12
    }
    "#);
}

#[test]
fn read_opt_output() {
    let got = DFTBPlus::read_output("testfiles/dftb+/opt").unwrap();
    assert_debug_snapshot!(got, @r"
    ProgramResult {
        energy: -4.0779379326,
        cart_geom: Some(
            [
                Atom {
                    atomic_number: 8,
                    x: 0.0,
                    y: -0.71603315,
                    z: 0.0,
                    weight: None,
                },
                Atom {
                    atomic_number: 1,
                    x: 0.0,
                    y: -0.14200298,
                    z: 0.77844804,
                    weight: None,
                },
                Atom {
                    atomic_number: 1,
                    x: -0.0,
                    y: -0.14200298,
                    z: -0.77844804,
                    weight: None,
                },
            ],
        ),
        time: 0.05,
    }
    ");
}

#[test]
fn read_single_output() {
    let got = DFTBPlus::read_output("testfiles/dftb+/single").unwrap();
    assert_debug_snapshot!(got, @r"
    ProgramResult {
        energy: -3.9798793068,
        cart_geom: None,
        time: 0.03,
    }
    ");
}
