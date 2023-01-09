//! build the big match on config parameters that changes how `run` or `resume`
//! get called.

extern crate proc_macro;
use proc_macro::TokenStream;
use std::fmt::Write;

#[proc_macro]
pub fn dispatch(_item: TokenStream) -> TokenStream {
    let coord_types = ["Normal", "Cart", "Sic"];
    let coord_builders = [
        "Normal::findiff(config.findiff)",
        "Cart",
        "Sic::new(Intder::load_file(\"intder.in\"))",
    ];
    let programs = ["Mopac", "Molpro"];
    let queues = ["Pbs", "Slurm"];
    let mut s = String::from(
        "fn dispatch(config: &Config, checkpoint: bool) -> (Spectro, Output) {
    let m = (config.coord_type, config.program, config.queue, checkpoint);
match m {
",
    );
    for (i, coord) in coord_types.iter().enumerate() {
        for program in programs {
            for queue in queues {
                for (val, fun) in [(true, "resume"), (false, "run")] {
                    let resume = if val {
                        format!(
                            "<{coord} as CoordType<::std::io::Stdout, \
			     {queue}, {program}>>::Resume::load(\"res.chk\")"
                        )
                    } else {
                        "config".to_owned()
                    };
                    write!(
                        s,
                        "(
config::CoordType::{coord},
config::Program::{program},
config::Queue::{queue},
{val},
) => <{coord} as CoordType<_, _, {program}>>::{fun}(
{},
&mut std::io::stdout(),
&queue!({queue}, config),
{resume}
),",
                        coord_builders[i],
                    )
                    .unwrap();
                }
            }
        }
    }
    write!(s, "_ => panic!(\"unsupported combination\"),}}}}").unwrap();
    s.parse().unwrap()
}
