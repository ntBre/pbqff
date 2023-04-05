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
    let queues = ["Pbs", "Slurm", "Local"];
    let mut s = String::from(
        "fn dispatch(config: &Config, args: Args) -> (Spectro, Output) {
    let m = (config.coord_type, config.program, config.queue, args.checkpoint);
match m {
",
    );
    for (i, coord) in coord_types.iter().enumerate() {
        for program in programs {
            for queue in queues {
                for (val, fun) in [(true, "resume"), (false, "run")] {
                    let mut resume = String::from("config,");
                    if val {
                        write!(
                            resume,
                            "<{coord} as CoordType<::std::io::Stdout, \
			     {queue}, {program}>>::Resume::load(\"res.chk\")"
                        )
                        .unwrap();
                    }
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
&queue!({queue}, config, args.no_del),
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
