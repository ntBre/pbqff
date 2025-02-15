use std::{error::Error, path::Path};

use clap::Parser;
use spectro::{Spectro, SpectroFinish};

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Optional name to operate on
    #[arg(value_parser)]
    infile: String,

    /// Writes the output in JSON format for use by other programs
    #[arg(short, long, value_parser, default_value_t = false)]
    json: bool,

    /// finish a run by loading a [SpectroFinish] from a JSON file
    #[arg(short, long, value_parser, default_value_t = false)]
    finish: bool,

    /// convert a [SpectroFinish] to the format used in Gaussian output files
    #[arg(short, long, value_parser, default_value_t = false)]
    to_gauss: bool,

    /// print more information. currently this means the full listing of
    /// vibrational states
    #[arg(
        short,
        long,
        value_parser,
        default_value_t = false,
        conflicts_with = "json"
    )]
    verbose: bool,

    /// Dump normal coordinate FCs to stdout
    #[arg(short, long, value_parser, default_value_t = false)]
    dump_fcs: bool,
}

fn main() -> Result<(), Box<dyn Error>> {
    let cfg = Args::parse();
    if cfg.to_gauss {
        let sf = SpectroFinish::load(&cfg.infile)?;
        sf.to_gauss();
        return Ok(());
    }
    let (got, spectro) = if cfg.finish {
        let SpectroFinish {
            mut spectro,
            freq,
            f3qcm,
            f4qcm,
            irreps,
            lxm,
            lx,
        } = SpectroFinish::load(&cfg.infile)?;
        spectro.verbose = cfg.verbose;
        let g = spectro.finish(freq, f3qcm, f4qcm, irreps, lxm, lx);
        (g, spectro)
    } else {
        let mut spectro = Spectro::load(&cfg.infile);
        spectro.verbose = cfg.verbose;
        spectro.dump_fcs = cfg.dump_fcs;
        let infile = Path::new(&cfg.infile);
        let dir = infile.parent().unwrap_or_else(|| Path::new("."));
        let g = spectro.run_files(
            dir.join("fort.15"),
            dir.join("fort.30"),
            dir.join("fort.40"),
        );
        (g, spectro)
    };
    if cfg.json {
        let data = serde_json::to_string_pretty(&got)?;
        println!("{data}");
    } else {
        spectro.write_output(&mut std::io::stdout(), &got)?;
    }
    Ok(())
}
