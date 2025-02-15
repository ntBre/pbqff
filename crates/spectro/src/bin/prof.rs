use spectro::Spectro;

fn main() -> Result<(), std::io::Error> {
    let spectro = Spectro::load("spectro/testfiles/c3h2/spectro.in");
    let got = spectro.run_files(
        "spectro/testfiles/c3h2/fort.15",
        "spectro/testfiles/c3h2/fort.30",
        "spectro/testfiles/c3h2/fort.40",
    );
    spectro.write_output(&mut std::io::stdout(), &got)?;
    Ok(())
}
