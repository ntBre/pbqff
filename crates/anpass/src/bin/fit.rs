use anpass::Anpass;

fn main() {
    let anpass = Anpass::load_file("testfiles/c3h2.in");
    anpass.fit();
}
