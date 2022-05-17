use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct Summary {
    pub lx: Vec<f64>,
    pub harm: Vec<f64>,
    pub fund: Vec<f64>,
    pub corr: Vec<f64>,
    pub rots: Vec<Vec<f64>>,
    pub deltas: Vec<f64>,
    pub phis: Vec<f64>,
    pub rhead: Vec<String>,
    pub ralpha: Vec<f64>,
    pub requil: Vec<f64>,
    pub fermi: Vec<String>,
    pub zpt: f64,
    pub lin: bool,
    pub imag: bool,
}

impl Summary {
    pub fn new(filename: &str) -> Self {
        let output = std::process::Command::new("/home/brent/go/bin/summarize")
            .arg("-json")
            .arg(filename)
            .output()
            .unwrap();
        let data = String::from_utf8(output.stdout).unwrap().to_lowercase();
        serde_json::from_str(&data).unwrap()
    }
}
