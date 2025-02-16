use std::fmt::Display;
use std::time::Duration;

#[derive(Default)]
pub(crate) struct Timer {
    pub(crate) writing_input: Duration,
    pub(crate) writing_script: Duration,
    pub(crate) submitting_script: Duration,
    pub(crate) reading: Duration,
    pub(crate) sleeping: Duration,
    pub(crate) removing: Duration,
}

impl Display for Timer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:.1} s reading ok, {:.1} s writing input, {:.1} s writing script, \
	     {:.1} s submitting, {:.1} s sleeping, {:.1} s removing",
            self.reading.as_millis() as f64 / 1000.0,
            self.writing_input.as_millis() as f64 / 1000.0,
            self.writing_script.as_millis() as f64 / 1000.0,
            self.submitting_script.as_millis() as f64 / 1000.0,
            self.sleeping.as_millis() as f64 / 1000.0,
            self.removing.as_millis() as f64 / 1000.0,
        )
    }
}
