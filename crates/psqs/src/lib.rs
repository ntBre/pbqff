pub mod geom;
pub mod program;
pub mod queue;

#[cfg(test)]
mod tests;

/// from [StackOverflow](https://stackoverflow.com/a/45145246)
#[macro_export]
macro_rules! string {
    // match a list of expressions separated by comma:
    ($($str:expr),*) => ({
        // create a Vec with this list of expressions,
        // calling String::from on each:
        vec![$(String::from($str),)*] as Vec<String>
    });
}

/// call `rayon::ThreadPoolBuilder` to set `num_threads` to `n`. Discards the
/// error returned by `build_global` if the thread pool has already been
/// initialized
pub fn max_threads(n: usize) {
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global();
}
