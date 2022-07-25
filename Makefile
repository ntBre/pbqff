TESTFLAGS = -- --nocapture --test-threads=1

SHORT = 0

ifeq (${SHORT}, 0)
TESTFLAGS += --include-ignored
endif

test:
	RUST_BACKTRACE=1 cargo test ${TESTFLAGS} ${ARGS}

BASE = /home/brent/Projects/rust-pbqff
ELAND_DEST = 'eland:programs/rust-pbqff/.'
eland:
# see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
	RUSTFLAGS='-C target-feature=+crt-static' \
	cargo build --release --target x86_64-unknown-linux-gnu
	scp -C ${BASE}/target/x86_64-unknown-linux-gnu/release/rust-pbqff ${ELAND_DEST}
