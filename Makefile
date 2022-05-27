TESTFLAGS = -- --nocapture --test-threads=1

test:
	RUST_BACKTRACE=1 cargo test ${TESTFLAGS} ${ARGS}
