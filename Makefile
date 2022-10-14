TESTFLAGS = -- --nocapture --test-threads=1

SHORT = 0

ifeq (${SHORT}, 0)
TESTFLAGS += --include-ignored
endif

test:
	RUST_BACKTRACE=1 cargo test ${TESTFLAGS} ${ARGS}

clean:
	cargo clean

cover:
	cargo tarpaulin --color=never --skip-clean ${TESTFLAGS} ${ARGS}

build:
    # see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
	RUSTFLAGS='-C target-feature=+crt-static' \
	cargo build --bin rust-pbqff --release --target x86_64-unknown-linux-gnu

BASE = /home/brent/Projects/rust-pbqff
TARGET = target/x86_64-unknown-linux-gnu/release/rust-pbqff
ELAND_DEST = 'eland:programs/rust-pbqff/.'
WOODS_DEST = 'woods:Programs/rpbqff/rpbqff'

eland: build
	scp -C ${BASE}/${TARGET} ${ELAND_DEST}

woods: build
	scp -C ${BASE}/${TARGET} ${WOODS_DEST}

profile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out	\
		--collect-jumps=yes --simulate-cache=yes		\
		${BASE}/target/release/$(1)

profile.cart:
	$(call profile,cart)

profile.build_points:
	$(call profile,build_points)
