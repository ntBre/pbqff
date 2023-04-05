TESTFLAGS = -- --nocapture --test-threads=1

SHORT = 0

ifeq (${SHORT}, 0)
TESTFLAGS += --include-ignored
endif

test:
	RUST_BACKTRACE=1 cargo test ${TESTFLAGS} ${ARGS}

clean:
	cargo clean

clippy:
	cargo clippy --all

cover:
	cargo tarpaulin --color=never --skip-clean ${TESTFLAGS} ${ARGS}

TARGET =
ifeq (${DEBUG}, 1)
	TARGET += target/x86_64-unknown-linux-gnu/debug/rust-pbqff
else
	TARGET += target/x86_64-unknown-linux-gnu/release/rust-pbqff
endif

build:
ifeq (${DEBUG}, 1)
    # see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
	RUSTFLAGS='-C target-feature=+crt-static' \
	cargo build --features vers --bin rust-pbqff --target x86_64-unknown-linux-gnu
else
    # see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
	RUSTFLAGS='-C target-feature=+crt-static' \
	cargo build --features vers --bin rust-pbqff --release --target x86_64-unknown-linux-gnu
endif


BASE = /home/brent/Projects/rust-pbqff
ELAND_DEST = 'eland:programs/rust-pbqff/.'
WOODS_DEST = 'woods:bin/rpbqff'

eland: build
	scp -C ${BASE}/${TARGET} ${ELAND_DEST}

woods: build docs
	scp -C ${BASE}/${TARGET} ${WOODS_DEST}

woods.alpha: build
	scp -C ${BASE}/${TARGET} ${WOODS_DEST}.alpha

docs: man/rpbqff.1
	scp -C $? 'woods:man/man1/.'
	date > docs

%.pdf: %.1
	groff -Tpdf $? -mman > $@

man/rpbqff.1: man/rpbqff.head testfiles/test.toml man/rpbqff.tail
	cat $^ > $@

scripts: qffbuddy/qffbuddy*
	scp -C $? 'woods:bin/'
	date > scripts

profile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out	\
		--collect-jumps=yes --simulate-cache=yes		\
		${BASE}/target/release/$(1)

profile.cart:
	$(call profile,cart)

profile.build_points:
	$(call profile,build_points)
