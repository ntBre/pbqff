TESTFLAGS = -- --nocapture

SHORT = 0

ifeq (${SHORT}, 0)
TESTFLAGS += --include-ignored
endif

test:
	RUST_BACKTRACE=1 cargo test ${TESTFLAGS} ${ARGS}

clean:
	cargo clean

clippy:
	cargo +nightly clippy --all-features --all-targets --workspace

cover:
	cargo tarpaulin --color=never --skip-clean ${TESTFLAGS} ${ARGS}

doc:
	cargo doc --no-deps ${ARGS}

TARGET =
ifeq (${DEBUG}, 1)
	TARGET += target/x86_64-unknown-linux-gnu/debug/pbqff
else
	TARGET += target/x86_64-unknown-linux-gnu/release/pbqff
endif

build:
ifeq (${DEBUG}, 1)
    # see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
	RUSTFLAGS='-C target-feature=+crt-static' \
	cargo build --target x86_64-unknown-linux-gnu
else
    # see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
	RUSTFLAGS='-C target-feature=+crt-static' \
	cargo build --release --target x86_64-unknown-linux-gnu
endif


install:
	cargo install --path . --bin pbqff

src := $(shell find . -name '*.rs')

target/release/pbqff: $(src)
	cargo build --release

PREFIX := /usr/bin
MANDIR := /usr/local/share/man/man1
install.full: man/rpbqff.1 target/release/pbqff qffbuddy/qffbuddy.py
	sudo ln -sf $(realpath target/release/pbqff) $(PREFIX)/pbqff
	sudo ln -sf $(realpath qffbuddy/qffbuddy.py) $(PREFIX)/qffbuddy
	sudo cp $< $(MANDIR)/pbqff.1
	touch $@

ELAND_DEST = 'eland:programs/rust-pbqff/.'

.woods.mk:
	./configure

hash := $(shell git rev-parse HEAD | cut -c 1-7)

include .woods.mk
WOODS_DEST = ${WOODS}':bin/rpbqff'
ALPHA = .$(hash)

eland: build
	scp -C ${TARGET} ${ELAND_DEST}

woods: build docs
	scp -C ${TARGET} ${WOODS_DEST}

woods.alpha: build
	scp -C ${TARGET} ${WOODS_DEST}${ALPHA}

docs: man/rpbqff.1
	scp -C $? ${WOODS}':man/man1/.'
	date > docs

%.pdf: %.1
	groff -Tpdf $? -mman > $@

man/config.man: src/config.rs scripts/config.awk
	scripts/config.awk $< > $@

man/rpbqff.1: man/rpbqff.head man/config.man man/example.man testfiles/test.toml man/rpbqff.tail
	cat $^ > $@

scripts: qffbuddy/qffbuddy*
	scp -C $? ${WOODS}':bin/'
	date > scripts

profile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out	\
		--collect-jumps=yes --simulate-cache=yes		\
		target/release/$(1)

profile.cart:
	$(call profile,cart)

profile.build_points:
	$(call profile,build_points)

memprofile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
                heaptrack -o /tmp/heaptrack.pbqff.%p.zst target/release/$(1)

memprofile.cart:
	$(call memprofile,cart)
