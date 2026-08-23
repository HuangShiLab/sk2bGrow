# sk2bGrow — development tasks.
.PHONY: help build release test test-rust test-py fmt lint clean smoke sim install

PY ?= python3
CARGO ?= cargo

help:
	@echo "build      debug build of both crates"
	@echo "release    optimised build"
	@echo "install    editable install of the Python layer"
	@echo "test       everything (Rust + Python)"
	@echo "smoke      full-stack run: index -> profile -> output.tsv"
	@echo "sim        reproduce the design report's section 5 simulations"
	@echo "fmt/lint   rustfmt + clippy + ruff"

build:
	$(CARGO) build --workspace

release:
	$(CARGO) build --release --workspace

install:
	$(PY) -m pip install -e python/

test: test-rust test-py

test-rust:
	$(CARGO) test --workspace

test-py:
	PYTHONPATH=python $(PY) -m pytest tests/python -q

smoke: build
	./scripts/smoke.sh

sim:
	PYTHONPATH=python $(PY) -m sk2bgrow.cli simulate a --reps 150
	PYTHONPATH=python $(PY) -m sk2bgrow.cli simulate b --reps 150

fmt:
	$(CARGO) fmt --all
	-$(PY) -m ruff format python tests

lint:
	$(CARGO) clippy --workspace --all-targets -- -D warnings
	-$(PY) -m ruff check python tests

clean:
	$(CARGO) clean
	rm -rf benches/work/* .pytest_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
