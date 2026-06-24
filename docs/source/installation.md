# Installation

CRISPRware has two parts: a **Python environment**, and the **`crispr-ots`** off-target scanner (a
standalone Rust binary). Install the env once, then put `crispr-ots` on your `PATH`, either a prebuilt
binary or one you build.

## 1. Python environment

```bash
git clone https://github.com/ericmalekos/crisprware && cd crisprware
conda env create -f environment.yml        # `conda` -> `micromamba` if you prefer
conda activate crisprware
pip install .
```

## 2. The crispr-ots scanner

### Prebuilt binary

Download the build for your platform from the
[Releases page](https://github.com/ericmalekos/crisprware/releases), extract it, and put it on your
`PATH`:

```bash
# example: Linux x86-64
curl -L -o crispr-ots.tar.gz \
  https://github.com/ericmalekos/crisprware/releases/latest/download/crispr-ots-x86_64-unknown-linux-gnu.tar.gz
tar xzf crispr-ots.tar.gz
install -m 755 crispr-ots ~/.local/bin/     # any directory on your PATH
crispr-ots --version
```

Targets published per release: `x86_64-unknown-linux-gnu`, `x86_64-unknown-linux-musl` (fully static,
no glibc dependency, for old distros), `aarch64-apple-darwin` (Apple Silicon), `x86_64-apple-darwin`
(Intel Mac), `x86_64-pc-windows-msvc`.

### Build it yourself

With a [Rust toolchain](https://rustup.rs):

```bash
cd crispr-ots
cargo install --path crates/crispr-cli      # -> ~/.cargo/bin/crispr-ots
# or build in place: cargo build --release -p crispr-cli  (-> crispr-ots/target/release/crispr-ots)
```

The build defaults to a **portable x86-64 baseline** that runs on any x86-64 CPU. For a faster binary
on modern hardware, opt in:

```bash
RUSTFLAGS="-C target-cpu=x86-64-v3" cargo install --path crates/crispr-cli   # Haswell 2013+
# or -C target-cpu=native to tune for the build machine
```

```{warning}
A `x86-64-v3` binary raises SIGILL on pre-2013 CPUs (and building with v3 on such a host fails too).
Use the baseline default on older hardware, e.g. cluster login nodes.
```

## Verify

```bash
crisprware -h          # Python CLI
crispr-ots --version   # scanner on PATH
```

## Docker

```bash
docker pull ericmalekos/crisprware:latest
docker run crisprware -h
```
