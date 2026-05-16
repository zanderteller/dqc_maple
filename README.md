# DQC for Maple

C++ implementation of **Dynamic Quantum Clustering** (DQC), packaged as a shared library that Maple can call through a worksheet module. DQC originated at Quantum Insights (QI), where it was known as "Dynamic Quantum Clustering"; the same algorithmic core lives at the heart of the open-source **DQM** package.

> **Looking for the Python interface?** See [zanderteller/dqm](https://github.com/zanderteller/dqm) — the current open-source DQM release.

This repository contains:

- `maple/CoreDQC.mw` — Maple worksheet that defines the `CoreDQC` module and binds it to the shared library
- `cpp/` — C++ source for the shared library, with a Makefile-driven build

---

## Quick start (macOS, Apple Silicon)

```bash
xcode-select --install         # one-time, if you don't already have it
brew install openblas          # installs OpenBLAS and libomp
cd cpp && make                 # builds cpp/build/dqc_maple.dylib
mkdir -p ~/maple/toolbox/CoreDQC/bin ~/maple/toolbox/CoreDQC/lib
cp build/dqc_maple.dylib ~/maple/toolbox/CoreDQC/bin/
cp ../maple/CoreDQC.mw ~/maple/toolbox/CoreDQC/lib/
```

Then open `~/maple/toolbox/CoreDQC/lib/CoreDQC.mw` in Maple and evaluate the worksheet. Inside that session:

```maple
with(CoreDQC):
SecondaryCountC();        # should return 0
```

For details on each step, and for other platforms, read on.

---

## Build

### Requirements by platform

| Platform                | Status | Toolchain | Dependencies |
|---|---|---|---|
| **macOS, Apple Silicon** | verified | Xcode CLT (`clang`) | Maple 2022+, Homebrew, `openblas` (pulls in `libomp`) |
| **macOS, Intel**         | should work, not verified | Xcode CLT | Maple 2022+, Homebrew, `openblas` |
| **Linux**                | should work, not verified | `gcc`, `make` | Maple, `libopenblas-dev`, `liblapacke-dev`, `libomp-dev` |
| **Windows**              | not covered by the Makefile | — | The Visual Studio `.sln`/`.vcxproj` files in `cpp/dqc_core/` and `cpp/dqc_maple/` are a possible starting point but aren't currently maintained. |

### macOS, Apple Silicon

```bash
xcode-select --install         # one-time
brew install openblas          # also installs libomp transitively
cd cpp
make                           # or `make all`
```

Produces:
- `cpp/build/dqc_maple.dylib` — the library Maple loads
- `cpp/build/secondary`, `testmatrix`, `testharness`, `testharnessdist` — standalone test executables (optional, useful for debugging)

### macOS, Intel

Identical steps to Apple Silicon — the Makefile auto-detects the CPU architecture and picks the right defaults (Intel Homebrew prefix `/usr/local`, and Maple's `bin.APPLE_UNIVERSAL_OSX` x86_64 binaries):

```bash
xcode-select --install
brew install openblas
cd cpp && make
```

If you have a non-Homebrew install of OpenBLAS or libomp, you can override on the command line: `make OPENBLASROOT=... LIBOMPROOT=...`.

### Linux

```bash
sudo apt install build-essential libopenblas-dev liblapacke-dev libomp-dev
cd cpp
make MAPLEROOT=/path/to/your/maple
```

The Makefile defaults to `MAPLEROOT=/usr/local/maple/current`. Set `MAPLEROOT=` to the directory that contains `bin.X86_64_LINUX/libmaplec.so` and `extern/include/maplec.h`. (Equivalent packages exist on Fedora/RHEL: `openblas-devel`, `lapack-devel`, `libomp-devel` or similar — names vary slightly by distribution.)

### Build outputs

The only built artifact you need for the Maple install (alongside `maple/CoreDQC.mw` from the repo) is the shared library:

- macOS: `cpp/build/dqc_maple.dylib`
- Linux: `cpp/build/dqc_maple.so`
- (Windows: `dqc_maple.dll`, if you build via Visual Studio)

The other binaries in `cpp/build/` are standalone test/diagnostic tools; you can ignore them.

---

## Install in Maple

The Maple worksheet hardcodes the dylib's expected location to a standard Maple toolbox path:

```
~/maple/toolbox/CoreDQC/
├── bin/dqc_maple.dylib    (or .so / .dll, depending on platform)
└── lib/CoreDQC.mw          (or CoreDQC.mla after compiling — see below)
```

Copy both files into place:

```bash
mkdir -p ~/maple/toolbox/CoreDQC/bin ~/maple/toolbox/CoreDQC/lib
cp cpp/build/dqc_maple.dylib ~/maple/toolbox/CoreDQC/bin/    # adjust extension on Linux/Windows
cp maple/CoreDQC.mw ~/maple/toolbox/CoreDQC/lib/
```

### Getting `with(CoreDQC):` to work

A `.mw` worksheet is a Maple source file — Maple's `with(...)` command will *not* find a raw `.mw` automatically. You have two options:

**Option A — per-session (quickest):**  
Open `~/maple/toolbox/CoreDQC/lib/CoreDQC.mw` in Maple and evaluate the whole worksheet (Edit → Execute → Worksheet, or use the toolbar). That defines the `CoreDQC` module in the current session. `with(CoreDQC):` works for the rest of that session, but you have to redo it each time Maple starts.

**Option B — one-time install (recommended):**  
Compile the worksheet to a `.mla` archive that Maple discovers automatically. After opening `CoreDQC.mw` and evaluating it, run this in the same worksheet:

```maple
LibraryTools[Save](CoreDQC, cat(kernelopts(homedir), "/maple/toolbox/CoreDQC/lib/CoreDQC.mla"));
```

(Or hardcode your path instead of `kernelopts(homedir)`.)

From that point on, any Maple session can do `with(CoreDQC):` directly — Maple finds `CoreDQC.mla` in the standard toolbox `lib/` directory without you opening the worksheet.

You can keep the `.mw` alongside the `.mla` (for reference / future edits) or delete it once the `.mla` is in place.

### Verify

In any Maple session after installing:

```maple
with(CoreDQC):
SecondaryCountC();
```

Expected result: `0` (no remote-secondary compute nodes are registered — that's the empty default state). If `define_external` couldn't locate or load the dylib, you'll get a clear error message pointing at the path Maple was looking for.

---

## Repository layout

```
.
├── README.md
├── LICENSE
├── maple/
│   └── CoreDQC.mw        # Maple worksheet — the user-facing entry point
└── cpp/
    ├── Makefile          # top-level (delegates to subdirs)
    ├── dqc_core/         # platform-agnostic DQC C++ core + diagnostic test executables
    │   ├── Makefile
    │   └── *.cpp / *.h
    └── dqc_maple/        # Maple bindings on top of dqc_core
        ├── Makefile
        └── dqc_maple.cpp
```

---

## License

[MIT](./LICENSE).
