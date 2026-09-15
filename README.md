# FiberStack

Mature tools for single-molecule Fiber-seq analysis, one repository each, pinned together here
as git submodules with a single `uv` workspace so that one environment runs all of them and the
glue between them has a home.

| tool | what it does | depends on |
| --- | --- | --- |
| [`fiberstore/`](fiberstore) | BAM → compact per-fiber store (Parquet + in-memory index); millisecond region queries returning m6A, 5mC, nucleosome and MSP coordinates per molecule | numpy, pyarrow (+ pysam to build) |
| [`hier-hmm/`](hier-hmm) | two-level HMM segmenting each molecule into open / linker / nucleosome, and footprints inside open runs | numpy, pyyaml, matplotlib |
| `fiberstack/` | code that needs two tools at once — currently `to_hmm`: a fiberstore window → hier-hmm input array / `.npz` | the above |

Data flow:

```
BAM ──fiberstore build──▶ store/ ──FiberStore.query──▶ per-fiber arrays
                                        │
                                        └─fiberstack.to_hmm──▶ (n_fiber, n_pos) NaN/0/1 ──hier-hmm──▶ chromatin states
```

## Setup

```bash
git clone --recurse-submodules <this repo>
cd FiberStack
uv sync --extra build --extra test      # one .venv with every tool installed editable
uv run pytest fiberstore hier-hmm tests
```

An analysis project installs what it needs from here (`uv pip install -e ../FiberStack/fiberstore`)
and imports it; it never copies the code.

## Using the glue

```python
from fiberstore import FiberStore
from fiberstack import run_window, save_window

fs = FiberStore("treg_store")
seq = fasta.fetch("chr1", 100_000_000, 100_001_000)          # any way of getting the sequence
res, d = run_window(fs, "chr1", 100_000_000, 100_001_000, seq, min_frac=0.88)
res.labels                                                   # (n_fiber, n_bp), see hier-hmm
save_window("win.npz", fs, "chr1", 100_000_000, 100_001_000, seq)   # then: hier-hmm run win.npz
```

Positions whose reference base is not A/T are non-callable (NaN); positions outside a fiber's
aligned span are NaN with `coverage = False`.

## What goes in here

A tool graduates into FiberStack when it has its own repository with a `pyproject.toml`, tests
that run without any project's data, a README stating its input/output contract, and no paths
or parameters specific to one analysis. Exploratory code stays in the analysis project that
produced it.

## Adding a tool

```bash
git submodule add https://github.com/mchikina/<tool>.git <tool>
# add "<tool>" to [tool.uv.workspace].members and [tool.uv.sources] in pyproject.toml
uv lock
```

## License

MIT.
