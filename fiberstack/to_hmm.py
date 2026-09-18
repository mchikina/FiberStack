"""fiberstore -> hier-hmm bridge.

hier-hmm wants one (n_fiber, n_pos) array per window: NaN = not callable, 0 = callable
unmethylated, 1 = methylated; optional `coverage` (bool, same shape), `positions`,
`center`, and any length-n_fiber array as per-fiber metadata.  `window_matrix` builds
that from a FiberStore; `save_window` writes the .npz `hier-hmm run` reads.
"""
import numpy as np

CALLABLE = "AT"


def window_matrix(fs, chrom, s, e, seq=None, min_frac=0.88, center=None, columns=()):
    """Fibers covering >= min_frac of [s, e) as a hier-hmm input dict.

    seq: reference sequence of [s, e) (upper/lower case); positions whose base is not
    A/T are marked non-callable unless the fiber carries a call there (read/reference
    mismatch), in which case the call is kept.  With seq=None every covered position is callable,
    which dilutes methylation rates by ~2x — pass the sequence when you can.
    columns: extra store columns to carry along as per-fiber metadata.
    """
    n_pos = e - s
    t = fs.query(chrom, s, e, min_frac=min_frac,
                 columns=list(dict.fromkeys(["start", "end", "reverse", "movie", "zmw", "m6a", *columns])))
    n = len(t["row"])
    cov = np.zeros((n, n_pos), bool)
    m6a = np.full((n, n_pos), np.nan, np.float32)
    if seq is not None:
        if len(seq) != n_pos:
            raise ValueError(f"seq has {len(seq)} bases, window has {n_pos}")
        callable_pos = np.frombuffer(seq.upper().encode(), np.uint8)
        callable_pos = np.isin(callable_pos, np.frombuffer(CALLABLE.encode(), np.uint8))
    else:
        callable_pos = np.ones(n_pos, bool)
    for i in range(n):
        a, b = max(int(t["start"][i]) - s, 0), min(int(t["end"][i]) - s, n_pos)
        cov[i, a:b] = True
        row = m6a[i]
        row[a:b] = np.where(callable_pos[a:b], 0.0, np.nan)
        p = t["m6a"][i] - s
        p = p[(p >= 0) & (p < n_pos)]
        row[p] = 1.0
    out = {"m6a": m6a, "coverage": cov, "positions": np.arange(s, e),
           "center": np.int64((s + e) // 2 if center is None else center),
           "qname": np.array([fs.qname(m, z) for m, z in zip(t["movie"], t["zmw"])]),
           "start": t["start"], "end": t["end"], "reverse": t["reverse"], "row": t["row"],
           "n_m6a": np.array([len(p) for p in t["m6a"]])}   # whole-read m6A count
    for c in columns:
        if c not in out:
            out[c] = t[c]
    return out


def save_window(path, fs, chrom, s, e, seq=None, **kw):
    """window_matrix -> .npz for `hier-hmm run`.  Only 1-D per-fiber arrays ride along."""
    d = window_matrix(fs, chrom, s, e, seq, **kw)
    n = d["m6a"].shape[0]
    keep = {k: v for k, v in d.items()
            if k in ("m6a", "coverage", "positions", "center")
            or (isinstance(v, np.ndarray) and v.ndim == 1 and len(v) == n)}
    np.savez_compressed(path, **keep)
    return keep


def run_window(fs, chrom, s, e, seq=None, cfg=None, **kw):
    """Segment one window with hier-hmm; returns (result, input dict)."""
    from hier_hmm import load_config, prepare, run
    cfg = cfg or load_config()
    d = window_matrix(fs, chrom, s, e, seq, **kw)
    n = d["m6a"].shape[0]
    meta = {k: v for k, v in d.items() if isinstance(v, np.ndarray) and v.ndim == 1 and len(v) == n}
    ds = prepare(d["m6a"], cfg, coverage=d["coverage"], meta=meta,
                 extras={"positions": d["positions"], "center": d["center"]})
    return run(ds, cfg), d
