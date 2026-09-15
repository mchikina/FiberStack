import numpy as np
import pytest

pysam = pytest.importorskip("pysam")

from fiberstore import FiberStore, build           # noqa: E402
from fiberstore.testing import synthetic_bam       # noqa: E402
from fiberstack.to_hmm import save_window, window_matrix  # noqa: E402


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    d = tmp_path_factory.mktemp("s")
    truth, ref = synthetic_bam(d / "syn.bam", n_reads=300, seed=3)
    build(str(d / "syn.bam"), str(d / "fs"), workers=2, chunk=10_000, log=None)
    return FiberStore(str(d / "fs")), truth, ref


def test_window_matrix_matches_store(store):
    fs, truth, ref = store
    c, s, e = "chrT", 12_000, 13_000
    d = window_matrix(fs, c, s, e, seq=ref[c][s:e], min_frac=0.88)
    q = fs.query(c, s, e, min_frac=0.88)
    assert d["m6a"].shape == (len(q["row"]), e - s)
    seq = np.frombuffer(ref[c][s:e].encode(), np.uint8)
    at = np.isin(seq, np.frombuffer(b"AT", np.uint8))
    for i in range(len(q["row"])):
        row = d["m6a"][i]
        a, b = max(int(q["start"][i]) - s, 0), min(int(q["end"][i]) - s, e - s)
        assert d["coverage"][i, a:b].all() and not d["coverage"][i, :a].any()
        assert np.isnan(row[:a]).all() and np.isnan(row[b:]).all()
        inside = np.zeros(e - s, bool)
        inside[a:b] = True
        v = row[inside & ~at]                 # a mismatch base can still carry a call
        assert np.all(np.isnan(v) | (v == 1))
        called = np.zeros(e - s, bool)
        p = q["m6a"][i] - s
        called[p[(p >= 0) & (p < e - s)]] = True
        assert np.array_equal(row == 1, called)
        assert np.array_equal(row == 0, inside & at & ~called)
    assert d["positions"][0] == s and d["center"] == (s + e) // 2
    assert d["qname"][0].endswith("/ccs")


def test_save_window_loads_in_hier_hmm(store, tmp_path):
    hier_hmm = pytest.importorskip("hier_hmm")
    fs, _, ref = store
    c, s, e = "chrT", 12_000, 13_000
    kept = save_window(tmp_path / "w.npz", fs, c, s, e, seq=ref[c][s:e], min_frac=0.5)
    cfg = hier_hmm.load_config()
    ds = hier_hmm.load_npz(str(tmp_path / "w.npz"), cfg)
    assert ds.n_pos == e - s
    assert 0 < len(ds.fibers) <= kept["m6a"].shape[0]
    assert "qname" in ds.meta and "zmw" not in kept  # zmw not requested via columns
