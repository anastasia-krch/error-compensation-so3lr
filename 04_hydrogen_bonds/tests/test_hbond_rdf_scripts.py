#!/usr/bin/env python3
"""Small no-Lion test of the HCR/HCW groups and RDF pooling normalisation."""
from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / "scripts" / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    writer = load("write_lionanalysis_hbond_rdfs")
    pooler = load("pool_hbond_rdfs")
    assert writer.ids(751, 1000).split() == [str(i) for i in range(751, 1001)]
    assert writer.ids(1001, 1500).split() == [str(i) for i in range(1001, 1501)]
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "curve.txt"
        path.write_text("\n".join([
            "#r RDF iRDF PMF(kTunits) rawcount rawcountpertimestepperatom # from hcr to oacc; nSteps = 2 size(hcr) = 250 size(oacc) = 1000 numberdens. = 0.5",
            "0.01 0 0 0 3 0", "0.03 0 0 0 6 0", "",
        ]))
        radius, raw, steps, density = pooler.read_curve(path)
        result_r, result_g = pooler.pool([path], n_donors=250)
        shell = 4 * np.pi / 3 * ((result_r + 0.01) ** 3 - (result_r - 0.01) ** 3)
        assert np.array_equal(radius, result_r)
        assert np.array_equal(raw, [3.0, 6.0]) and steps == 2.0 and density == 0.5
        assert np.allclose(result_g, raw / (250 * shell * steps * density))
    print("HCR/HCW RDF config and pooling tests passed")


if __name__ == "__main__":
    main()
