#!/usr/bin/env python3
"""Shared pieces for the charge chapter: ion fragments and AMS result parsing.

Imported by the sibling compute_*.py scripts, the same way 01_density's
build_density_table.py imports extract_density. Nothing here is needed to plot: the
notebook reads only the committed CSVs.

All methods, including SO3LR through AMS's ASE engine, are read from AMS output.  Heavy
dependencies are therefore limited to ASE/scipy for preparing clusters; plotting and
post-processing never import JAX, mlff, or so3lr directly.
"""
from __future__ import annotations

import pathlib
import re

import numpy as np

ANION_SIZE = 15      # NTf2- : C F3 S N O2 S O2 C F3
CATION_SIZE = 19     # EMIM+
SR_CUTOFF_A = 4.5    # the model's short-range cutoff
N_LAYERS = 3         # message passes; receptive field is the product, 13.5 A
LR_CUTOFF_A = 12.0   # where the model force-shifts its electrostatics to zero
DEBYE_PER_E_BOHR = 2.541746473


def detect_ion_fragments(atoms):
    """Group atom indices into per-ion fragments by PBC-aware bond detection.

    Raises unless the frame decomposes into exactly the two expected fragment sizes, so a
    mis-detected topology cannot pass silently into a charge average.
    """
    from ase.neighborlist import natural_cutoffs, neighbor_list
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    n = len(atoms)
    i, j = neighbor_list("ij", atoms, cutoff=natural_cutoffs(atoms, mult=1.25))
    _, labels = connected_components(
        coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)), directed=False)

    fragments: dict[int, list] = {}
    for label in np.unique(labels):
        index = np.where(labels == label)[0]
        fragments.setdefault(len(index), []).append(index)
    if sorted(fragments) != [ANION_SIZE, CATION_SIZE]:
        raise RuntimeError(f"expected fragment sizes {[ANION_SIZE, CATION_SIZE]}, "
                           f"found {sorted(fragments)}")
    return fragments[ANION_SIZE], fragments[CATION_SIZE]


def unwrap_fragment(positions, fragment, cell_diagonal):
    """Make one fragment whole across the periodic boundary, anchored on its first atom."""
    reference = positions[fragment[0]]
    delta = positions[fragment] - reference
    delta -= np.round(delta / cell_diagonal) * cell_diagonal
    return reference + delta


def parse_ams_output(path: pathlib.Path) -> dict | None:
    """Dipole, energy and atomic charges from one AMS stdout log.

    AMS sends ADF's readable output to stdout, so the scheduler log IS the data file --
    `ams.results/ams.log` holds only the terse progress trace.

    Returns None unless the run terminated normally, so a crashed or half-written job can
    never contribute a number to a table. The MDC table is absent for force-field engines,
    which is not a failure: those keys come back None and the caller must not ask for them.
    """
    if not path.exists():
        return None
    lines = path.read_text(errors="replace").splitlines()
    if not any("NORMAL TERMINATION" in line for line in lines):
        return None

    # Take the LAST occurrence of each block: it belongs to the final geometry.
    symbols, mdc = [], []
    heads = [n for n, line in enumerate(lines) if "MDC-m" in line and "Level" in line]
    if heads:
        for line in lines[heads[-1] + 2:]:
            parts = line.split()
            if len(parts) < 5 or not parts[0].isdigit():
                break
            symbols.append(parts[1])
            mdc.append([float(parts[2]), float(parts[3]), float(parts[4])])

    mulliken, mulliken_symbols = [], []
    heads = [n for n, line in enumerate(lines) if line.strip() == "Atomic charges"]
    if heads:
        for line in lines[heads[-1] + 3:]:
            parts = line.split()
            if len(parts) != 3 or not parts[0].isdigit():
                break
            mulliken_symbols.append(parts[1])
            mulliken.append(float(parts[2]))
    if not symbols:
        symbols = mulliken_symbols        # force fields: the only source of symbols
    if not symbols:
        return None

    energy = np.nan
    for line in lines:
        if line.strip().startswith("Energy (hartree)"):
            try:
                energy = float(line.split()[-1])
            except ValueError:
                pass
    magnitudes = []
    for line in lines:
        match = re.search(r"magnitude\s*(?:\(Debye\))?\s*:\s*([-+0-9.eEdD]+)",
                          line, flags=re.IGNORECASE)
        if match:
            magnitudes.append(float(match.group(1).replace("D", "E").replace("d", "e")))
    result = dict(
        symbols=symbols,
        energy_hartree=energy,
        dipole_D=magnitudes[-1] if magnitudes else np.nan,
        mdc_q=np.array([row[2] for row in mdc]) if mdc else None,
        # ADF's "Atomic charges" are Mulliken charges.  For Engine ASE, the same
        # table holds the calculator-provided values (SO3LR partial charges).
        atomic_charges=(np.array(mulliken) if len(mulliken) == len(symbols) else None),
        mulliken=(np.array(mulliken) if len(mulliken) == len(symbols) else None))
    # Engine ASE prints charges rounded to four decimals.  The driver copies its RKF
    # beside the stdout log, giving the analysis the exact charge vector and energy
    # returned by AMS without treating human-readable output as a data format.
    rkf = path.with_suffix(".rkf")
    if rkf.is_file():
        try:
            from scm.plams import KFFile

            values = KFFile(str(rkf)).read_section("AMSResults")
            if "Charges" in values:
                result["atomic_charges"] = np.asarray(values["Charges"], dtype=float)
                result["mulliken"] = result["atomic_charges"]
            if "Energy" in values:
                result["energy_hartree"] = float(values["Energy"])
            if "DipoleMoment" in values:
                result["dipole_D"] = float(np.linalg.norm(values["DipoleMoment"]) * DEBYE_PER_E_BOHR)
        except Exception as exc:
            raise RuntimeError(f"could not read AMS ASE result {rkf}: {exc}") from exc
    return result


def index_ams_logs(logdir: pathlib.Path) -> dict[tuple[str, str], dict]:
    """Index every AMS result under `logdir` by (tag, level of theory).

    Identification must include the level: the same geometry is run at several
    functional/basis combinations, so matching on the tag alone silently returns whichever
    log the glob happened to reach first. Logs without an explicit `tag=`/`species=`
    provenance line are skipped rather than guessed at.
    """
    index: dict[tuple[str, str], dict] = {}
    for candidate in sorted(pathlib.Path(logdir).glob("*.out")):
        text = candidate.read_text(errors="replace")
        tag = re.search(r"(?:species|tag)=(\S+)", text)
        if tag is None:
            continue
        xc = re.search(r"xc=(\S+)", text)
        basis = re.search(r"basis=(\S+)", text)
        level = f"{xc.group(1) if xc else 'PBE'}/{basis.group(1) if basis else 'TZP'}"
        parsed = parse_ams_output(candidate)
        if parsed is not None:
            index[(tag.group(1), level)] = parsed
    return index
