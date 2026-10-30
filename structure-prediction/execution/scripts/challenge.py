"""The challenge compound set, and how run outputs are matched back to it.

One place owns the 20 blinded compounds, so ingest, submission building and gating
cannot disagree about which IDs exist or which SMILES belongs to each.
"""

from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEST_SET = ROOT / "data" / "cyp-challenge-TEST-BLINDED_structures.csv"

# The challenge expects exactly this many structures in a submission.
EXPECTED_COUNT = 20


def use_utf8_stdout() -> None:
    """Print UTF-8 regardless of the console code page.

    These scripts report chemistry: angstrom signs, dashes, the degree sign. A Windows
    console defaults to cp1252, which renders them as replacement characters even though
    every file this project writes is explicitly UTF-8. Called from each script's
    `main()` rather than on import, so importing the module has no side effects.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8")  # type: ignore[union-attr]
        except (AttributeError, OSError):
            pass  # a redirected or unusual stream; plain output is still fine


@dataclass(frozen=True)
class Compound:
    name: str  # e.g. "OCNT-2395436" — also the required PDB filename stem
    smiles: str


@lru_cache(maxsize=1)
def compounds() -> tuple[Compound, ...]:
    """The blinded structure-track compounds, in file order."""
    with open(TEST_SET, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    missing = {"Molecule_Name", "SMILES"} - set(rows[0] if rows else {})
    if missing:
        raise SystemExit(f"{TEST_SET.name}: missing column(s) {sorted(missing)}")
    out = tuple(Compound(r["Molecule_Name"].strip(), r["SMILES"].strip()) for r in rows)
    if len(out) != EXPECTED_COUNT:
        raise SystemExit(
            f"{TEST_SET.name}: expected {EXPECTED_COUNT} compounds, found {len(out)}"
        )
    return out


def smiles_map() -> dict[str, str]:
    return {c.name: c.smiles for c in compounds()}


def names() -> list[str]:
    return [c.name for c in compounds()]


def match_name(stem: str) -> str | None:
    """Recover a compound ID from an arbitrary filename stem.

    Tools name their output in their own style — the tutorial's Boltz-2 files are
    `cyp3a4_ocnt-2312554_model_0`, and docking exports may add a pose or variant
    suffix. Rather than guessing a convention, this searches the stem for any known
    compound ID, case-insensitively. Returns None if no ID matches, or if more than one
    does (which would make the intended compound ambiguous).
    """
    upper = stem.upper()
    hits = {name for name in names() if name in upper}
    if len(hits) != 1:
        # Try again allowing a missing separator, e.g. "OCNT2312554".
        loose = {
            name
            for name in names()
            if re.sub(r"[^A-Z0-9]", "", name) in re.sub(r"[^A-Z0-9]", "", upper)
        }
        hits = loose if len(loose) == 1 else hits
    return hits.pop() if len(hits) == 1 else None


def variant_tag(stem: str, name: str) -> str:
    """What a filename says about which variant of `name` it holds.

    `OCNT-2312554_R` -> `R`; `cyp3a4_OCNT-2312554_S_model_0` -> `S_MODEL_0`. Used to tell
    a stereoisomer pair apart from two unrelated poses of the same compound, so the
    choice between them can be stated rather than guessed. Returns "" when the stem is
    just the compound ID.
    """
    upper = stem.upper()
    index = upper.find(name)
    if index < 0:
        return ""
    return upper[index + len(name) :].strip("_-. ")


# A single undefined tetrahedral centre is expanded to these two configurations by
# scripts/enumerate_stereoisomers.py.
STEREO_TAGS = ("R", "S")


def match_files(paths: list[Path]) -> tuple[dict[str, Path], list[Path], dict[str, list[Path]]]:
    """Map compound ID -> file, plus the files that matched nothing or collided.

    Returns (matched, unmatched, duplicates). A duplicate means two files claim the
    same compound; the caller must resolve it rather than this picking silently, since
    choosing the wrong pose would quietly change what gets scored.
    """
    by_name: dict[str, list[Path]] = {}
    unmatched: list[Path] = []
    for path in sorted(paths):
        name = match_name(path.stem)
        if name is None:
            unmatched.append(path)
        else:
            by_name.setdefault(name, []).append(path)
    matched = {n: p[0] for n, p in by_name.items() if len(p) == 1}
    duplicates = {n: p for n, p in by_name.items() if len(p) > 1}
    return matched, unmatched, duplicates
