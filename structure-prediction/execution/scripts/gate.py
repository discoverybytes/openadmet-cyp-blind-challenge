"""Local gate for a submission zip: official validator + PoseBusters as the scorer runs it.

A submission is "ready" only once this passes (CLAUDE.md conventions). Two independent
reasons a pose scores zero on the leaderboard:

1. The zip or a PDB is malformed -> `vendor/.../structure_validation.py` rejects it.
2. The pose fails PoseBusters -> the scorer sets LDDT-PLI to 0 for that structure.

Step 2 is reproduced here rather than approximated, because the exact call matters: the
ligand's bond orders come from the challenge SMILES (not from the PDB, which has none),
and the protein context handed to PoseBusters includes *every* non-LIG residue, so the
heme is part of the clash and distance checks. That is mirrored from
`score_single_structure` in the vendored reference scorer (upstream lines ~723-776).

What this cannot tell us: the organisers' `POSEBUSTERS_MAX_FAILURES` is not public, so a
pose failing one check may or may not be zeroed. Every failure is therefore reported by
name and compound, and the caller decides (brief: Success criteria 2, open question).
"""

from __future__ import annotations

import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

_VENDOR = Path(__file__).resolve().parent.parent / "vendor"
if str(_VENDOR) not in sys.path:
    sys.path.insert(0, str(_VENDOR))

from openadmet_cyp_tutorial.posebusters_utils import (  # noqa: E402
    assign_bond_orders_from_template,
    bonded_ligand_copies,
    dedupe_ligand_copies,
)
from openadmet_cyp_tutorial.structure_validation import (  # noqa: E402
    validate_structure_submission,
)


# The numeric columns worth keeping from a full PoseBusters report. A relative distance
# is (actual distance) / (sum of radii); the check fails below CLASH_CUTOFF, so these
# say how much margin a pose has rather than just whether it passed. Measured
# thresholds and the reasoning behind each are in
# docs/findings/posebusters-crystal-reference.md.
CLASH_CUTOFF = 0.75

DETAIL_COLUMNS = {
    "rel_protein": "most_extreme_relative_distance_protein",
    "rel_organic": "most_extreme_relative_distance_organic_cofactors",
    "rel_inorganic": "most_extreme_relative_distance_inorganic_cofactors",
    "dist_protein": "smallest_distance_protein",
    "dist_organic": "smallest_distance_organic_cofactors",
    "dist_inorganic": "smallest_distance_inorganic_cofactors",
    "energy_ratio": "energy_ratio",
}


@dataclass
class PoseReport:
    """PoseBusters outcome for one compound's pose."""

    molecule_name: str
    failed_checks: list[str] = field(default_factory=list)
    n_checks: int = 0
    ligand_source: str = ""  # "template" (bond orders from SMILES) or "fallback"
    note: str = ""
    detail: dict[str, float] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not self.failed_checks and not self.note.startswith("ERROR")

    @property
    def tightest_margin(self) -> float | None:
        """Smallest relative distance across the three contact checks.

        Below `CLASH_CUTOFF` the pose is failing one of them; just above it, the pose is
        close enough that a small coordinate change could flip the verdict.
        """
        vals = [
            self.detail[k]
            for k in ("rel_protein", "rel_organic", "rel_inorganic")
            if self.detail.get(k) is not None and self.detail[k] == self.detail[k]
        ]
        return min(vals) if vals else None


@dataclass
class GateReport:
    """Everything the caller needs to decide whether to upload."""

    zip_path: Path
    validator_ok: bool = False
    validator_errors: list[str] = field(default_factory=list)
    poses: list[PoseReport] = field(default_factory=list)
    sha256: str = ""

    @property
    def posebusters_ok(self) -> bool:
        return all(p.ok for p in self.poses)

    @property
    def ok(self) -> bool:
        return self.validator_ok and self.posebusters_ok

    @property
    def failing(self) -> list[PoseReport]:
        return [p for p in self.poses if not p.ok]


def _quiet():
    """Silence RDKit/MDAnalysis chatter that would bury the actual findings."""
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.warning")
    warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


class _ForgivingTempdir:
    """`tempfile` stand-in whose TemporaryDirectory tolerates cleanup failures.

    The vendored validator extracts each PDB into a `TemporaryDirectory` and opens it
    with MDAnalysis, which keeps the file handle open. On Windows the directory then
    cannot be deleted, and the resulting `WinError 32` is swallowed by the validator's
    outer `except`, which reports the whole submission invalid — after every check has
    actually run and passed. The validator is vendored unmodified, so the fix belongs
    here: this changes only how the temporary directory is cleaned up, never a single
    validation check.
    """

    def __init__(self, real):
        self._real = real

    def TemporaryDirectory(self, *args, **kwargs):  # noqa: N802 - mirrors stdlib
        kwargs.setdefault("ignore_cleanup_errors", True)
        return self._real.TemporaryDirectory(*args, **kwargs)

    def __getattr__(self, item):
        return getattr(self._real, item)


def run_validator(
    zip_path: str | Path,
    expected_ids: set[str],
    expected_ligand_smiles: dict[str, str],
) -> tuple[bool, list[str]]:
    """Run the official validator with ID and ligand-connectivity checks enabled."""
    _quiet()
    from openadmet_cyp_tutorial import structure_validation

    original = structure_validation.tempfile
    structure_validation.tempfile = _ForgivingTempdir(original)  # type: ignore[assignment]
    try:
        return validate_structure_submission(
            Path(zip_path),
            expected_ids=expected_ids,
            expected_ligand_smiles=expected_ligand_smiles,
            require_lig_resname=True,
            check_numbering=True,
        )
    finally:
        structure_validation.tempfile = original  # type: ignore[assignment]


def _as_float(row, column: str) -> float:
    """Read one numeric cell as a float, with NaN for anything missing.

    A missing PoseBusters value can arrive as float NaN *or* as pandas' `NA` (for
    example the water-distance columns when a structure has no waters). The usual
    `x == x` NaN test raises `TypeError: boolean value of NA is ambiguous` on the
    latter, so conversion is attempted rather than guarded.
    """
    if column not in row:
        return float("nan")
    try:
        return float(row[column])
    except (TypeError, ValueError):
        return float("nan")


def _n_official_checks(buster) -> int:
    """How many leading columns of a full report are the scorer's actual verdicts.

    PoseBusters builds a full report as `test_columns + extra_columns` (see
    `posebusters.py::_make_table`), where `test_columns` are the config's
    `chosen_binary_test_output` entries — exactly the columns a default (non-full)
    report returns, and the only ones the challenge scorer reads. The extras are
    diagnostics, and some are booleans whose `False` means *good* (for example
    `most_extreme_clash_protein`). Counting those as failures would invent failures that
    the scorer never sees, so the verdict is taken from this many leading columns and
    the rest are used only as numbers.
    """
    return sum(
        len(module.get("chosen_binary_test_output", []))
        for module in buster.config["modules"]
    )


def _protein_context(pdb_path: Path):
    """Every non-LIG residue, combined — the `mol_cond` the scorer passes PoseBusters.

    Mirrors the scorer: it splits by PDB residue, drops `LIG`, and combines the rest.
    The heme stays in, so Fe contacts are inside the clash and minimum-distance checks.
    """
    from rdkit import Chem

    pred = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
    if pred is None:
        return None, None
    fragments = Chem.SplitMolByPDBResidues(pred)
    fragments.pop("LIG", None)
    protein = None
    for frag in fragments.values():
        protein = frag if protein is None else Chem.CombineMols(protein, frag)
    return pred, protein


def check_pose(pdb_path: str | Path, molecule_name: str, smiles: str) -> PoseReport:
    """Run PoseBusters on one complex the way the challenge scorer does."""
    from posebusters import PoseBusters
    from rdkit import Chem

    _quiet()
    pdb_path = Path(pdb_path)
    rep = PoseReport(molecule_name=molecule_name)

    pred, protein = _protein_context(pdb_path)
    if pred is None:
        rep.note = "ERROR: RDKit could not parse the PDB"
        return rep
    if protein is None:
        rep.note = "ERROR: no non-LIG residues found (no protein context)"
        return rep

    template = Chem.MolFromSmiles(smiles)
    if template is None:
        rep.note = f"ERROR: could not parse SMILES {smiles!r}"
        return rep

    # Preferred path: one bonded Mol per ligand copy, bond orders from the SMILES.
    copies = bonded_ligand_copies(str(pdb_path), smiles)
    if copies:
        rep.ligand_source = "template"
    else:
        # The scorer's fallback: dedupe to one representative copy and try the
        # template once. If that fails too, PoseBusters sees single bonds only.
        ligand = Chem.SplitMolByPDBResidues(pred).get("LIG")
        if ligand is None:
            rep.note = "ERROR: no LIG residue"
            return rep
        ligand = dedupe_ligand_copies(ligand)
        try:
            ligand = assign_bond_orders_from_template(template, ligand)
            rep.ligand_source = "fallback-template"
        except Exception as exc:  # noqa: BLE001
            rep.ligand_source = "fallback-single-bonds"
            rep.note = f"bond-order assignment failed ({exc})"
        copies = [ligand]

    buster = PoseBusters(config="dock")
    n_official = _n_official_checks(buster)
    # `full_report` adds the numeric columns behind each verdict, so the report can show
    # how much margin a pose has instead of only whether it passed.
    # The scorer keeps the best copy, so a compound passes if ANY copy passes.
    best: tuple[int, list[str], dict[str, float]] | None = None
    for copy in copies:
        result = buster.bust(
            mol_pred=copy, mol_true=template, mol_cond=protein, full_report=True
        )
        row = result.iloc[0]
        official = result.iloc[:, :n_official]
        bool_cols = [c for c in official.columns if official[c].dtype == bool]
        bools = official.iloc[0][bool_cols]
        failed = sorted(bools.index[~bools].tolist())
        detail = {key: _as_float(row, col) for key, col in DETAIL_COLUMNS.items()}
        if best is None or len(failed) < len(best[1]):
            best = (len(bool_cols), failed, detail)
    if best is None:
        # Only reachable if the ligand had no copies at all, which the branches above
        # are meant to prevent; reported rather than raised so one odd complex cannot
        # abort a 20-compound gate run.
        rep.note = "ERROR: no ligand copy could be checked"
        return rep
    rep.n_checks = best[0]
    rep.failed_checks = best[1]
    rep.detail = best[2]
    return rep


def sha256_of(path: str | Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gate(
    zip_path: str | Path,
    complexes: dict[str, Path],
    smiles: dict[str, str],
) -> GateReport:
    """Validate the zip, then PoseBusters every complex in it.

    Args:
        zip_path: the archive that would be uploaded.
        complexes: {molecule_name: path to that compound's complex PDB}.
        smiles: {molecule_name: challenge SMILES}.
    """
    report = GateReport(zip_path=Path(zip_path))
    report.sha256 = sha256_of(zip_path)
    report.validator_ok, report.validator_errors = run_validator(
        zip_path, set(smiles), smiles
    )
    for name in sorted(complexes):
        report.poses.append(check_pose(complexes[name], name, smiles[name]))
    return report


def format_report(report: GateReport) -> str:
    """Render a gate report as the markdown committed next to each submission."""
    lines = [
        f"# Gate report — `{report.zip_path.name}`",
        "",
        f"- SHA-256: `{report.sha256}`",
        f"- Validator: {'PASS' if report.validator_ok else 'FAIL'}",
        f"- PoseBusters: {'PASS' if report.posebusters_ok else 'FAIL'}"
        f" ({len(report.failing)} of {len(report.poses)} compounds failing)",
        f"- **Verdict: {'READY' if report.ok else 'NOT READY'}**",
        "",
    ]
    if report.validator_errors:
        lines += ["## Validator errors", ""]
        lines += [f"- {e}" for e in report.validator_errors] + [""]
    def fmt(v: float | None, places: int = 2) -> str:
        if v is None or v != v:
            return "—"
        return f"{v:.{places}f}"

    lines += [
        "## PoseBusters per compound",
        "",
        "Relative distances are (closest contact) / (sum of radii); below"
        f" **{CLASH_CUTOFF}** that check fails. `Fe` uses covalent radii, which is why"
        " real coordination passes; `heme` uses van der Waals radii and is the check"
        " that off-axis coordination trips. See"
        " `docs/findings/posebusters-crystal-reference.md`.",
        "",
        "A `Fe···lig` distance near 2.0-2.2 Å means the pose is coordinating the heme"
        " iron; 3 Å or more means it is not.",
        "",
        "| Compound | Failed | Fe···lig (Å) | rel. protein | rel. heme | rel. Fe | energy ratio | Ligand bonds | Note |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | :--- | :--- |",
    ]
    for p in report.poses:
        failed = ", ".join(f"`{c}`" for c in p.failed_checks) if p.failed_checks else "—"
        d = p.detail
        lines.append(
            f"| {p.molecule_name} | {failed} | {fmt(d.get('dist_inorganic'))} |"
            f" {fmt(d.get('rel_protein'))} |"
            f" {fmt(d.get('rel_organic'))} | {fmt(d.get('rel_inorganic'))} |"
            f" {fmt(d.get('energy_ratio'), 1)} | {p.ligand_source} |"
            f" {p.note or '—'} |"
        )
    lines.append("")
    tight = [
        p
        for p in report.poses
        if p.ok
        and p.tightest_margin is not None
        and p.tightest_margin < CLASH_CUTOFF + 0.05
    ]
    if tight:
        lines += [
            "### Passing, but close to a clash cutoff",
            "",
            "These poses pass now; a small coordinate change would fail them.",
            "",
        ]
        lines += [
            f"- {p.molecule_name}: tightest relative distance {p.tightest_margin:.3f}"
            for p in tight
        ]
        lines.append("")
    if not report.posebusters_ok:
        lines += [
            "Any failed check means the challenge scorer may set that structure's",
            "LDDT-PLI to 0. The organisers' tolerance is not public, so each failure",
            "needs explicit acceptance before upload.",
            "",
        ]
    return "\n".join(lines)
