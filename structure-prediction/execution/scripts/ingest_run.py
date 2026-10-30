"""Check and summarise a docking or co-folding run handed over for analysis.

Run this first on any new `runs/<run-id>/`. It answers: are all 20 compounds present,
is each complex usable (protein, heme, exactly one ligand), and what will need fixing
when the submission is built. It changes nothing except writing a draft `run.yaml`.

Why a separate step from building the submission: a run can be wrong in ways that are
cheap to fix at the source (a missing heme, a dropped compound, two poses for one
compound) and expensive to discover later. This reports those in one pass so the run
can be corrected or annotated before a submission slot is involved.

Usage (from the project root):
    .venv\\Scripts\\python scripts\\ingest_run.py runs\\20261005_tab5_2v0m-geometry
    .venv\\Scripts\\python scripts\\ingest_run.py runs\\<run-id> --ligand-resname MOL
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

import challenge
import pdbtools

STRUCTURE_SUFFIXES = (".pdb", ".cif", ".mmcif")


@dataclass
class ComplexSummary:
    """What one handed-over complex contains, and what is wrong with it."""

    name: str
    path: Path
    n_protein: int = 0
    n_heme: int = 0
    n_ligand: int = 0
    ligand_resname: str = ""
    ligand_residues: int = 0
    chains: str = ""
    # blocking: the submission cannot be built until this is resolved
    problems: list[str] = field(default_factory=list)
    # non-blocking: the builder normalises these automatically
    notes: list[str] = field(default_factory=list)


@dataclass
class LigandChoice:
    """Which records are the pose, and what the caller should be told about it."""

    atoms: list = field(default_factory=list)
    resname: str = ""
    problems: list[str] = field(default_factory=list)  # blocking
    notes: list[str] = field(default_factory=list)  # informational


def designate_ligand(cx: pdbtools.Complex, override: str | None) -> LigandChoice:
    """Decide which records are the ligand.

    Preference order: an explicit `--ligand-resname`, then an existing `LIG` residue,
    then a single unambiguous leftover group. Ambiguity is reported rather than
    resolved, because silently picking one of two candidate groups could submit a
    cofactor or a crystallisation additive as the ligand.
    """
    if override:
        picked = [a for a in cx.all_kept if a.resname == override.upper()]
        if not picked:
            return LigandChoice(problems=[f"no residue named {override.upper()}"])
        return LigandChoice(atoms=picked, resname=override.upper())

    if cx.ligand:
        return LigandChoice(atoms=cx.ligand, resname=pdbtools.LIGAND_RESNAME)

    groups: dict[str, list] = {}
    for atom in cx.other:
        groups.setdefault(atom.resname, []).append(atom)
    # Single-atom groups are ions, not the ligand we are predicting.
    candidates = {k: v for k, v in groups.items() if len(v) > 1}
    if len(candidates) == 1:
        resname, picked = next(iter(candidates.items()))
        return LigandChoice(
            atoms=picked,
            resname=resname,
            notes=[f"ligand taken from residue {resname}"],
        )
    if not candidates:
        return LigandChoice(problems=["no ligand candidate found"])
    sizes = ", ".join(f"{k} ({len(v)} atoms)" for k, v in sorted(candidates.items()))
    return LigandChoice(
        problems=[f"ambiguous ligand: {sizes} — pass --ligand-resname"]
    )


def summarise(
    name: str, path: Path, override: str | None, allow_no_heme: bool = False
) -> ComplexSummary:
    cx = pdbtools.read_complex(path)
    s = ComplexSummary(name=name, path=path)
    choice = designate_ligand(cx, override)
    ligand, resname = choice.atoms, choice.resname
    s.n_protein = len(cx.protein)
    s.n_heme = len(cx.heme)
    s.n_ligand = len(ligand)
    s.ligand_resname = resname
    s.chains = "".join(sorted(cx.chains())) or "—"

    s.problems.extend(choice.problems)
    s.notes.extend(choice.notes)
    if not cx.protein:
        s.problems.append("no protein records")
    if ligand:
        residues = {a.residue_key for a in ligand}
        s.ligand_residues = len(residues)
        if len(residues) > 1:
            # The validator allows at most one LIG residue; two copies of the pose
            # would be merged into one residue by the builder, which is wrong.
            s.problems.append(
                f"ligand spans {len(residues)} residues — expected one pose"
            )
    if not cx.heme:
        # Not fatal for the validator, but it changes the prediction: the scorer's
        # reference has a heme, and a pose placed in an empty porphyrin cavity is a
        # different structure. Blocking unless asked for, because the most likely cause
        # is an export that dropped the cofactor rather than a deliberate choice: a
        # docking tool's pose export may keep only the protein ATOM records.
        # See docs/findings/posebusters-crystal-reference.md.
        (s.notes if allow_no_heme else s.problems).append(
            "no heme (cofactor absent from the export)"
        )
    if cx.after_end:
        s.notes.append(f"{len(cx.after_end)} records after END (invisible to RDKit)")
    if cx.solvent:
        s.notes.append(f"{len(cx.solvent)} solvent atoms (will be dropped)")
    if resname and resname != pdbtools.LIGAND_RESNAME:
        s.notes.append(f"ligand is {resname}, not LIG (will be renamed)")
    if len(cx.chains()) > 3:
        s.notes.append(f"{len(cx.chains())} chains (will be reassigned to A/B/C)")

    dupe_names = 0
    for key in cx.residues():
        names = [a.name for a in cx.all_kept if a.residue_key == key]
        if len(names) != len(set(names)):
            dupe_names += 1
    if dupe_names:
        s.notes.append(f"{dupe_names} residues with duplicate atom names (will be fixed)")
    serials = [a.serial for a in cx.all_kept]
    if len(serials) != len(set(serials)):
        s.notes.append("duplicate atom serials (will be renumbered)")
    return s


def draft_run_yaml(run_dir: Path, summaries: list[ComplexSummary]) -> Path:
    """Write a `run.yaml` skeleton for the user to correct.

    Deliberately a draft: the settings that actually produced a run live in the
    docking or co-folding session, not in its output, so the fields that matter for
    reproducibility cannot be inferred here and are left as TODO.
    """
    path = run_dir / "run.yaml"
    if path.exists():
        return path
    heme = sum(1 for s in summaries if s.n_heme)
    body = f"""\
# Draft written by scripts/ingest_run.py on {date.today().isoformat()}.
# Correct every TODO: a leaderboard score is only evidence if the run is reproducible.
run_id: {run_dir.name}
process: TODO  # tab3-apf | tab5-vina | boltz2
variant: TODO  # the one thing this run changes vs the run it is compared to
date: {date.today().isoformat()}

molecule_2_commit: TODO  # git -C <docking tool> rev-parse --short HEAD at run time
receptor: TODO  # e.g. 2V0M chain A + HEM
settings:
  TODO: see the run log

complexes:
  found: {len(summaries)} of {challenge.EXPECTED_COUNT}
  with_heme: {heme}

notes: |
  TODO what this run was testing, and anything unusual about it.
"""
    path.write_text(body, encoding="utf-8")
    return path


def main() -> int:
    challenge.use_utf8_stdout()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="runs/<run-id>")
    parser.add_argument(
        "--ligand-resname",
        help="residue name holding the pose, if it is not LIG and not unambiguous",
    )
    parser.add_argument(
        "--allow-no-heme",
        action="store_true",
        help="accept complexes without a heme (a deliberate apo-cavity run, "
        "e.g. the tutorial's Boltz-2 baseline) instead of treating it as an error",
    )
    args = parser.parse_args()

    # Resolved so paths print relative to the project whether the caller passed a
    # relative or an absolute run directory.
    run_dir: Path = args.run_dir.resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"not a directory: {run_dir}")
    search = run_dir / "complexes" if (run_dir / "complexes").is_dir() else run_dir
    files = [p for p in search.rglob("*") if p.suffix.lower() in STRUCTURE_SUFFIXES]
    if not files:
        raise SystemExit(f"no {'/'.join(STRUCTURE_SUFFIXES)} files under {search}")

    matched, unmatched, duplicates = challenge.match_files(files)
    print(f"{run_dir.name}: {len(files)} structure file(s) in {search}")
    print(
        f"matched {len(matched)} of {challenge.EXPECTED_COUNT} compounds"
        + (f", {len(duplicates)} with several candidates" if duplicates else "")
        + "\n"
    )

    # The expected case: a compound whose stereocentre the challenge left undefined was
    # docked as both enantiomers, so one has to be chosen for the submission. Kept
    # separate from a genuine ambiguity, and from `missing` — these compounds are
    # present, they just need a decision.
    stereo = {
        n: p
        for n, p in duplicates.items()
        if {challenge.variant_tag(x.stem, n) for x in p} >= set(challenge.STEREO_TAGS)
    }
    other = {n: p for n, p in duplicates.items() if n not in stereo}

    if duplicates:
        if stereo:
            print(
                f"STEREOISOMER PAIRS — {len(stereo)} compound(s) docked as both "
                "enantiomers:"
            )
            for name, paths in sorted(stereo.items()):
                tags = ", ".join(
                    sorted(challenge.variant_tag(p.stem, name) for p in paths)
                )
                print(f"  {name}: {tags}")
            print(
                "  Only one configuration can be submitted per compound ID. Choose at"
                "\n  build time, e.g. "
                + " ".join(f"--pick {n}=R" for n in sorted(stereo))
                + "\n"
            )
        if other:
            print("AMBIGUOUS — more than one file per compound:")
            for name, paths in sorted(other.items()):
                print(f"  {name}: {', '.join(p.name for p in paths)}")
            print()
    if unmatched:
        print(f"UNMATCHED — no compound ID in the filename ({len(unmatched)}):")
        for p in unmatched[:10]:
            print(f"  {p.name}")
        if len(unmatched) > 10:
            print(f"  ... and {len(unmatched) - 10} more")
        print()

    # A compound with candidate files is not missing; it is awaiting a choice.
    missing = [
        n for n in challenge.names() if n not in matched and n not in duplicates
    ]
    if missing:
        print(f"MISSING — no file for {len(missing)} compound(s):")
        print("  " + ", ".join(missing))
        print()

    summaries = [
        summarise(name, matched[name], args.ligand_resname, args.allow_no_heme)
        for name in sorted(matched)
    ]

    print(f"{'Compound':<16} {'chains':<7} {'prot':>5} {'heme':>5} {'lig':>4}  issues")
    for s in summaries:
        issues = "; ".join(s.problems) or ("ok" if not s.notes else "")
        print(
            f"{s.name:<16} {s.chains:<7} {s.n_protein:>5} {s.n_heme:>5} "
            f"{s.n_ligand:>4}  {issues}"
        )

    noted = [s for s in summaries if s.notes]
    if noted:
        print("\nThe submission builder will normalise these:")
        for s in noted:
            for note in s.notes:
                print(f"  {s.name}: {note}")

    blocking = [s for s in summaries if s.problems]
    yaml_path = draft_run_yaml(run_dir, summaries)
    try:
        shown = yaml_path.relative_to(challenge.ROOT)
    except ValueError:  # a run directory outside the project
        shown = yaml_path
    print(f"\nrun manifest: {shown}")
    if yaml_path.read_text(encoding="utf-8").count("TODO"):
        print("  (draft — fill in the TODOs before this run's score is used as evidence)")

    if blocking or missing or duplicates or unmatched:
        counts = [
            f"{len(blocking)} complex(es) with problems" if blocking else "",
            f"{len(missing)} missing" if missing else "",
            f"{len(stereo)} awaiting a stereo choice" if stereo else "",
            f"{len(other)} ambiguous" if other else "",
            f"{len(unmatched)} unmatched file(s)" if unmatched else "",
        ]
        print("\nNOT READY to build: " + ", ".join(c for c in counts if c) + ".")
        if stereo and not (blocking or missing or other or unmatched):
            print("  Everything else is fine — pass --pick at build time.")
        return 1
    print(f"\nREADY to build: all {len(summaries)} complexes usable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
