"""Assemble a run into a submission zip and gate it before upload.

Takes a run that `ingest_run.py` accepted, rewrites each complex into the format the
challenge requires, zips the 20 files, and runs the gate (official validator +
PoseBusters as the scorer calls it). Nothing is uploaded: the script prints the zip
path and its verdict, and John uploads through the portal (decision 0001).

Normalisation is deliberate rather than best-effort. Every rewrite here corresponds to
a rule that would otherwise cost the whole structure its score (brief: Risk 5):
`LIG` as the ligand residue name, at most three chains, unique atom serials and names,
no records stranded after `END`, and no solvent.

Compounds with an undefined stereocentre arrive as an `_R`/`_S` pair (see
`scripts/enumerate_stereoisomers.py`), and only one configuration can be submitted per
compound ID. `--pick` states which, so the choice is explicit and reproducible rather
than decided by sort order.

Usage (from the project root):
    .venv\\Scripts\\python scripts\\build_submission.py runs\\20261005_tab5_2v0m-geometry
    .venv\\Scripts\\python scripts\\build_submission.py runs\\<run-id> --label v2 --allow-no-heme
    .venv\\Scripts\\python scripts\\build_submission.py runs\\<run-id> \\
        --pick OCNT-2312554=R --pick OCNT-2395415=S --pick OCNT-2312966=R
"""

from __future__ import annotations

import argparse
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import challenge
import gate
import ingest_run
import pdbtools

ROOT = challenge.ROOT
SUBMISSIONS = ROOT / "submissions"


def normalise(
    name: str, source: Path, out_dir: Path, override: str | None
) -> tuple[Path, list[str]]:
    """Write one compound's complex in submission format.

    Returns the path and any notes about what was changed.
    """
    cx = pdbtools.read_complex(source)
    choice = ingest_run.designate_ligand(cx, override)
    if not choice.atoms:
        raise SystemExit(f"{name}: {'; '.join(choice.problems) or 'no ligand found'}")
    notes: list[str] = []

    # Ligand hydrogens are dropped. The official validator rebuilds the ligand with
    # `MolFromPDBFile` and no CONECT records, so bonds are perceived from interatomic
    # distances; with explicit hydrogens present that perception misassigns bonds and
    # the submission is rejected as "Ligand connectivity does not match expected
    # SMILES" (measured on 8 of 20 refined complexes). Nothing is lost: the scorer
    # selects `rname=LIG and ele!=H` for both LDDT-PLI and BiSyRMSD, and PoseBusters'
    # distance checks ignore hydrogens, so they are never scored.
    ligand = [a for a in choice.atoms if a.element != "H"]
    if len(ligand) != len(choice.atoms):
        notes.append(f"dropped {len(choice.atoms) - len(ligand)} ligand hydrogens")
    if not ligand:
        raise SystemExit(f"{name}: ligand has no heavy atoms")

    # Some exports write the element from the atom name's first letter, which turns the
    # heme iron's FE into F. Left alone, the submitted file declares a fluorine where the
    # iron belongs and every consumer that reads the element believes it.
    heme, fixed = pdbtools.fix_heme_iron_element(cx.heme)
    if fixed:
        notes.append(f"corrected {fixed} heme iron element (was not FE)")

    # The heme is kept as its own chain; anything else non-protein (ions, additives,
    # co-folding artefacts) is left out, since the challenge expects protein + heme +
    # ligand and extra chains breach the three-chain limit.
    ligand_keys = {a.residue_key for a in choice.atoms}
    extras = [
        a
        for a in cx.other
        if a.residue_key not in ligand_keys and a.resname != choice.resname
    ]
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
    header = [
        f"REMARK 300 {name} for the OpenADMET CYP structure track",
        f"REMARK 300 source {source.name}",
        f"REMARK 300 normalised by scripts/build_submission.py at {stamp}",
    ]
    if extras:
        dropped = sorted({a.resname for a in extras})
        notes.append(f"dropped non-protein residues: {', '.join(dropped)}")
        header.append(f"REMARK 300 dropped non-protein residues: {', '.join(dropped)}")
    if notes:
        header.append(f"REMARK 300 {'; '.join(notes)}")
    path = pdbtools.write_complex(
        out_dir / f"{name}.pdb",
        protein=cx.protein,
        heme=heme,
        ligand=ligand,
        header=header,
    )
    return path, notes


def resolve_choices(
    duplicates: dict[str, list[Path]], choices: dict[str, str]
) -> tuple[dict[str, Path], list[str]]:
    """Apply `--pick NAME=TAG` to compounds that have more than one candidate file.

    The three compounds with an undefined stereocentre arrive as an `_R`/`_S` pair, so
    exactly one configuration has to be submitted. Requiring the choice on the command
    line keeps it in the shell history and in the run manifest, and makes the submission
    reproducible — unlike resolving it by sort order.
    """
    resolved: dict[str, Path] = {}
    unresolved: list[str] = []
    for name, paths in sorted(duplicates.items()):
        tags = {challenge.variant_tag(p.stem, name): p for p in paths}
        wanted = choices.get(name)
        if wanted and wanted.upper() in tags:
            resolved[name] = tags[wanted.upper()]
        else:
            shown = ", ".join(f"{t or '(none)'} -> {p.name}" for t, p in sorted(tags.items()))
            hint = (
                f"--pick {name}=R  (or =S)"
                if set(tags) >= set(challenge.STEREO_TAGS)
                else f"--pick {name}=<tag>"
            )
            unresolved.append(f"{name}: {shown}\n      choose with {hint}")
    return resolved, unresolved


def main() -> int:
    challenge.use_utf8_stdout()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="runs/<run-id>")
    parser.add_argument(
        "--label", default="", help="suffix for the zip name, e.g. a variant tag"
    )
    parser.add_argument(
        "--pick",
        action="append",
        default=[],
        metavar="NAME=TAG",
        help="which variant to submit for a compound with several candidate files, "
        "e.g. --pick OCNT-2312554=R for the R enantiomer (repeatable)",
    )
    parser.add_argument("--ligand-resname", help="residue name holding the pose")
    parser.add_argument(
        "--allow-no-heme",
        action="store_true",
        help="build even if complexes have no heme (a deliberate apo-cavity run)",
    )
    parser.add_argument(
        "--skip-posebusters",
        action="store_true",
        help="validator only — for a fast format check, never for an upload decision",
    )
    args = parser.parse_args()

    run_dir: Path = args.run_dir.resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"not a directory: {run_dir}")
    search = run_dir / "complexes" if (run_dir / "complexes").is_dir() else run_dir
    files = [
        p
        for p in search.rglob("*")
        if p.suffix.lower() in ingest_run.STRUCTURE_SUFFIXES
    ]
    matched, unmatched, duplicates = challenge.match_files(files)

    choices: dict[str, str] = {}
    for item in args.pick:
        name, _, tag = item.partition("=")
        if not tag:
            raise SystemExit(f"--pick needs NAME=TAG, got {item!r}")
        if name not in set(challenge.names()):
            raise SystemExit(f"--pick: unknown compound {name!r}")
        choices[name] = tag
    picked, unresolved = resolve_choices(duplicates, choices)
    matched.update(picked)

    # A compound with candidate files is awaiting a choice, not missing.
    missing = [
        n for n in challenge.names() if n not in matched and n not in duplicates
    ]
    if missing or unresolved:
        print(f"Cannot build {run_dir.name}:")
        if missing:
            print(f"  missing {len(missing)}: {', '.join(missing)}")
        for line in unresolved:
            print(f"  several candidates for {line}")
        print("  run scripts/ingest_run.py for the full picture.")
        return 1
    if unmatched:
        print(f"note: ignoring {len(unmatched)} file(s) with no compound ID\n")

    staged = run_dir / "submission"
    for old in staged.glob("*.pdb"):
        old.unlink()
    results = {
        name: normalise(name, matched[name], staged, args.ligand_resname)
        for name in challenge.names()
    }
    built = {name: path for name, (path, _) in results.items()}
    all_notes = sorted({note for _, notes in results.values() for note in notes})
    if all_notes:
        print("normalisation applied: " + "; ".join(all_notes))
    if picked:
        # Each normalised PDB records its source filename in a REMARK, so the chosen
        # configuration stays readable in the submitted file itself.
        print("variant chosen per compound:")
        for name in sorted(picked):
            print(f"  {name}: {picked[name].name}")

    no_heme = [n for n, p in built.items() if not pdbtools.read_complex(p).heme]
    if no_heme and not args.allow_no_heme:
        print(
            f"Cannot build: {len(no_heme)} complex(es) have no heme "
            f"({', '.join(no_heme[:5])}{'...' if len(no_heme) > 5 else ''}).\n"
            "  The scorer's reference has one, so this changes the prediction.\n"
            "  Re-export with the cofactor, or pass --allow-no-heme deliberately."
        )
        return 1

    SUBMISSIONS.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.label}" if args.label else ""
    # Named after the run, which already starts with its date; the build timestamp and
    # SHA-256 live in the gate report, so a second date here would only be redundant.
    zip_path = SUBMISSIONS / f"{run_dir.name}{suffix}.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name in challenge.names():
            zf.write(built[name], arcname=f"{name}.pdb")
    print(f"built {zip_path.relative_to(ROOT)} ({len(built)} complexes)")

    smiles = challenge.smiles_map()
    if args.skip_posebusters:
        ok, errors = gate.run_validator(zip_path, set(smiles), smiles)
        print(f"validator: {'PASS' if ok else 'FAIL'}")
        for e in errors:
            print(f"  {e}")
        print("PoseBusters skipped — not an upload decision.")
        return 0 if ok else 1

    print("gating (validator + PoseBusters on 20 complexes, a few minutes)...")
    report = gate.gate(zip_path, built, smiles)
    report_path = zip_path.with_suffix(".gate.md")
    report_path.write_text(gate.format_report(report), encoding="utf-8")
    zip_path.with_suffix(".sha256").write_text(
        f"{report.sha256}  {zip_path.name}\n", encoding="utf-8"
    )

    print(f"\nvalidator:   {'PASS' if report.validator_ok else 'FAIL'}")
    for e in report.validator_errors[:10]:
        print(f"  {e}")
    print(
        f"PoseBusters: {'PASS' if report.posebusters_ok else 'FAIL'}"
        f" ({len(report.failing)} of {len(report.poses)} failing)"
    )
    for p in report.failing:
        print(f"  {p.molecule_name}: {', '.join(p.failed_checks) or p.note}")
    print(f"\nreport: {report_path.relative_to(ROOT)}")
    print(f"sha256: {report.sha256}")

    if report.ok:
        print(f"\nREADY — upload {zip_path.relative_to(ROOT)} via the challenge portal.")
        return 0
    print(
        "\nNOT READY. Any failed PoseBusters check may zero that structure's LDDT-PLI;"
        "\nthe organisers' tolerance is unknown. Accept each failure explicitly, or fix"
        "\nthe run and rebuild."
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
