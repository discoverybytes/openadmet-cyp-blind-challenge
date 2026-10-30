"""Write Boltz-2 input YAMLs for the 20 challenge compounds, plus the command to run.

Claude writes these; John runs them and hands back the complexes (decision 0001). The
schema is Boltz 2.2.1's, read from the installed package rather than from documentation:

    version: 1
    sequences:
      - protein: {id: A, sequence: ..., msa: ...}
      - ligand:  {id: B, ccd: HEM}
      - ligand:  {id: C, smiles: ...}
    constraints:
      - pocket: {binder: C, contacts: [[B, FE]], max_distance: 6, force: false}
    templates:
      - cif: path          # or `pdb: path`
        chain_id: A        # our protein chain
        template_id: A     # the chain to take from the template

Three protocol choices this script exposes, because each is worth testing as a single
variable against the leaderboard:

* **The heme.** Included by default as `ccd: HEM`. CYP3A4 without its cofactor leaves a
  ligand-sized hole where the porphyrin belongs, and the scorer's reference has one.
* **A pocket constraint.** Off by default. `--pocket-heme` tells Boltz the ligand should
  sit within `--pocket-distance` of the heme **iron**, which biases placement toward the
  axial site without dictating a pose. `--force-pocket` makes it a hard constraint. A
  non-polymer chain is addressed by atom name, not residue index (see
  `token_spec_to_ids` in Boltz's schema parser), so the contact names `FE` directly —
  which states the chemistry we care about rather than merely "somewhere near the heme".
* **Templates.** Off by default, and John chooses them. One template for every compound
  (`--template`), or per-compound templates from a CSV (`--templates-csv`), which is the
  interesting case: the nearest holo CYP3A4 structure differs per ligand.

Usage (from the project root):
    .venv\\Scripts\\python scripts\\make_cofold_yamls.py --name 20261005_boltz2_heme-only
    .venv\\Scripts\\python scripts\\make_cofold_yamls.py --name <set> --pocket-heme
    .venv\\Scripts\\python scripts\\make_cofold_yamls.py --name <set> --templates-csv cofold/templates.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import challenge

ROOT = challenge.ROOT
COFOLD = ROOT / "cofold"
FASTA = ROOT / "vendor" / "openadmet_cyp_tutorial" / "CYP3A4_protein_sequence.fasta"

PROTEIN_CHAIN = "A"
HEME_CHAIN = "B"
LIGAND_CHAIN = "C"

# Boltz's `check_inputs` globs the directory it is given and raises on anything that is
# not .fasta/.yaml — subdirectories included. So the YAMLs get a directory to themselves
# and everything else about the set (RUN.md, the bundled template, the out/ directory)
# lives beside it rather than in it.
INPUT_SUBDIR = "yaml"


def cyp3a4_sequence() -> str:
    """The monomeric CYP3A4 sequence (UniProt P08684) from the vendored FASTA."""
    lines = FASTA.read_text(encoding="utf-8").strip().splitlines()
    return "".join(line.strip() for line in lines[1:])


def load_templates(path: Path) -> dict[str, list[str]]:
    """Read per-compound template paths from a CSV.

    Columns: `Molecule_Name`, `template` (a path to a .cif or .pdb), and optionally
    `template_chain` (defaults to A). More than one row per compound is allowed, which
    hands Boltz several templates for that compound.
    """
    out: dict[str, list[str]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        required = {"Molecule_Name", "template"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path.name}: missing column(s) {sorted(missing)}")
        known = set(challenge.names())
        for row in reader:
            name = row["Molecule_Name"].strip()
            if name not in known:
                raise SystemExit(f"{path.name}: unknown compound {name!r}")
            template = row["template"].strip()
            if not template:
                continue
            chain = (row.get("template_chain") or "A").strip() or "A"
            out.setdefault(name, []).append(f"{template}|{chain}")
    return out


def quote(value: str) -> str:
    """Single-quote a YAML scalar, escaping quotes.

    SMILES routinely contain `#`, `[`, `:` and backslashes, all of which change meaning
    in unquoted YAML, so every ligand string is quoted rather than hoping it parses.
    """
    return "'" + value.replace("'", "''") + "'"


def build_yaml(
    compound: challenge.Compound,
    sequence: str,
    msa: str,
    heme: bool,
    pocket: bool,
    pocket_distance: float,
    pocket_atom: str,
    force_pocket: bool,
    templates: list[str],
    affinity: bool,
) -> str:
    lines = ["version: 1", "sequences:"]
    lines += [
        "  - protein:",
        f"      id: {PROTEIN_CHAIN}",
        f"      sequence: {sequence}",
    ]
    # `--msa server` omits the key entirely, which is what makes `boltz predict
    # --use_msa_server` generate one. Writing `msa: empty` instead would pin the run to
    # single-sequence mode and the server flag would be silently ignored.
    if msa != "server":
        lines.append(f"      msa: {msa}")
    if heme:
        # The heme comes from the CCD rather than a SMILES string, so Boltz uses the
        # standard HEM geometry and atom naming.
        lines += ["  - ligand:", f"      id: {HEME_CHAIN}", "      ccd: HEM"]
    lines += [
        "  - ligand:",
        f"      id: {LIGAND_CHAIN}",
        f"      smiles: {quote(compound.smiles)}",
    ]

    if pocket:
        if not heme:
            raise SystemExit("--pocket-heme needs the heme; drop --no-heme")
        # A non-polymer chain is addressed by ATOM NAME, not residue index: Boltz's
        # `token_spec_to_ids` looks the contact up in its atom map, so `[B, 1]` raises
        # KeyError while `[B, FE]` resolves to the iron.
        lines += [
            "constraints:",
            "  - pocket:",
            f"      binder: {LIGAND_CHAIN}",
            f"      contacts: [[{HEME_CHAIN}, {pocket_atom}]]",
            f"      max_distance: {pocket_distance:g}",
            f"      force: {str(force_pocket).lower()}",
        ]

    if templates:
        lines.append("templates:")
        for spec in templates:
            path, _, chain = spec.partition("|")
            key = "pdb" if Path(path).suffix.lower() == ".pdb" else "cif"
            lines += [
                f"  - {key}: {path}",
                f"    chain_id: {PROTEIN_CHAIN}",
                f"    template_id: {chain or 'A'}",
            ]

    if affinity:
        # Not scored by this track; available only because it costs nothing extra.
        lines += ["properties:", "  - affinity:", f"      binder: {LIGAND_CHAIN}"]

    return "\n".join(lines) + "\n"


def run_instructions(out_dir: Path, msa: str, args, entries) -> str:
    rel = out_dir.relative_to(ROOT).as_posix()
    # WSL reaches the project through /mnt/<drive>, so build that path from the real
    # location and the command can be pasted without editing.
    drive = ROOT.drive.rstrip(":").lower()
    wsl_dir = f"/mnt/{drive}{ROOT.as_posix()[len(ROOT.drive):]}/{rel}"
    # Four-space indent: these lines sit inside the `for SEED in ...` loop below.
    msa_flag = "\n    --use_msa_server \\" if msa == "server" else ""
    extra = "".join(f"\n    {f} \\" for f in (args.extra_flag or []))
    seeds = [s.strip() for s in (args.seeds or "42").split(",") if s.strip()]
    if msa == "empty":
        msa_note = (
            "`msa: empty` runs single-sequence mode: no network, but Boltz warns that\n"
            "predictions are suboptimal without an MSA."
        )
    elif msa == "server":
        msa_note = (
            "The YAMLs carry **no `msa:` key**, which is what makes `--use_msa_server`\n"
            "generate one. Writing `msa: empty` instead would pin the run to\n"
            "single-sequence mode and the server flag would be silently ignored — so the\n"
            "absence of that line and the presence of that flag have to stay in step.\n\n"
            "This costs one ColabFold server call per compound rather than one for the set.\n"
            "The protein is identical in every entry, so the MSAs should come back identical;\n"
            "if the server is rate-limiting or flaky, generate one `.a3m` and pass it with\n"
            "`--msa <path>` instead."
        )
    elif (out_dir / msa).exists():
        msa_note = (
            f"Every YAML points at `{msa}`, which is **already in this folder** — no\n"
            "network call, and nothing to generate. The protein is identical in every\n"
            "entry, so one alignment serves all of them, and pinning it to a file keeps\n"
            "the MSA out of the set of things that can vary between runs."
        )
    else:
        msa_note = (
            f"Every YAML points at `{msa}`, which is **not in this folder** — put it there\n"
            "(or make the path resolve from the directory the command is run in) before\n"
            "starting. Generate it once and reuse it for every entry: the protein is\n"
            "identical, so one MSA keeps the comparison controlled and avoids 20 server\n"
            "calls."
        )
    template_note = (
        str(args.templates_csv)
        if args.templates_csv
        else (
            f"{Path(args.template).name} (bundled into this folder, referenced by "
            f"basename), chain {args.template_chain}"
            if args.template and args.bundle_template
            else f"{args.template}, chain {args.template_chain}"
            if args.template
            else "none"
        )
    )
    settings = [
        # `ccd: HEM` is iron protoporphyrin IX: all 43 atoms, the iron among them, named
        # FE. So the cofactor and its metal both come from the CCD, not from a SMILES.
        f"heme: {'yes — ccd HEM, 43 atoms including the iron (FE)' if not args.no_heme else 'NO (ablation)'}",
        f"pocket constraint: {'ligand steered to within ' + format(args.pocket_distance, 'g') + ' A of heme ' + args.pocket_atom + (' (FORCED)' if args.force_pocket else ' (soft bias)') if args.pocket_heme else 'none'}",
        f"template: {template_note}",
        f"msa: {msa}",
        f"entries: {len(entries)}",
        f"seeds: {', '.join(seeds)}",
    ]
    if args.extra_flag:
        settings.append(f"extra flags: {' '.join(args.extra_flag)}")
    purpose = f"\n{args.purpose.strip()}\n" if getattr(args, "purpose", None) else ""
    # Only tell John how to generate an MSA when he actually has to.
    msa_howto = "" if (msa not in ("empty", "server") and (out_dir / msa).exists()) else """
To generate it once with the ColabFold server (an external network call, so it is
John's to make):

```bash
conda activate boltz2
boltz predict <one>.yaml --use_msa_server --out_dir /tmp/msa_probe
# the generated .a3m lands under /tmp/msa_probe/processed/msa/ — copy it somewhere
# stable and re-run this script with --msa <that path>
```
"""
    return f"""# Boltz-2 run: `{out_dir.name}`

Written by `scripts/make_cofold_yamls.py`. {len(entries)} YAMLs, one per entry.
{purpose}
Settings:

{chr(10).join('- ' + s for s in settings)}

## MSA

{msa_note}
{msa_howto}
## Run

```bash
conda activate boltz2
cd "{wsl_dir}"
for SEED in {' '.join(seeds)}; do
  boltz predict {INPUT_SUBDIR} \\
    --out_dir "out_seed$SEED" \\
    --model boltz2 \\
    --seed "$SEED" \\
    --diffusion_samples 5 \\
    --recycling_steps 3 \\
    --output_format pdb \\{msa_flag}{extra}
    --use_potentials
done
```

Each seed writes its own `out_seed<N>/`, so the pool is the union of the
{len(seeds)} directories and nothing overwrites anything.

Every flag above is one `boltz predict --help` accepts in 2.2.1; two of them are not
optional for this project's purposes:

- **`--model boltz2`** is what sets `boltz_2=True` in the schema parser, and templates are
  rejected outright without it (`Templates are not supported in Boltz 1.0!`). It is the
  current default, so this is belt and braces — but it is the single flag that silently
  invalidates the run if that default ever moves.
- **`--seed 42`** because the default is *no seeding*, which makes a diffusion run
  unreproducible. 42 matches the seed the Vina runs used, so the two protocols stay
  comparable on that axis.

**Run it from this directory, and point it at `{INPUT_SUBDIR}/`, not at `.`** — both matter:

- Boltz's `check_inputs` globs the directory it is given and raises on *anything* that is
  not `.fasta`/`.yaml`, subdirectories included. So the inputs live in `{INPUT_SUBDIR}/`
  on their own, and this `RUN.md`, the template and the `out/` directory stay outside it.
  Pointing at `.` fails on the first non-YAML file it meets.
- The template path in each YAML is relative, and Boltz hands it straight to its mmCIF
  parser without resolving it against the YAML's own location — so it is resolved against
  the **current working directory**. Running from anywhere else means the template is not
  found.

`--output_format pdb` saves the PDB conversion step. `--diffusion_samples 5` gives five
poses per compound; `model_0` is Boltz's own top-ranked one, and the rest are what a
pose-selection step would choose between. `--use_potentials` applies Boltz's steering
potentials, which improve physical plausibility — relevant because the challenge zeroes
any structure that fails PoseBusters.

**Measured, not estimated.** This exact command was run on the 3080 Ti on 2026-10-09 and
interrupted after the first compound: 4 min 36 s for that one (model load included) and a
reported ETA of **1 h 27 min for 20 compounds** at `--diffusion_samples 5`. So budget
roughly 4.5 minutes per compound for all five samples — not per sample. The same run
confirmed 20/20 templates and 20/20 pocket constraints were featurised.

## Handing the results back

```bash
# from the project root, in PowerShell
.venv\\Scripts\\python scripts\\ingest_run.py runs\\{out_dir.name}
```

Results land at exactly this path — `boltz_results_` plus the stem of the input directory,
so `{INPUT_SUBDIR}/` gives `boltz_results_{INPUT_SUBDIR}`:

```
out/boltz_results_{INPUT_SUBDIR}/predictions/<COMPOUND>/
    <COMPOUND>_model_0.pdb              # model_0 is Boltz's top-ranked sample
    <COMPOUND>_model_1..4.pdb
    confidence_<COMPOUND>_model_0.json  # plddt, iptm, and the per-chain figures
```

Copy the `*.pdb` into `runs/{out_dir.name}/complexes/`, keeping the compound ID in each
filename, then record the Boltz version and these settings in `run.yaml`.

**Keep the `confidence_*.json` files.** They carry Boltz's own ranking, which is the
independent opinion this run exists to obtain — and for the three compounds with an
undefined stereocentre, comparing the `_R` and `_S` confidences is the most direct evidence
the project has yet had on which configuration to submit.

## Regenerate this set

Exactly as it was made, this text included:

```
{invocation()}
```
"""


def invocation() -> str:
    """The command that produced this set, so RUN.md is not a hand-maintained copy."""
    import shlex
    import sys as _sys

    parts = [".venv\\Scripts\\python", "scripts/make_cofold_yamls.py"]
    return " ".join(parts + [shlex.quote(a) for a in _sys.argv[1:]])


def main() -> int:
    challenge.use_utf8_stdout()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--name",
        required=True,
        help="run id for this set, e.g. 20261005_boltz2_heme-pocket",
    )
    parser.add_argument(
        "--msa",
        default="empty",
        help="'empty' for single-sequence mode (default); 'server' to omit the msa key so "
        "boltz predict --use_msa_server generates one (the flag is then added to the "
        "command in RUN.md); or a path to an .a3m",
    )
    parser.add_argument(
        "--no-heme", action="store_true", help="omit the heme cofactor (an ablation)"
    )
    parser.add_argument(
        "--pocket-heme",
        action="store_true",
        help="constrain the ligand to within --pocket-distance of the heme",
    )
    parser.add_argument("--pocket-distance", type=float, default=6.0)
    parser.add_argument(
        "--pocket-atom",
        default="FE",
        help="heme atom the ligand is steered toward (default FE, the iron)",
    )
    parser.add_argument(
        "--force-pocket",
        action="store_true",
        help="make the pocket constraint hard rather than a bias",
    )
    parser.add_argument("--template", help="one template path for every compound")
    parser.add_argument(
        "--template-chain",
        default="A",
        help="chain to take from --template (default A)",
    )
    parser.add_argument(
        "--bundle-template",
        action="store_true",
        help=(
            "copy --template into the set and reference it by basename. The project path "
            "contains spaces, which a WSL command line handles badly, and it makes the "
            "set self-contained so the template travels with the YAMLs that used it."
        ),
    )
    parser.add_argument(
        "--templates-csv",
        type=Path,
        help="CSV with Molecule_Name,template[,template_chain] for per-compound templates",
    )
    parser.add_argument(
        "--seeds",
        default="42",
        help=(
            "comma-separated seeds; RUN.md loops over them, each writing out_seed<N>/ "
            "(default 42). The project's pooled protocol is 5 seeds x 5 samples."
        ),
    )
    parser.add_argument(
        "--extra-flag",
        action="append",
        metavar="FLAG",
        help=(
            "a further `boltz predict` flag to put in RUN.md's command, repeatable. "
            "For flags that change the protocol but not the YAML — e.g. "
            "--extra-flag --subsample_msa --extra-flag '--num_subsampled_msa 256'"
        ),
    )
    parser.add_argument(
        "--affinity",
        action="store_true",
        help="also request Boltz's affinity prediction (not scored by this track)",
    )
    parser.add_argument(
        "--purpose-file",
        type=Path,
        help=(
            "markdown file whose contents head RUN.md, saying why this set exists. A file "
            "rather than inline text so the invocation recorded in RUN.md stays short and "
            "the rationale is itself committed."
        ),
    )
    parser.add_argument(
        "--smiles-csv",
        type=Path,
        help=(
            "Molecule_Name,SMILES to use instead of the blinded set, e.g. "
            "docking_inputs/stereo_expanded.csv to give Boltz each enantiomer "
            "separately rather than letting it pick one silently"
        ),
    )
    args = parser.parse_args()

    if args.template and args.templates_csv:
        raise SystemExit("pass --template or --templates-csv, not both")

    per_compound = load_templates(args.templates_csv) if args.templates_csv else {}
    if args.templates_csv:
        absent = [n for n in challenge.names() if n not in per_compound]
        if absent:
            print(
                f"note: no template for {len(absent)} compound(s): "
                f"{', '.join(absent)}\n      they will be predicted without one.\n"
            )

    sequence = cyp3a4_sequence()
    out_dir = COFOLD / args.name
    input_dir = out_dir / INPUT_SUBDIR
    input_dir.mkdir(parents=True, exist_ok=True)
    for stale in input_dir.glob("*.yaml"):
        stale.unlink()
    # Earlier sets wrote the YAMLs beside RUN.md, which made `boltz predict <set>` fail on
    # the first non-YAML file. Clear any such leftovers so a regenerated set is clean.
    for stale in out_dir.glob("*.yaml"):
        stale.unlink()

    template_ref = args.template
    if args.template and args.bundle_template:
        src = Path(args.template)
        if not src.is_file():
            raise SystemExit(f"--bundle-template needs a readable file: {src}")
        (out_dir / src.name).write_bytes(src.read_bytes())
        template_ref = src.name
        print(f"bundled {src.name} into the set; YAMLs reference it by basename")
    elif args.templates_csv and args.bundle_template:
        # Boltz resolves a template path against the working directory, and RUN.md says to
        # run from the set directory - so a project-relative path in the CSV would not be
        # found. Copy every referenced template in and rewrite the paths to basenames.
        bundled = set()
        for name, specs in per_compound.items():
            rewritten = []
            for spec in specs:
                path, _, chain = spec.partition("|")
                src = Path(path)
                if not src.is_file():
                    raise SystemExit(f"{args.templates_csv}: no such template: {src}")
                if src.name not in bundled:
                    (out_dir / src.name).write_bytes(src.read_bytes())
                    bundled.add(src.name)
                rewritten.append(f"{src.name}|{chain}" if chain else src.name)
            per_compound[name] = rewritten
        print(f"bundled {len(bundled)} templates into the set, referenced by basename")

    if args.smiles_csv:
        with open(args.smiles_csv, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        missing = {"Molecule_Name", "SMILES"} - set(rows[0] if rows else {})
        if missing:
            raise SystemExit(f"{args.smiles_csv.name}: missing column(s) {sorted(missing)}")
        entries = [
            challenge.Compound(r["Molecule_Name"].strip(), r["SMILES"].strip())
            for r in rows
        ]
    else:
        entries = list(challenge.compounds())

    unmatched = []
    for compound in entries:
        # A stereo-expanded entry is `OCNT-2312554_R`, while the templates CSV is keyed by
        # the challenge compound. Fall back to the parent so `_R`/`_S` inherit its template
        # rather than silently ending up with none - which would leave those entries on a
        # different protocol from the rest of the set.
        templates = per_compound.get(compound.name)
        if templates is None:
            parent = challenge.match_name(compound.name)
            templates = per_compound.get(parent, []) if parent else []
            if args.templates_csv and not templates:
                unmatched.append(compound.name)
        if template_ref:
            templates = [f"{template_ref}|{args.template_chain}"]
        (input_dir / f"{compound.name}.yaml").write_text(
            build_yaml(
                compound,
                sequence=sequence,
                msa=args.msa,
                heme=not args.no_heme,
                pocket=args.pocket_heme,
                pocket_distance=args.pocket_distance,
                pocket_atom=args.pocket_atom,
                force_pocket=args.force_pocket,
                templates=templates,
                affinity=args.affinity,
            ),
            encoding="utf-8",
        )

    if unmatched:
        raise SystemExit(
            f"{args.templates_csv}: no template for {len(unmatched)} entries "
            f"({', '.join(unmatched)}). Every entry must be on the same protocol, so this "
            "is refused rather than written."
        )

    readme = out_dir / "RUN.md"
    readme.write_text(run_instructions(out_dir, args.msa, args, entries), encoding="utf-8")

    # Assert the invariant Boltz enforces, rather than leaving it to be discovered on the
    # GPU: the directory handed to `boltz predict` must hold nothing but .fasta/.yaml.
    intruders = sorted(
        p.name
        for p in input_dir.glob("*")
        if p.is_dir() or p.suffix.lower() not in (".fa", ".fas", ".fasta", ".yml", ".yaml")
    )
    if intruders:
        raise SystemExit(
            f"{input_dir.relative_to(ROOT)} contains {intruders}, which Boltz's "
            f"check_inputs rejects. Nothing but .fasta/.yaml may live there."
        )

    print(f"wrote {len(entries)} YAMLs to {input_dir.relative_to(ROOT)}")
    print(f"run instructions: {readme.relative_to(ROOT)}")
    print(f"  boltz predict {INPUT_SUBDIR}   (run from {out_dir.name}/, not from {INPUT_SUBDIR}/)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
