"""Convert a Boltz-2 prediction directory into project complexes, and collect its scores.

Two things have to be repaired on the way in.

**Residue names.** Boltz's PDB writer hardcodes `"LIG"` for every HETATM record
(`boltz/data/write/pdb.py:100`), discarding the component name, so the heme and the ligand
arrive with the *same* residue name. The challenge requires the ligand to be `LIG` and
reads the heme as context, and this project has already lost a leaderboard result to the
scorer failing to match a mis-named heme. The chain assignment is recoverable because it is
ours: the YAMLs declare protein `A`, heme `B`, ligand `C`. That is asserted rather than
assumed — chain B must carry 43 atoms including an iron, or the file is refused.

(`--output_format mmcif` would have avoided this: the mmCIF writer keeps the component
name. Worth preferring next time; it does not justify re-running 1 h 45 min of GPU.)

**Frame.** Boltz predicts in its own frame, unrelated to the docking receptor's. Nothing is
done about that here — `scripts/compare_runs_poses.py` superposes on the heme, which is
exactly the invariant that makes the two comparable.

Also writes `boltz_confidence.csv`: the per-model scores from the `confidence_*.json`
files, including the ligand-heme interface confidence, which is Boltz's own opinion on the
coordination this project has been trying to get independent evidence for.

Usage (from the project root):
    .venv\\Scripts\\python scripts\\ingest_boltz_poses.py ^
        "cofold\\<set>\\out\\boltz_results_yaml\\predictions" ^
        --out runs\\20261009_boltz_9bvb-template-stereo
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import challenge
import pdbtools

HEME_ATOMS = 43
PROTEIN_CHAIN, HEME_CHAIN, LIGAND_CHAIN = "A", "B", "C"

# Chain order in the YAML, which is how Boltz indexes chains in pair_chains_iptm.
CHAIN_INDEX = {PROTEIN_CHAIN: "0", HEME_CHAIN: "1", LIGAND_CHAIN: "2"}


def split_chains(path: Path) -> dict[str, list[pdbtools.Atom]]:
    by_chain: dict[str, list[pdbtools.Atom]] = {}
    for raw in pdbtools.pdb_text(path).splitlines():
        if raw[:6] not in ("ATOM  ", "HETATM") or len(raw) < 54:
            continue
        atom = pdbtools.Atom(raw.ljust(80))
        by_chain.setdefault(atom.chain, []).append(atom)
    return by_chain


def confidence_row(path: Path) -> dict[str, float | str]:
    """Flatten a confidence JSON, keeping the interfaces that bear on the pose."""
    d = json.loads(path.read_text(encoding="utf-8"))
    pair = d.get("pair_chains_iptm") or {}
    lig, heme, prot = (
        CHAIN_INDEX[LIGAND_CHAIN],
        CHAIN_INDEX[HEME_CHAIN],
        CHAIN_INDEX[PROTEIN_CHAIN],
    )

    def cell(a: str, b: str):
        try:
            return float(pair[a][b])
        except (KeyError, TypeError, ValueError):
            return ""

    return {
        "confidence_score": d.get("confidence_score", ""),
        "ligand_iptm": d.get("ligand_iptm", ""),
        "complex_plddt": d.get("complex_plddt", ""),
        "complex_iplddt": d.get("complex_iplddt", ""),
        # The two that matter for a pose: how sure Boltz is about where the ligand sits
        # relative to the heme, and relative to the protein.
        "iptm_ligand_heme": cell(lig, heme),
        "iptm_ligand_protein": cell(lig, prot),
    }


def main() -> int:
    challenge.use_utf8_stdout()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("predictions", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--tag",
        default="",
        help=(
            "short label distinguishing this run's samples from others staged alongside "
            "them, e.g. s21 for a seed-21 run. Pose names become <entry>_<tag>m<rank> and "
            "the confidence CSV is suffixed to match, so several runs can share one "
            "run folder and be pooled."
        ),
    )
    args = parser.parse_args()

    complexes = args.out / "complexes"
    complexes.mkdir(parents=True, exist_ok=True)
    smiles = challenge.smiles_map()

    rows, written, refused = [], 0, []
    for entry in sorted(p for p in args.predictions.iterdir() if p.is_dir()):
        parent = challenge.match_name(entry.name)
        if not parent:
            refused.append(f"{entry.name}: no challenge compound in the name")
            continue
        for model in sorted(entry.glob("*_model_*.pdb")):
            rank = model.stem.rsplit("_model_", 1)[1]
            stem = f"{entry.name}_{args.tag}m{rank}" if args.tag else f"{entry.name}_m{rank}"
            chains = split_chains(model)

            heme = chains.get(HEME_CHAIN, [])
            ligand = chains.get(LIGAND_CHAIN, [])
            protein = chains.get(PROTEIN_CHAIN, [])
            heavy_heme = [a for a in heme if a.element != "H"]
            if len(heavy_heme) != HEME_ATOMS or not any(
                a.element == "FE" for a in heavy_heme
            ):
                refused.append(
                    f"{stem}: chain {HEME_CHAIN} has {len(heavy_heme)} heavy atoms and "
                    f"{'an' if any(a.element == 'FE' for a in heavy_heme) else 'no'} iron; "
                    f"expected {HEME_ATOMS} with an iron, so the chain assignment is not "
                    f"what the YAML declared"
                )
                continue
            if not ligand or not protein:
                refused.append(f"{stem}: missing chain B or C")
                continue

            pdbtools.write_complex(
                complexes / f"{stem}.pdb",
                protein=protein,
                # Boltz wrote the heme as LIG; the challenge reads it as HEM.
                heme=pdbtools.set_resname(heavy_heme, "HEM"),
                ligand=pdbtools.set_resname(
                    [a for a in ligand if a.element != "H"], "LIG"
                ),
                header=[f"REMARK   3 BOLTZ {entry.name} model {rank}"],
            )
            written += 1

            conf = entry / f"confidence_{entry.name}_model_{rank}.json"
            row = {
                "name": stem,
                "parent": parent,
                "variant": challenge.variant_tag(stem, parent),
                "model": rank,
                "pool_tag": args.tag,
                "n_ligand_atoms": len([a for a in ligand if a.element != "H"]),
                "expected_heavy_atoms": "",
            }
            if parent in smiles:
                from rdkit import Chem, RDLogger

                RDLogger.DisableLog("rdApp.*")
                mol = Chem.MolFromSmiles(smiles[parent])
                row["expected_heavy_atoms"] = mol.GetNumHeavyAtoms() if mol else ""
            row.update(confidence_row(conf) if conf.exists() else {})
            rows.append(row)

    if rows:
        name = f"boltz_confidence_{args.tag}.csv" if args.tag else "boltz_confidence.csv"
        out_csv = args.out / name
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {out_csv}")

    print(f"wrote {written} complexes to {complexes}")
    mismatched = [
        r
        for r in rows
        if r["expected_heavy_atoms"] not in ("", r["n_ligand_atoms"])
    ]
    print(
        f"ligand heavy-atom count matches the challenge SMILES: "
        f"{len(rows) - len(mismatched)}/{len(rows)}"
    )
    for m in mismatched[:5]:
        print(f"  MISMATCH {m['name']}: {m['n_ligand_atoms']} vs {m['expected_heavy_atoms']}")
    for r in refused:
        print(f"  refused {r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
