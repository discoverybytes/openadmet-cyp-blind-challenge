"""How much pose choice a run actually offers, measured against what it submitted.

Re-selecting from an all-poses export can only help if the alternatives differ from the
pose already chosen. This measures that directly: for each compound, the symmetry-corrected
RMSD of every alternative pose against the one the submitted run exported, plus the widest
separation within the compound's own pose set.

If the alternatives all sit within the noise of the submitted pose, the export contains no
better submission whatever it is scored on, and the limit is sampling rather than
selection. If some differ by several angstrom, there is a real decision to argue.

The submitted pose is identified by matching coordinates against the submitted run's own
exported complex, not by trusting a filename convention: the two compounds whose poses were
corrected late have no suffix-free entry in the export, so a name-based guess would silently
compare the wrong things.

Usage (from the project root):
    .venv\\Scripts\\python scripts\\pose_set_diversity.py ^
        runs\\20261008_tab3_apf-9bvb-allposes --submitted runs\\20261007_tab3_apf-9bvb-v3
"""

from __future__ import annotations

import argparse
import statistics as st
from pathlib import Path

import challenge


def ligand_mol(path: Path, smiles: str, frame: Path | None = None):
    """The pose's ligand with bond orders from the challenge SMILES, hydrogens dropped.

    With `frame`, the ligand is first moved into that complex's frame by superposing the
    two hemes. This is **not** optional for every run: a docking run writes every pose
    against one receptor, but Boltz emits each prediction in its own arbitrary frame, so
    comparing two of its models directly measures the frame offset rather than the pose.
    Superposing is a no-op when the frames already agree, so it is always applied.
    """
    import sys

    vendor = Path(__file__).resolve().parent.parent / "vendor"
    if str(vendor) not in sys.path:
        sys.path.insert(0, str(vendor))
    from openadmet_cyp_tutorial.posebusters_utils import bonded_ligand_copies
    from rdkit import Chem

    copies = bonded_ligand_copies(str(path), smiles)
    if not copies:
        return None
    mol = Chem.RemoveHs(copies[0])
    if frame is None or frame == path:
        return mol

    import numpy as np

    import compare_runs_poses as crp

    here, there = crp.heme_frame(path), crp.heme_frame(frame)
    shared = sorted(set(here) & set(there))
    if len(shared) < 20:
        raise SystemExit(
            f"{path.name} and {frame.name} share only {len(shared)} named heme atoms, "
            f"so they cannot be brought into a common frame"
        )
    rotation, pc, qc, _ = crp.kabsch([here[k] for k in shared], [there[k] for k in shared])
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = np.array(list(conf.GetAtomPosition(i)))
        conf.SetAtomPosition(i, (((p - pc) @ rotation) + qc).tolist())
    return mol


def main() -> int:
    challenge.use_utf8_stdout()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--submitted", type=Path, required=True)
    args = parser.parse_args()

    from rdkit import RDLogger
    from rdkit.Chem import rdMolAlign

    RDLogger.DisableLog("rdApp.*")
    smiles = challenge.smiles_map()

    poses: dict[str, list[Path]] = {}
    for path in sorted((args.run / "complexes").glob("*.pdb")):
        parent = challenge.match_name(path.stem)
        if parent:
            poses.setdefault(parent, []).append(path)

    print(
        f"{'compound':<16}{'n':>3}{'submitted':>18}{'nearest alt':>12}"
        f"{'furthest alt':>14}{'set width':>11}"
    )
    widths, furthest_all = [], []
    for parent in challenge.names():
        members = poses.get(parent, [])
        if not members:
            continue
        submitted_path = args.submitted / "complexes"
        # Find which export matches the submitted complex, by coordinates.
        target = None
        for cand in sorted(submitted_path.glob(f"{parent}*.pdb")):
            ref = ligand_mol(cand, smiles[parent])
            if ref is None:
                continue
            for m in members:
                mol = ligand_mol(m, smiles[parent])
                if mol is None or mol.GetNumAtoms() != ref.GetNumAtoms():
                    continue
                if rdMolAlign.CalcRMS(mol, ref) < 0.05:
                    target = (m, ref)
                    break
            if target:
                break
        if target is None:
            print(f"{parent:<16}{len(members):>3}   could not match the submitted pose")
            continue

        chosen, ref = target
        mols = {}
        for m in members:
            mol = ligand_mol(m, smiles[parent])
            if mol is not None and mol.GetNumAtoms() == ref.GetNumAtoms():
                mols[m] = mol
        alts = {k: v for k, v in mols.items() if k != chosen}
        if not alts:
            print(f"{parent:<16}{len(members):>3}{chosen.stem:>18}  (no alternatives)")
            continue
        d = {k: rdMolAlign.CalcRMS(v, ref) for k, v in alts.items()}
        near = min(d.values())
        far = max(d.values())
        # Widest separation between any two poses in the set, submitted one included.
        keys = list(mols)
        width = max(
            rdMolAlign.CalcRMS(mols[keys[i]], mols[keys[j]])
            for i in range(len(keys))
            for j in range(i + 1, len(keys))
        )
        widths.append(width)
        furthest_all.append(far)
        print(
            f"{parent:<16}{len(members):>3}{chosen.stem:>18}"
            f"{near:>11.2f} Å{far:>13.2f} Å{width:>10.2f} Å"
        )

    print(
        f"\nacross {len(widths)} compounds: median set width {st.median(widths):.2f} Å, "
        f"max {max(widths):.2f} Å"
        f"\nfurthest alternative from the submitted pose: median "
        f"{st.median(furthest_all):.2f} Å, max {max(furthest_all):.2f} Å"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
