"""Pick each compound's pose as the centre of its largest binding mode.

**Why this rather than the medoid.** The medoid minimises mean RMSD to every other sample,
which is the right estimator for one diffuse cloud and the wrong one for a pool split across
modes: the most central point can land inside a minority cluster that happens to lie between
the others. `OCNT-2313300` is exactly that - its medoid sat in a cluster of 3 while the
plurality held 7 - and the alt-modes submission, which moved it to the plurality pose among two
other changes, gained +0.336 of summed LDDT-PLI.

This is a **general rule applied uniformly to all 20**, not a per-compound correction. It
differs from the medoid wherever a pool is multimodal and agrees with it everywhere else, so it
is one protocol variable and one slot.

**Ties are broken explicitly.** Where two modes hold the same number of samples - `OCNT-2395342`
splits 8 against 8 - "largest" is undefined and picking by cluster label would be arbitrary.
The tie goes to the **tighter** mode, the one whose members agree more closely with each other,
on the grounds that a mode the model keeps reproducing precisely is better evidence than one it
reproduces loosely. Stated here because an unstated tiebreak is how an arbitrary choice gets
mistaken for a rule.

Modes are average-linkage clusters of symmetry-corrected, heme-superposed RMSD cut at --cut A,
the same 2 A the project uses for "different binding mode" everywhere else.

Usage (from the project root):
    .venv\\Scripts\\python scripts\\plurality_selection.py runs\\<run-id>
    .venv\\Scripts\\python scripts\\plurality_selection.py runs\\<run-id> --csv picks.csv
"""

from __future__ import annotations

import argparse
import csv
import itertools
from pathlib import Path

import challenge
import numpy as np
import pose_set_diversity as psd

STEREO = frozenset({"OCNT-2312554", "OCNT-2395415", "OCNT-2312966"})


def main() -> int:
    challenge.use_utf8_stdout()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--config", default="R")
    parser.add_argument("--cut", type=float, default=2.0)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()

    from rdkit import RDLogger
    from rdkit.Chem import rdMolAlign
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    RDLogger.DisableLog("rdApp.*")
    smiles = challenge.smiles_map()

    groups: dict[str, list[Path]] = {}
    for path in sorted((args.run / "complexes").glob("*.pdb")):
        parent = challenge.match_name(path.stem)
        if not parent:
            continue
        tag = challenge.variant_tag(path.stem, parent)
        if parent in STEREO and not tag.startswith(args.config.upper()):
            continue
        groups.setdefault(parent, []).append(path)

    print(f"{'compound':<15}{'modes':>6}{'sizes':>16}{'tie?':>6}"
          f"{'same as medoid?':>17}")
    rows, changed = [], 0
    for parent in challenge.names():
        paths = groups.get(parent, [])
        if len(paths) < 5:
            continue
        ref = paths[0]
        mols = {}
        for p in paths:
            m = psd.ligand_mol(p, smiles[parent], frame=ref)
            if m is not None:
                mols[p] = m
        keys = list(mols)
        n = len(keys)
        D = np.zeros((n, n))
        for i, j in itertools.combinations(range(n), 2):
            D[i, j] = D[j, i] = rdMolAlign.CalcRMS(mols[keys[i]], mols[keys[j]])

        labels = fcluster(
            linkage(squareform(D, checks=False), method="average"),
            t=args.cut, criterion="distance",
        )
        medoid = keys[int(np.argmin(D.mean(axis=1)))]

        # candidate modes, largest first; tie broken by tightness (mean within-mode RMSD)
        stats = []
        for lab in sorted(set(labels)):
            members = [i for i in range(n) if labels[i] == lab]
            sub = D[np.ix_(members, members)]
            tight = float(sub[np.triu_indices(len(members), 1)].mean()) if len(members) > 1 else 0.0
            centre = keys[members[int(np.argmin(sub.mean(axis=1)))]]
            stats.append((len(members), -tight, lab, centre))
        stats.sort(reverse=True)
        best_size = stats[0][0]
        tied = sum(1 for s in stats if s[0] == best_size) > 1
        pick = stats[0][3]

        same = pick is medoid
        changed += not same
        sizes = " ".join(str(s[0]) for s in stats)
        print(f"{parent:<15}{len(stats):>6}{sizes:>16}{('yes' if tied else ''):>6}"
              f"{('yes' if same else 'NO'):>17}")
        rows.append(
            {
                "compound": parent,
                "pick": pick.stem,
                "variant": challenge.variant_tag(pick.stem, parent),
                "mode_size": best_size,
                "n_modes": len(stats),
                "tie_broken_on_tightness": str(tied),
                "medoid": medoid.stem,
                "differs_from_medoid": str(not same),
            }
        )

    print(f"\n{len(rows)} compounds | differs from the medoid in {changed}")
    if args.csv and rows:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.csv}")
    print("\nPicks for build_submission.py:\n")
    print(" \\\n".join(f"  --pick {r['compound']}={r['variant']}" for r in rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
