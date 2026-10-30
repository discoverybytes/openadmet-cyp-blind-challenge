# execution

The code that produced and checked the submitted archive. This folder is its own project
root, so the scripts resolve `data/` and `vendor/` relative to it and run from a clone
without edits.

```
pip install rdkit numpy pandas MDAnalysis posebusters gemmi    # python 3.12
# boltz == 2.2.1 goes in a separate env with CUDA
```

## The pipeline, in order

| Step | Script |
|---|---|
| Write one Boltz-2 YAML per entry, plus the run command | `scripts/make_cofold_yamls.py` |
| Stage the predictions, repairing the haem/ligand residue names | `scripts/ingest_boltz_poses.py` |
| Cluster the 25 poses and pick the biggest cluster's centre | `scripts/plurality_selection.py` |
| Normalise, zip, and check the archive | `scripts/build_submission.py` |

`gate.py` is the check `build_submission.py` runs: the official `structure_validation.py`
with the expected IDs and SMILES, then PoseBusters in `dock` configuration with protein
**and** haem as context, which is how the scorer calls it. An archive is "ready" only when
both pass with zero failures. `challenge.py`, `pdbtools.py`, `ingest_run.py` and
`pose_set_diversity.py` are the shared helpers those four import.

Every command is in the report, one section up: [../README.md](../README.md).

## vendor/

`vendor/openadmet_cyp_tutorial/` is the official validator and the PoseBusters call pattern,
vendored unmodified from [OpenADMET/CYP-Challenge-Tutorial](https://github.com/OpenADMET/CYP-Challenge-Tutorial)
at commit `832ae6d`. Apache-2.0; its `LICENSE` travels with it. Vendored rather than
reimplemented so the local check is the same code the scorer runs.

## What is not here

This is a curated copy of a private working repository: the lab notebook, the run manifests,
the project's own findings documents and the rest of its analysis scripts are not published.
The four steps above are the whole protocol — nothing withheld is needed to reproduce it.
