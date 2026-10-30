# Vendored: OpenADMET CYP-Challenge-Tutorial

Source: <https://github.com/OpenADMET/CYP-Challenge-Tutorial>, commit
`832ae6d6447b07c2c71059de34d07f1130388247` (2026-09-25). Licence: Apache-2.0 (see `LICENSE`).
These files are unmodified copies. Do not edit them. To refresh them, re-download at a new
commit and update this file and decision 0003.

| File | Upstream path | Used for |
|---|---|---|
| `structure_validation.py` | `validation/structure_validation.py` | Gate step 1: the official zip validator |
| `posebusters_utils.py` | `evaluation/posebusters_utils.py` | Gate step 2: building the ligand exactly as the scorer does before PoseBusters runs |
| `CYP3A4_protein_sequence.fasta` | `inputs/CYP3A4_protein_sequence.fasta` | Boltz-2 YAML sequence (UniProt P08684) |
| `reference/evaluate_predictions.py` | `evaluation/evaluate_predictions.py` | Reference only (not importable, see below). `scripts/gate.py` mirrors `score_single_structure`'s PoseBusters block (upstream lines ~723–776) |
| `reference/config.py` | `evaluation/config.py` | Reference only |

`reference/evaluate_predictions.py` imports `POSEBUSTERS_MAX_FAILURES` and
`BISYRMSD_NAN_PENALTY` from its `config.py`, but the public `config.py` doesn't define them.
The public evaluation code therefore lags the challenge backend (brief: Risk 8), and the real
tolerance is unknown.
