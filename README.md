# OpenADMET Blind Challenge: Predicting CYP Inhibition

Entry by **discoverybytes** (John Proudfoot).

The challenge ran three tracks. Only the structure track was entered.

| Track | Entered | Result |
|---|---|---|
| [Structure Prediction](structure-prediction/) | yes | LDDT-PLI **0.4548** |
| Direct Inhibition (regression) | no | — |
| Time-Dependent Inhibition (classification) | no | — |

**[structure-prediction/](structure-prediction/)** holds the method report, the submitted
archive with the checks it passed, and the code that built it.

The task: predict 20 CYP3A4–ligand complexes from ligand SMILES alone, scored against
cryo-EM structures that stayed blinded until after the deadline.

Challenge: <https://huggingface.co/spaces/openadmet/cyp-challenge>
