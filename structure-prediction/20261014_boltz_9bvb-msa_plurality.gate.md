# Gate report — `20261014_boltz_9bvb-msa_plurality.zip`

- SHA-256: `cc1adc45f513681c8c546967123c12de23a044d93dd21fb8346cc2b7b6606734`
- Validator: PASS
- PoseBusters: PASS (0 of 20 compounds failing)
- **Verdict: READY**

## PoseBusters per compound

Relative distances are (closest contact) / (sum of radii); below **0.75** that check fails. `Fe` uses covalent radii, which is why real coordination passes; `heme` uses van der Waals radii and is the check that off-axis coordination trips. See `docs/findings/posebusters-crystal-reference.md`.

A `Fe···lig` distance near 2.0-2.2 Å means the pose is coordinating the heme iron; 3 Å or more means it is not.

| Compound | Failed | Fe···lig (Å) | rel. protein | rel. heme | rel. Fe | energy ratio | Ligand bonds | Note |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | :--- | :--- |
| OCNT-2312554 | — | 2.83 | 0.78 | 0.92 | 1.39 | 2.5 | template | — |
| OCNT-2312966 | — | 3.00 | 0.86 | 0.98 | 1.44 | 3.8 | template | — |
| OCNT-2313038 | — | 2.83 | 0.89 | 0.94 | 1.39 | 3.1 | template | — |
| OCNT-2313300 | — | 2.84 | 0.97 | 0.90 | 1.40 | 1.8 | template | — |
| OCNT-2313397 | — | 2.83 | 0.82 | 0.95 | 1.40 | 1.3 | template | — |
| OCNT-2313429 | — | 3.43 | 0.78 | 1.06 | 1.65 | 1.5 | template | — |
| OCNT-2315379 | — | 3.49 | 0.82 | 0.99 | 1.68 | 1.2 | template | — |
| OCNT-2395328 | — | 3.93 | 0.78 | 1.00 | 1.89 | 2.2 | template | — |
| OCNT-2395342 | — | 3.37 | 0.79 | 1.07 | 1.62 | 2.0 | template | — |
| OCNT-2395415 | — | 2.85 | 0.85 | 0.90 | 1.40 | 3.8 | template | — |
| OCNT-2395436 | — | 2.87 | 0.94 | 0.91 | 1.42 | 2.3 | template | — |
| OCNT-2395457 | — | 2.83 | 0.92 | 0.95 | 1.40 | 2.2 | template | — |
| OCNT-2395460 | — | 4.04 | 0.91 | 1.00 | 1.94 | 2.0 | template | — |
| OCNT-2395464 | — | 2.85 | 0.85 | 0.92 | 1.40 | 1.4 | template | — |
| OCNT-2395503 | — | 2.84 | 0.97 | 0.93 | 1.40 | 2.6 | template | — |
| OCNT-2395645 | — | 3.64 | 0.81 | 0.97 | 1.75 | 1.5 | template | — |
| OCNT-2395712 | — | 3.48 | 0.78 | 1.12 | 1.67 | 1.1 | template | — |
| OCNT-2395789 | — | 2.87 | 0.85 | 0.90 | 1.41 | 2.1 | template | — |
| OCNT-2535422 | — | 3.47 | 0.91 | 0.95 | 1.67 | 1.2 | template | — |
| OCNT-2535519 | — | 2.83 | 0.87 | 0.89 | 1.40 | 1.4 | template | — |

### Passing, but close to a clash cutoff

These poses pass now; a small coordinate change would fail them.

- OCNT-2312554: tightest relative distance 0.779
- OCNT-2313429: tightest relative distance 0.780
- OCNT-2395328: tightest relative distance 0.781
- OCNT-2395342: tightest relative distance 0.788
- OCNT-2395712: tightest relative distance 0.785
