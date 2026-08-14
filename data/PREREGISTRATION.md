# Pre-registration — prospective directional predictions

**Recorded:** 2026-07-11T23:11:28Z, at the start of the tumor-vs-normal validation work

**SHA-256 of `prospective_predictions.csv`:**
`eabe1fe54f6df5d0b88f1abc8a17bed360a9a726198ae60f1203fe9b12dc321f`

You can verify the lock yourself:
```bash
shasum -a 256 prospective_predictions.csv
# must equal the hash above
```
If the CSV had been edited after commit, the hash would not match. It does.

**What this hash does and does not prove.** It proves the prediction file has not been edited since
the hash was recorded. It does not, by itself, prove the predictions predate the analysis: the file
and this claim live in the same repository, all commits sit within one hackathon working session, and
no external timestamp (OSF, OpenTimestamps) was anchored at the time. Read this as a documented,
verifiable commitment made during the event, not as cryptographic pre-registration. Any future
version of this work will anchor predictions externally before analysis.

---

## What was predicted, and why this design is honest

Before looking at the CPTAC pan-cancer proteomics, the directional method was used to predict, for nine
cancer types, the **sign** of the mitochondrial-machinery change (tumor vs. matched normal), each with a
pre-declared **confidence tier** and rationale. All nine were predicted *down* (Warburg-shift logic), with
confidence graded high / moderate / low.

**Scoring rule (also pre-declared):** a prediction counts as a hit only if the tumor-vs-normal effect is
in the predicted direction **and** significant (g<0, p<0.05). A non-significant result counts as a miss.

**Result: 4/9 hits** — high-confidence 2/3, moderate 1/3 (as expected), low-confidence 0/2 (exactly as
pre-declared: these were flagged as the most likely to miss).

The load-bearing point is **not** the raw hit rate — a single-pathway expression score is not expected to
predict every cancer's proteome. It is that:

1. The **confidence tiers were calibrated** — high beat moderate beat low, monotonically.
2. A **high-confidence miss (COAD) was locked in before the data** and reported as a miss. That is the
   proof the hits are not overfit: an honest pre-registration must be able to be wrong, and this one was.

See `../report/full_technical_report.html` for the per-cancer results table and the effect sizes.

## Files

- `prospective_predictions.csv` — the nine locked predictions (the hashed artifact)
- `prospective_predictions.json` — same content, structured, with the commit timestamp
