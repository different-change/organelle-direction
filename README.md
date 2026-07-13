# Organellome Signature — a directional readout of cell biology

**Read the *direction* of a cell's organelle machinery — is it building or clearing — and, from gene
expression alone, you can predict which cancers are addicted to their mitochondria and which respond to
the drugs that target them. Direction beats amount, and no single molecular layer sees it.**

Built end-to-end with Claude Science, on public data only.

---

## The one idea

Two cells can hold the *same amount* of an organelle for opposite reasons: one is **building** it, the other
has **blocked its clearance**. To the cell those are opposite states; to a standard "how much is there"
readout they look identical — and the genes that would tell them apart (selective-autophagy receptors) are
not even in a conventional organelle signature.

We score the two opposing programs separately and subtract them:

```
net direction  =  biogenesis  −  selective degradation
```

One signed number per sample, read from ordinary transcriptomes, proteomes, or phosphoproteomes. Positive =
net building; negative = net clearing. It is a **set-point** (which program dominates now), validated against
timecourses where the true direction is known — not an organelle headcount, and not a measured rate.

---

## What the method delivers (the well-powered spine)

| Result | Number | Data |
|---|---|---|
| **Direction → genetic vulnerability** | ρ = −0.35, p = 1.7×10⁻³², n = 1,066 cell lines | DepMap CRISPR |
| ↳ survives dependency-burden control | partial r = −0.32 | |
| ↳ lineage-robust | negative in 17/18 lineages | |
| ↳ mitochondria-specific | mito→mito −0.35 vs mito→ribosome −0.07 | |
| **Direction → drug response** | IACS-010759 at 3rd percentile of 1,514 drugs; MitoQ #1 | DepMap PRISM |
| ↳ survives a proliferation confound | mito-drug class MWU p = 0.017 | |
| **Direction beats amount** | AUROC 0.666 (direction) vs 0.613 (raw expression) | on the vulnerability endpoint |
| **Generalizes, graded by organelle** | mito −0.35 / ER −0.15 (survives controls) / lysosome null | DepMap |
| **Cross-modal agreement in human tumor** | mito-biogenesis RNA g=−1.52 · protein g=−3.14 | CPTAC ccRCC |

Effect sizes are modest (ρ 0.15–0.35; AUROC 0.60–0.67) — as expected for a single-pathway expression score
predicting a functional phenotype across highly heterogeneous cell lines. The value is a **real, specific,
lineage-robust, correctly-signed** signal, not a high-accuracy point predictor.

---

## Factchecks and provenance

Every headline here survived an adversarial check, and we report the ones that didn't:

- **Reported a null at power.** A weak survival hint in CPTAC kidney cancer (Cox p=0.04, 21 deaths) did **not**
  replicate when retested at 8× the power in TCGA-KIRC (508 tumors, 168 deaths; Cox p=0.21, log-rank p=0.65).
  We show the flat survival curves — positive controls (stage HR=1.93, p=5.6×10⁻²¹) confirm the pipeline works.
- **Demoted its own statistic.** Swapped an over-generous equal-variance test (p=0.0005) for the correct
  Welch test (p=0.0094) — chose the weaker, right number.
- **Refuted its own headline.** Flagged a raw r=0.88 "coupling" as a shared-timecourse artifact and re-led with
  the mechanism (TFEB-dependence) instead.
- **Walked back a replication claim.** A second drug platform (GDSC) looked like independent replication;
  the confound test showed it was mostly a general drug-sensitivity axis, so we downgraded it to
  "weak directional consistency." PRISM survives the same test; GDSC does not.
- **Pre-registered its predictions.** Nine blind directional predictions were committed and SHA-256-hashed
  *before* looking at the data (see `data/PREREGISTRATION.md`). 4/9 hit — including a high-confidence miss
  (COAD) locked in before data, which is the proof the hits aren't overfit.

---

## Repository layout

```
report/
  report.html                  ← the judging document: 7-beat spine, standalone-readable figures
  full_technical_report.html   ← complete technical record (35 figures, all analyses)
demo/
  index.html                   ← interactive walkthrough (open in a browser; press P for Presenter Mode)
figures/                       ← the six 300-DPI figures, each readable cold
scorer/
  run_dynamics.py              ← run the directional score on your own expression matrix
  score_organelle_dynamics.py
  README.md                    ← scorer usage
data/
  genesets/                    ← curated biogenesis / selective-clearance modules (5 organelles × 4 organisms)
  scores/                      ← result tables behind every figure
  prospective_predictions.*    ← the hash-locked blind predictions
  PREREGISTRATION.md           ← human-readable pre-registration record
```

## Run the scorer

```bash
cd scorer
python run_dynamics.py --expression your_matrix.csv --organelle mitochondria --organism human
```

The score is `rank_pct(biogenesis genes).mean − rank_pct(selective-clearance genes).mean` per sample.
Gene modules and their provenance are in `data/genesets/`.

## Data sources (all public)

DepMap 24Q2 (expression, CRISPR gene-effect, PRISM drug response) · TCGA-KIRC (cBioPortal) ·
CPTAC ccRCC (RNA, protein, phosphoproteome) · MoTrPAC · GEO series (accessions in the technical report).

## How Claude Science was used

The entire pipeline — gene-set curation, scoring, every confound control, all figures, both reports, and the
adversarial self-refereeing above — was built in Claude Science with a scientist in the loop.
