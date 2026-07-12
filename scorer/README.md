# Organelle Dynamics Kit — directional organelle-dynamics scoring

A small, reusable engine that scores **organelle build-up vs tear-down** from a
gene-expression (or protein-abundance) matrix. It separates two directional gene
modules — **biogenesis** (import / assembly / proliferation machinery) and
**degradation** (selective-autophagy machinery) — and returns a signed
**net_direction** per sample, so the *sign* distinguishes a shift toward organelle
proliferation from a shift toward clearance. A generic "organelle activity"
signature cannot do this; both directions move many of the same genes.

**Sixteen modules across seven organelles and four species ship in this kit** (the
engine is organelle-agnostic — add more by dropping in gene-set CSVs with the same
columns). Five organelles are now built out in **human** (mitochondrion, lysosome, ER,
lipid droplet, ferritin/iron), each validated scout-first against a public dataset with a
known-direction positive control.

| organelle | species | biogenesis | degradation | selective | clearance arm |
|---|---|---|---|---|---|
| peroxisome | yeast | 30 | 20 | 1 (Atg36) | **works** |
| peroxisome | *Arabidopsis* | 26 | 17 | 0 | empty (no receptor) |
| peroxisome | rice | 27 | 23 | 0 | empty (no receptor) |
| mitochondrion | yeast | 33 | 10 | 2 (Atg32/33) | **works** |
| mitochondrion | *Arabidopsis* | 27 | 14 | 0 | empty (no receptor) |
| mitochondrion | rice | 35 | 23 | 0 | empty (no receptor) |
| **mitochondrion** | **human** | **83** (machinery) | **20** | **6 (PINK1/PRKN/NIX/BNIP3/FUNDC1)** | **works** |
| **lysosome** | **human** | **79** (CLEAR) | **29** | **6 (LGALS3/8/9, TRIM16, UBE2QL1, FBXO27)** | **works** |
| **ER** | **human** | **73** (machinery) | **25** | **7 (FAM134B/TEX264/CCPG1/RTN3/ATL3/…)** | works (mRNA-underpowered) |
| **lipid droplet** | **human** | **31** (machinery) | **0** | **0** | empty (no LD receptor)† |
| **ferritin / iron** | **human** | **5** (FTH1/FTL/…) | **13** | **1 (NCOA4)** | partial (post-translational)‡ |
| chloroplast | *Arabidopsis* | 37 | 26 | 7 candidates | effectively empty* |
| chloroplast | rice | 31 | 17 | 0 | empty (no receptor) |
| ER | yeast | 29 | 10 | 2 (Atg39/40) | **works** |
| ER | *Arabidopsis* | 21 | 14 | 0 | empty (no receptor) |
| ER | rice | 32 | 23 | 0 | empty (no receptor) |

The clearance arm is only as good as its selective receptor. The pattern is the
platform's central lesson, and the full grid makes it unmistakable: **the
selective-degradation arm is populated wherever a genuine selective-autophagy receptor
exists** (yeast peroxisome Atg36, mitophagy Atg32/33, ER-phagy Atg39/40; human mitophagy
PINK1/Parkin/NIX/BNIP3/FUNDC1, lysophagy galectins, ER-phagy FAM134B family) **and empty
where there is none** — every plant cell (no receptor orthologs) and the human lipid
droplet (lipophagy is bulk autophagy, no LD-selective receptor). The pattern is biological,
not a scorer defect.

†**Lipid droplet is biogenesis-only** — lipophagy has no clean LD-selective receptor, so
its clearance file is empty by construction, exactly as for plant organelles.
‡**Ferritin's ferritinophagy arm** has one genuine selective receptor (NCOA4) but is *not
mRNA-validatable*: NCOA4 is controlled post-translationally, so its transcript is flat even
when ferritinophagy is induced. The iron-storage biogenesis arm is validated (FAC loading →
FTH1/FTL up); the clearance arm is reported honestly as a partial.

**Coupled-organelle fingerprints (the real novelty).** Beyond single organelles, the kit's
directional scores can be read *jointly* across coupled organelles. The **mito–lysosome
axis** (mitophagy delivers cargo to lysosomes) shows strong coupling under shared
lysosomal-damage stress (r = +0.88), with lysosome recovery TFEB-dependent — a QC signature
no single-organelle score captures. A **lipid–mito–ER metabolic-disease fingerprint** and a
**ferritin–mito ferroptosis axis** are likewise built from the validated arms.

**Within-system net-direction integration.** In hiPSC→cardiomyocyte maturation the mito
net_direction *changes sign within one system* (iPSC +0.82 build → contractile CM −0.50
clearance; trend p = 0.013), the selective-mitophagy arm overtaking biogenesis — the
integrated directional read the two-system validation had lacked.

**Human mitochondrion — machinery vs regulators.** The primary biogenesis score is
`mito_biogenesis_machinery_human.csv` (import/MICOS/OXPHOS-assembly/mtDNA apparatus).
The master transcriptional regulators (PGC1α/PPARGC1A, NRF1, TFAM, ESRRA, GABPA) live in
a *separate* file `mito_regulators_human.csv` and are **not** in the primary score — because
PGC1α induction is the textbook definition of the exercise response, scoring exercise with
PGC1α inside biogenesis would be near-tautological. A general-anabolism negative control
(`ribosome_biogenesis_anabolism_human.csv`) is included to test that the mito score is
organelle-specific, not a generic growth signal.

*The 7 chloroplast "selective" candidates (ATI1/2, CV, PUB4, SP1/2, SPL2) are carried
but were **not** coordinately senescence-induced (validation d = −0.26, n.s.) — plants
have no clean chloroplast-selective receptor, so treat this arm as empty in practice.

## What's here
- `score_organelle_dynamics.py` — the organelle-agnostic scoring engine
  (pure function; empirical permutation null; single-sample rank scoring; no
  hard-coded gene lists).
- `run_dynamics.py` — `load_module(species)` convenience loader + demo.
- `genesets/` — resolver-backed gene modules with provenance:
  - `peroxisome_biogenesis_yeast.csv` / `peroxisome_pexophagy_yeast.csv`
  - `peroxisome_biogenesis_arabidopsis.csv` / `peroxisome_pexophagy_arabidopsis.csv`
  - `geneset_provenance.json` — GO terms used, resolver services, collision
    flags, and genes dropped in cross-taxa transfer, with reasons.

## Quick start
```python
import pandas as pd
from score_organelle_dynamics import score_organelle_dynamics
from run_dynamics import load_module

bio, deg, deg_sel = load_module("yeast", "peroxisome")
# other combos: ("arabidopsis","peroxisome"), ("yeast","mitochondrion"),
#   ("arabidopsis","chloroplast"), ("rice","chloroplast"), ("yeast","er")
expr = pd.read_csv("my_matrix.csv", index_col=0)      # genes x samples, IDs in module namespace
scores = score_organelle_dynamics(
    expr, bio, deg, degradation_selective_ids=deg_sel,
    method="rank", n_perm=1000)
print(scores[["biogenesis_score","degradation_score","net_direction"]])
```
Gene ID namespace: **SGD systematic names** (e.g. `YDR329C`) for yeast, **AGI
loci** (e.g. `AT3G47430`) for Arabidopsis. Pass `id_map=` to translate a matrix
in another namespace.

## How it was validated (built-in ground truth)
- **Yeast, oleate vs glucose (GSE5862):** oleate proliferates peroxisomes →
  biogenesis score up. Recovered: gene-level Mann-Whitney p = 3.5e-3, Cohen
  d = 1.17.
- **Yeast, rapamycin vs growing (GSE149016):** TORC1 inhibition induces
  autophagy → degradation up. Recovered: degradation p = 0.005; peroxisome-
  selective ATG36 score p = 5e-4; net_direction flips negative.
- **Arabidopsis transfer (AtGenExpress, 9 stresses, 248 samples):** the same
  unchanged scorer + orthology-mapped modules; the biogenesis arm
  **discriminates stresses** (Kruskal-Wallis roots p = 1.6e-6) — heat/osmotic
  raise it, cold/salt suppress it.

## Honest limits (tiered confidence)
- **mRNA is a proxy, not an organelle count.** Read scores as "consistent with a
  shift toward biogenesis / clearance," never as turnover rates. Corroborate
  against proteome/imaging where possible.
- **Biogenesis arm = high confidence** (conserved machinery, clean signal).
- **Clearance arm = lower confidence**, especially in plants: pexophagy
  machinery overlaps bulk autophagy, and Arabidopsis has no clean ATG36
  ortholog, so `degradation_selective_score` is empty there by construction and
  the full degradation score is bulk-autophagy-confounded.
- **Cross-taxa rule:** machinery conserves, regulators diverge — only conserved
  machinery was mapped by orthology; lineage-specific regulators (yeast
  Oaf1/Pip2/Adr1) were dropped.

## Extending it
The engine is organelle- and trait-agnostic. Point `load_module` (or the
function directly) at any pair of gene-set CSVs with the same columns
(`symbol`, `resolved_systematic_id`, `submodule`, `shared_with_bulk_autophagy`)
to score a different organelle, trait, or taxon.

All data public (GEO); all IDs resolver-backed. See `geneset_provenance.json`.
