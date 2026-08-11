"""
run_dynamics.py — load a peroxisome gene module + score an expression matrix.

Usage (as a library):
    from score_organelle_dynamics import score_organelle_dynamics
    from run_dynamics import load_module
    bio, deg, deg_sel = load_module("yeast")          # or "arabidopsis"
    scores = score_organelle_dynamics(expr, bio, deg,
                                      degradation_selective_ids=deg_sel,
                                      method="rank", n_perm=1000)

`expr` is a genes x samples pandas DataFrame whose index holds gene IDs in the
module namespace: SGD systematic names (e.g. YDR329C) for yeast, AGI loci
(e.g. AT3G47430) for Arabidopsis. Values may be raw or log; the scorer
log2(x+1)-normalizes if the matrix does not already look log-scaled.

The score is DIRECTIONAL: net_direction = biogenesis_score - degradation_score.
Positive => consistent with a shift toward peroxisome BIOGENESIS (build-up);
negative => consistent with a shift toward CLEARANCE. mRNA is a proxy, not an
organelle count — never read these as turnover rates.

Swap in a different organelle/trait by pointing load_module at other gene-set
CSVs with the same columns (symbol, resolved_systematic_id, submodule,
shared_with_bulk_autophagy). The engine is organelle-agnostic.
"""
import os
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))


def _find_genesets():
    """Locate the gene-set directory.

    Ships at <repo>/data/genesets; also supported is a scorer-local genesets/
    copy, so the kit works both inside this repo and when the scorer folder is
    vendored somewhere else on its own.
    """
    candidates = [
        os.path.join(_HERE, "genesets"),
        os.path.join(_HERE, os.pardir, "data", "genesets"),
        os.path.join(os.getcwd(), "data", "genesets"),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return os.path.abspath(c)
    raise FileNotFoundError(
        "Could not find the gene-set directory. Looked in:\n  "
        + "\n  ".join(os.path.abspath(c) for c in candidates)
        + "\nSet ORGANELLE_GENESETS to point at it."
    )


_GS = os.environ.get("ORGANELLE_GENESETS") or _find_genesets()

# tolerate the common spellings people actually type on the command line
_ORGANELLE_ALIASES = {
    "mitochondria": "mitochondrion",
    "mito": "mitochondrion",
    "mitochondrial": "mitochondrion",
    "peroxisomes": "peroxisome",
    "lysosomes": "lysosome",
    "chloroplasts": "chloroplast",
    "lipid-droplet": "lipid_droplet",
    "lipid droplet": "lipid_droplet",
    "ld": "lipid_droplet",
    "endoplasmic_reticulum": "er",
}
_SPECIES_ALIASES = {
    "arabidopsis_thaliana": "arabidopsis",
    "at": "arabidopsis",
    "hsapiens": "human",
    "homo_sapiens": "human",
    "scerevisiae": "yeast",
    "saccharomyces_cerevisiae": "yeast",
    "osativa": "rice",
    "oryza_sativa": "rice",
}

# keyed by (organelle, species) -> (biogenesis_csv, degradation_csv)
_FILES = {
    ("peroxisome", "yeast"):       ("peroxisome_biogenesis_yeast.csv",       "peroxisome_pexophagy_yeast.csv"),
    ("peroxisome", "arabidopsis"): ("peroxisome_biogenesis_arabidopsis.csv", "peroxisome_pexophagy_arabidopsis.csv"),
    ("mitochondrion", "yeast"):    ("mito_biogenesis_yeast.csv",             "mito_mitophagy_yeast.csv"),
    ("chloroplast", "arabidopsis"):("chloroplast_biogenesis_arabidopsis.csv","chloroplast_chlorophagy_arabidopsis.csv"),
    ("chloroplast", "rice"):       ("chloroplast_biogenesis_rice.csv",       "chloroplast_chlorophagy_rice.csv"),
    ("er", "yeast"):               ("er_biogenesis_yeast.csv",               "er_erphagy_yeast.csv"),
    # fill cells (orthology transfers onto public stress data; plant selective arms empty by construction)
    ("peroxisome", "rice"):        ("peroxisome_biogenesis_rice.csv",        "peroxisome_pexophagy_rice.csv"),
    ("mitochondrion", "arabidopsis"):("mito_biogenesis_arabidopsis.csv",     "mito_mitophagy_arabidopsis.csv"),
    ("mitochondrion", "rice"):     ("mito_biogenesis_rice.csv",              "mito_mitophagy_rice.csv"),
    ("er", "arabidopsis"):         ("er_biogenesis_arabidopsis.csv",         "er_erphagy_arabidopsis.csv"),
    ("er", "rice"):                ("er_biogenesis_rice.csv",                "er_erphagy_rice.csv"),
    # human mitochondrion: the disease-relevant capstone. Primary biogenesis score is MACHINERY-ONLY
    # (master regulators PGC1a/NRF1/TFAM kept in mito_regulators_human.csv, reported separately to
    # avoid the exercise-response circularity). The selective sub-score is genuinely populated here
    # (PINK1/PRKN/NIX/BNIP3/FUNDC1 — the arm empty in every plant), OPTN/CALCOCO2 flagged bulk-shared.
    ("mitochondrion", "human"):    ("mito_biogenesis_machinery_human.csv",   "mito_mitophagy_human.csv"),
    # human organelle expansion (all validated scout-first on public data with known-direction positive controls):
    # LYSOSOME — TFEB/CLEAR biogenesis + lysophagy (6 selective: LGALS3/8/9, TRIM16, UBE2QL1, FBXO27). The mito-lysosome
    #   coupled fingerprint (r=+0.88 under shared lysosomal-damage stress) is the coupled-organelle novelty.
    ("lysosome", "human"):         ("lysosome_biogenesis_human.csv",         "lysosome_lysophagy_human.csv"),
    # ER (human) — biogenesis validated HIGH (tunicamycin, 2 cell lines); 7 selective ER-phagy receptors
    #   (FAM134B/RETREG1, TEX264, CCPG1, RTN3, ATL3, RETREG2/3), clearance PARTIAL (mRNA-underpowered, post-translational).
    ("er", "human"):               ("er_biogenesis_human.csv",               "er_erphagy_human.csv"),
    # LIPID DROPLET (human) — BIOGENESIS-ONLY: lipophagy is bulk autophagy with no clean LD-selective receptor, so the
    #   clearance file is empty by construction (same 'empty selective arm' the platform records for plant organelles).
    ("lipid_droplet", "human"):    ("lipid_biogenesis_human.csv",            "lipid_lipophagy_human.csv"),
    # FERRITIN / IRON (human) — storage arm validated (FAC iron loading -> FTH1/FTL up); ferritinophagy arm has ONE
    #   selective receptor (NCOA4) but is NOT mRNA-validatable (NCOA4 post-translationally controlled) — honest partial.
    ("ferritin", "human"):         ("ferritin_biogenesis_human.csv",         "ferritin_ferritinophagy_human.csv"),
}

def load_module(species, organelle="peroxisome", id_col="resolved_systematic_id"):
    """Return (biogenesis_ids, degradation_ids, degradation_selective_ids).

    Keyed by (organelle, species). degradation_selective_ids = degradation genes
    with shared_with_bulk_autophagy == False (the organelle-SELECTIVE receptors,
    e.g. peroxisome ATG36 or mitochondrion ATG32/ATG33). This is empty by
    construction for peroxisome+arabidopsis (no clean ATG36 ortholog in plants) —
    a lower-confidence, bulk-autophagy-confounded clearance arm there. It is
    non-empty for mitochondrion+yeast (Atg32 is a genuine mitophagy receptor).
    """
    organelle = _ORGANELLE_ALIASES.get(str(organelle).lower(), str(organelle).lower())
    species = _SPECIES_ALIASES.get(str(species).lower(), str(species).lower())
    key = (organelle, species)
    if key not in _FILES:
        raise ValueError(f"(organelle, species) must be one of {list(_FILES)}")
    bio_f, deg_f = _FILES[key]
    bio = pd.read_csv(os.path.join(_GS, bio_f))
    deg = pd.read_csv(os.path.join(_GS, deg_f))
    def _flag(df):
        col = df.get("shared_with_bulk_autophagy")
        if col is None:
            return pd.Series([False]*len(df))
        return col.astype(str).str.lower().isin(["true","1","yes"])
    bio_ids = bio[id_col].dropna().astype(str).tolist()
    deg_ids = deg[id_col].dropna().astype(str).tolist()
    sel = deg.loc[~_flag(deg), id_col].dropna().astype(str).tolist()
    return bio_ids, deg_ids, sel


def _self_check():
    """Load every shipped module; show where the selective-clearance arm exists."""
    print(f"gene sets: {_GS}\n")
    print(f"{'organelle':14s}{'species':12s}{'biogenesis':>11s}{'degradation':>13s}{'selective':>11s}")
    print("-" * 61)
    for organelle, species in _FILES:
        b, d, s = load_module(species, organelle)
        print(f"{organelle:14s}{species:12s}{len(b):>11d}{len(d):>13d}{len(s):>11d}")
    print("\nA populated 'selective' column means a genuine selective-autophagy receptor")
    print("exists in that species; 0 means none is known — empty by construction, not a bug.")


def main(argv=None):
    import argparse

    combos = ", ".join(f"{o}/{s}" for o, s in _FILES)
    ap = argparse.ArgumentParser(
        prog="run_dynamics.py",
        description="Score directional organelle dynamics (biogenesis - clearance) "
                    "for every sample in an expression matrix.",
        epilog=f"available organelle/species combinations: {combos}",
    )
    ap.add_argument("--expression", "-e", metavar="CSV",
                    help="genes x samples matrix; first column = gene IDs in the module namespace")
    ap.add_argument("--organelle", default="peroxisome", help="default: peroxisome")
    ap.add_argument("--organism", "--species", dest="organism", default="yeast",
                    help="default: yeast")
    ap.add_argument("--out", "-o", metavar="CSV", help="write the score table here (default: stdout)")
    ap.add_argument("--method", choices=["rank", "zscore"], default="rank",
                    help="rank (default, works on a single sample) or zscore (needs >=2 samples)")
    ap.add_argument("--n-perm", type=int, default=1000,
                    help="permutations for the empirical null; 0 disables it (default: 1000)")
    ap.add_argument("--self-check", action="store_true",
                    help="load every shipped gene module and print its size, then exit")
    args = ap.parse_args(argv)

    if args.self_check or not args.expression:
        _self_check()
        if not args.self_check:
            print("\nNo --expression given. Pass a matrix to score one, e.g.:")
            print("  python run_dynamics.py --expression matrix.csv "
                  "--organelle mitochondrion --organism human")
        return 0

    from score_organelle_dynamics import score_organelle_dynamics

    bio, deg, deg_sel = load_module(args.organism, args.organelle)
    expr = pd.read_csv(args.expression, index_col=0)
    scores = score_organelle_dynamics(
        expr, bio, deg,
        degradation_selective_ids=deg_sel,
        method=args.method,
        n_perm=args.n_perm,
    )
    if args.out:
        scores.to_csv(args.out)
        print(f"wrote {len(scores)} sample scores to {args.out}")
    else:
        cols = ["biogenesis_score", "degradation_score", "net_direction",
                "biogenesis_emp_p", "degradation_emp_p"]
        print(scores[cols].round(4).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
