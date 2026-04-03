"""
convergence.py
--------------
Evidence convergence scoring and final report generation.

Each variant accumulates evidence points across all three layers.
A variant is flagged as a high-confidence functional candidate when
independent lines of evidence point at the same biological conclusion,
not just when a single model produces a large score.

Scoring rubric
--------------
Layer 1 — AlphaGenome
    +2  quantile score >= 0.90  (HIGH)
    +1  quantile score >= 0.70  (MODERATE)

Layer 2 — Cross-reference
    +2  GTEx eQTL found in kidney tissue
    +1  RegulomeDB rank 1x or 2x  (strong regulatory evidence)
    +1  RegulomeDB TF ChIP + DNase  (dual open-chromatin evidence)
    +1  RegulomeDB QTL flag

Layer 3 — Population genetics
    +1  HWE deviation in cases  (p < 0.001)
    +1  Allele frequency enriched in cases vs controls  (> 2%)

Maximum possible score : 9
Actionable threshold   : >= 4  (configurable via config.ACTIONABLE_THRESHOLD)
"""

import numpy as np
import pandas as pd

from config import ACTIONABLE_THRESHOLD


# ---------------------------------------------------------------------------
# Per-variant scoring
# ---------------------------------------------------------------------------

def _is_valid(value) -> bool:
    """Return False if value is None, NaN, or an empty string."""
    if value is None:
        return False
    if isinstance(value, float) and np.isnan(value):
        return False
    if isinstance(value, str) and value.strip() in ("", "nan"):
        return False
    return True


def compute_convergence_score(row: pd.Series) -> tuple[int, list[str]]:
    """
    Award evidence points for independent signals across all three layers.
    Returns (total_points, list_of_reason_strings).
    """
    points  = 0
    reasons = []

    # Layer 1 — AlphaGenome
    ag_impact = row.get("ag_impact", "NEUTRAL")
    if ag_impact == "HIGH":
        points += 2
        reasons.append(f"AlphaGenome HIGH (quantile={row['quantile_score']:.3f})")
    elif ag_impact == "MODERATE":
        points += 1
        reasons.append(f"AlphaGenome MODERATE (quantile={row['quantile_score']:.3f})")

    # Layer 2 — GTEx eQTL
    if row.get("gtex_eqtl_found") and _is_valid(row.get("gtex_max_nes")):
        points += 2
        reasons.append(
            f"GTEx eQTL in kidney "
            f"(NES={row['gtex_max_nes']:.3f}, "
            f"p={row['gtex_min_pvalue']:.2e}, "
            f"genes={row['gtex_genes']})"
        )

    # Layer 2 — RegulomeDB rank
    rank = str(row.get("regulome_rank", ""))
    if rank and rank[0] in ("1", "2"):
        points += 1
        reasons.append(
            f"RegulomeDB rank={rank} "
            f"(probability={row['regulome_probability']:.2f})"
        )

    # Layer 2 — dual chromatin evidence
    if row.get("regulome_has_chip") and row.get("regulome_has_dnase"):
        points += 1
        reasons.append("RegulomeDB: TF ChIP + DNase overlap")

    # Layer 2 — QTL flag
    if row.get("regulome_has_qtl"):
        points += 1
        reasons.append("RegulomeDB: QTL overlap")

    # Layer 3 — HWE deviation
    if row.get("hwe_deviation") and _is_valid(row.get("hwe_pvalue")):
        points += 1
        reasons.append(f"HWE deviation in cases (p={row['hwe_pvalue']:.4f})")

    # Layer 3 — AF enrichment in cases
    if row.get("af_enriched_cases") and _is_valid(row.get("case_af")) and _is_valid(row.get("control_af")):
        points += 1
        reasons.append(
            f"AF enriched in cases "
            f"(case={row['case_af']:.3f}, "
            f"control={row['control_af']:.3f})"
        )

    return points, reasons


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

def build_final_report(
    ag_summary: pd.DataFrame,
    l2_df: pd.DataFrame,
    l3_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge all three layers and compute convergence scores."""
    merged = ag_summary.merge(l2_df, on="variant_id", how="left")
    merged = merged.merge(l3_df,  on="variant_id", how="left")

    scores  = []
    reasons = []
    for _, row in merged.iterrows():
        s, r = compute_convergence_score(row)
        scores.append(s)
        reasons.append(" | ".join(r) if r else "No evidence flags")

    merged["convergence_score"]  = scores
    merged["evidence_summary"]   = reasons
    merged["action_recommended"] = merged["convergence_score"] >= ACTIONABLE_THRESHOLD

    return merged.sort_values("convergence_score", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Action report (stdout)
# ---------------------------------------------------------------------------

def _suggest_assay(row: pd.Series) -> str:
    suggestions = []
    if row.get("regulome_has_chip"):
        suggestions.append("EMSA or CUT&RUN to test TF binding")
    if row.get("regulome_has_dnase") or row.get("regulome_has_footprint"):
        suggestions.append("ATAC-seq on isogenic cell line pair")
    if row.get("gtex_eqtl_found"):
        suggestions.append("luciferase reporter assay")
    if not suggestions:
        suggestions.append("luciferase reporter or ATAC-seq depending on output type")
    return "; ".join(suggestions)


def _fmt(value, fmt=None, fallback="n/a"):
    """Format a value, returning fallback if it is None or NaN."""
    try:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return fallback
        return format(value, fmt) if fmt else str(value)
    except (TypeError, ValueError):
        return fallback


def print_action_report(report: pd.DataFrame):
    """Print a plain-text multi-layer triage report to stdout."""
    sep = "-" * 72
    print(f"\n{sep}")
    print("MULTI-LAYER TRIAGE REPORT  —  African Cohort Variants")
    print(sep)

    actionable = report[report["action_recommended"]]
    not_yet    = report[~report["action_recommended"]]

    print(f"  Total variants scored   : {len(report)}")
    print(f"  Actionable (score >= {ACTIONABLE_THRESHOLD}) : {len(actionable)}")
    print(f"  Not yet actionable      : {len(not_yet)}")
    print(sep)

    if not actionable.empty:
        print("\nACTIONABLE VARIANTS — multi-layer evidence converges:\n")
        for _, row in actionable.iterrows():
            print(
                f"  Variant          : {row['variant_id']}\n"
                f"  Convergence score: {row['convergence_score']} / 9\n"
                f"  AlphaGenome      : {row['ag_impact']}  "
                f"(quantile={_fmt(row['quantile_score'], '.3f')}, "
                f"output={_fmt(row['output_type'])}, "
                f"tissue={_fmt(row['biosample_name'])})\n"
                f"  GTEx eQTL        : "
                + ("yes — " + _fmt(row["gtex_genes"]) if row["gtex_eqtl_found"] else "not found")
                + f"\n  RegulomeDB       : rank={_fmt(row['regulome_rank'])}  "
                f"probability={_fmt(row['regulome_probability'], '.2f')}\n"
                f"  Evidence         : {row['evidence_summary']}\n"
                f"  Suggested assay  : {_suggest_assay(row)}\n"
            )

    if not not_yet.empty:
        print("NOT YET ACTIONABLE — insufficient convergent evidence:\n")
        for _, row in not_yet.iterrows():
            print(
                f"  {row['variant_id']}  "
                f"score={row['convergence_score']}  "
                f"AlphaGenome={row['ag_impact']}"
            )

    print(f"\n{sep}\n")