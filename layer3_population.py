"""
layer3_population.py
--------------------
Layer 3: Population genetics checks.

  - Hardy-Weinberg equilibrium test (chi-squared, 1 d.f.)
  - Case vs control allele-frequency enrichment flag

Optional VCF columns consumed:
    N_REF_REF, N_REF_ALT, N_ALT_ALT  — integer genotype counts in cases
    CASE_AF, CONTROL_AF               — allele frequencies

If these columns are absent the layer still runs, marking results as
insufficient_data rather than raising an error.
"""

import pandas as pd
from scipy import stats
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------------------

def hardy_weinberg_test(n_ref_ref: int, n_ref_alt: int, n_alt_alt: int) -> float:
    """
    Chi-squared Hardy-Weinberg equilibrium test.
    Returns the p-value.  Low p-value = deviation from HWE (potential signal).
    Returns 1.0 if counts are zero (no data).
    """
    n = n_ref_ref + n_ref_alt + n_alt_alt
    if n == 0:
        return 1.0

    p = (2 * n_ref_ref + n_ref_alt) / (2 * n)
    q = 1.0 - p

    expected_rr = p ** 2 * n
    expected_ra = 2 * p * q * n
    expected_aa = q ** 2 * n

    if any(e < 1e-9 for e in [expected_rr, expected_ra, expected_aa]):
        return 1.0

    chi2 = sum(
        (o - e) ** 2 / e
        for o, e in zip(
            [n_ref_ref, n_ref_alt, n_alt_alt],
            [expected_rr, expected_ra, expected_aa],
        )
    )
    return float(stats.chi2.sf(chi2, df=1))


# ---------------------------------------------------------------------------
# Public entry point for main.py
# ---------------------------------------------------------------------------

def run_layer3(vcf: pd.DataFrame) -> pd.DataFrame:
    """
    Run population genetics checks for each variant.

    Looks for optional genotype count columns in the VCF:
        N_REF_REF, N_REF_ALT, N_ALT_ALT  (integer counts in cases)
        CASE_AF, CONTROL_AF               (allele frequencies)

    If these columns are absent the layer still runs but marks
    results as insufficient_data rather than raising an error.
    """
    has_geno = {"N_REF_REF", "N_REF_ALT", "N_ALT_ALT"}.issubset(vcf.columns)
    has_af   = {"CASE_AF", "CONTROL_AF"}.issubset(vcf.columns)

    rows = []
    for _, row in tqdm(vcf.iterrows(), total=len(vcf), desc="Layer 3: Population genetics"):
        rsid = str(row["variant_id"])

        if has_geno:
            hwe_pval = hardy_weinberg_test(
                int(row["N_REF_REF"]),
                int(row["N_REF_ALT"]),
                int(row["N_ALT_ALT"]),
            )
            hwe_flag = hwe_pval < 0.001
        else:
            hwe_pval = None
            hwe_flag = None

        if has_af:
            case_af     = float(row["CASE_AF"])
            control_af  = float(row["CONTROL_AF"])
            af_diff     = case_af - control_af
            af_enriched = af_diff > 0.02
        else:
            case_af     = None
            control_af  = None
            af_diff     = None
            af_enriched = None

        rows.append({
            "variant_id":         rsid,
            "hwe_pvalue":         hwe_pval,
            "hwe_deviation":      hwe_flag,
            "case_af":            case_af,
            "control_af":         control_af,
            "af_difference":      af_diff,
            "af_enriched_cases":  af_enriched,
            "pop_data_available": has_geno or has_af,
        })

    return pd.DataFrame(rows)