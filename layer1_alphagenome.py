"""
layer1_alphagenome.py
---------------------
Layer 1: AlphaGenome regulatory prediction via the DeepMind API.

Scores each variant across all recommended scorers, filters results to
kidney-relevant ontology terms, classifies impact, and returns a
one-row-per-variant summary.
"""

import numpy as np
import pandas as pd
from tqdm import tqdm

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

from config import (
    KIDNEY_ONTOLOGY_TERMS,
    SEQUENCE_LENGTH_KEY,
    THRESHOLD_HIGH,
    THRESHOLD_MODERATE,
)


# ---------------------------------------------------------------------------
# Scorer construction
# ---------------------------------------------------------------------------

def build_scorers() -> list:
    """Return AlphaGenome's full recommended scorer set."""
    return list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def score_variants_alphagenome(
    dna_model,
    vcf: pd.DataFrame,
    sequence_length: int,
    scorers: list,
    organism=dna_client.Organism.HOMO_SAPIENS,
) -> pd.DataFrame:
    """
    Run AlphaGenome scoring for every variant in the VCF.
    Returns a tidy dataframe: one row per (variant x scorer x track).
    """
    results = []
    for _, row in tqdm(vcf.iterrows(), total=len(vcf), desc="Layer 1: AlphaGenome"):
        variant = genome.Variant(
            chromosome=str(row["CHROM"]),
            position=int(row["POS"]),
            reference_bases=str(row["REF"]),
            alternate_bases=str(row["ALT"]),
            name=str(row["variant_id"]),
        )
        interval = variant.reference_interval.resize(sequence_length)
        variant_scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=scorers,
            organism=organism,
        )
        results.append(variant_scores)
    return variant_scorers.tidy_scores(results)


# ---------------------------------------------------------------------------
# Filtering and classification
# ---------------------------------------------------------------------------

def filter_to_kidney(df: pd.DataFrame) -> pd.DataFrame:
    """Subset AlphaGenome scores to kidney ontology terms."""
    filtered = df[df["ontology_curie"].isin(KIDNEY_ONTOLOGY_TERMS)].copy()
    if filtered.empty:
        print(
            "  WARNING: No kidney-tissue tracks returned. "
            "Using full unfiltered scores."
        )
        return df.copy()
    return filtered


def classify_alphagenome_impact(df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify each row by absolute quantile score.
    HIGH >= 0.90 / MODERATE >= 0.70 / NEUTRAL below.
    """
    abs_q = df["quantile_score"].abs()
    df["ag_impact"] = np.select(
        [abs_q >= THRESHOLD_HIGH, abs_q >= THRESHOLD_MODERATE],
        ["HIGH", "MODERATE"],
        default="NEUTRAL",
    )
    return df


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def summarise_alphagenome(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse to one row per variant: keep the track with the
    highest absolute quantile score.
    """
    idx = df.groupby("variant_id")["quantile_score"].apply(
        lambda x: x.abs().idxmax()
    )
    cols = [
        "variant_id", "output_type", "biosample_name",
        "raw_score", "quantile_score", "ag_impact",
    ]
    summary = df.loc[idx, cols].copy()
    summary = summary.sort_values("quantile_score", key=abs, ascending=False)
    return summary.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Public entry point for main.py
# ---------------------------------------------------------------------------

def run_layer1(api_key: str, vcf: pd.DataFrame) -> pd.DataFrame:
    """
    Initialise the AlphaGenome model, score all variants, filter to kidney
    tissues, classify impact, and return a per-variant summary dataframe.
    """
    dna_model       = dna_client.create(api_key)
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[SEQUENCE_LENGTH_KEY]
    scorers         = build_scorers()

    df_all    = score_variants_alphagenome(dna_model, vcf, sequence_length, scorers)
    df_kidney = filter_to_kidney(df_all)
    df_kidney = classify_alphagenome_impact(df_kidney)

    return dna_model, df_kidney
