"""
utils.py
--------
Shared utility functions: VCF loading and GTEx variant ID formatting.
"""

from io import StringIO
from pathlib import Path

import pandas as pd

from config import EXAMPLE_VCF


def load_vcf(vcf_path: str | None) -> pd.DataFrame:
    """
    Load variants from a tab-separated VCF file.
    Required columns: variant_id, CHROM, POS, REF, ALT.
    CHROM must use the 'chr' prefix.  POS is 1-based (hg38).
    Falls back to the built-in example set if no path is given.
    """
    if vcf_path is None:
        print("No VCF file supplied. Using built-in example variants.")
        return pd.read_csv(StringIO(EXAMPLE_VCF), sep="\t")

    path = Path(vcf_path)
    if not path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    df = pd.read_csv(path, sep="\t", comment="#")
    required = {"variant_id", "CHROM", "POS", "REF", "ALT"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"VCF is missing required columns: {missing}")
    return df


def gtex_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """
    Build the GTEx variant ID format expected by the API.
    Example:  chr1_153925321_G_A_b38
    """
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{chrom}_{pos}_{ref}_{alt}_b38"