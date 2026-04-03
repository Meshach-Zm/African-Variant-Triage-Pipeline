"""
layer2_crossref.py
------------------
Layer 2: Cross-reference lookups.

  2a — GTEx v8 significant eQTL lookup via bulk download files.
       GTEx API v2 is gene-centric and cannot be queried by variant
       alone, so we use pre-computed significant eQTL flat files instead.
       Files are downloaded once (~50-200 MB each) and cached locally.

  2b — RegulomeDB v2 regulatory annotation (live API).
       Response structure confirmed via curl:
         data["regulome_score"]["ranking"]           — rank string e.g. "1a", "7"
         data["regulome_score"]["probability"]       — probability string
         data["features"]["ChIP"]                    — bool
         data["features"]["Chromatin_accessibility"] — bool (DNase equivalent)
         data["features"]["QTL"]                     — bool
         data["features"]["Footprint"]               — bool

Returns one row per variant with all Layer 2 annotation columns.
"""

import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

from config import API_PAUSE, GTEX_KIDNEY_TISSUES, REGULOME_API
from utils import gtex_variant_id


# ---------------------------------------------------------------------------
# GTEx bulk file URLs (Google Cloud Storage, GTEx v8)
# ---------------------------------------------------------------------------
GTEX_BULK_URLS = {
    "Kidney_Cortex": (
        "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/"
        "GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.signif_variant_gene_pairs.txt.gz"
    ),
    "Kidney_Medulla": (
        "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/"
        "GTEx_Analysis_v8_eQTL/Kidney_Medulla.v8.signif_variant_gene_pairs.txt.gz"
    ),
}

GTEX_CACHE_DIR = Path("gtex_cache")


# ---------------------------------------------------------------------------
# Layer 2a — GTEx bulk file lookup
# ---------------------------------------------------------------------------

def _download_gtex_file(tissue: str, cache_dir: Path) -> Path | None:
    """Download a GTEx significant eQTL file and cache it locally."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    url   = GTEX_BULK_URLS[tissue]
    local = cache_dir / f"{tissue}.signif_pairs.txt.gz"

    if local.exists():
        return local

    print(f"  Downloading GTEx {tissue} eQTL file (~50-200 MB, one-time)...")
    try:
        response = requests.get(url, stream=True, timeout=300)
        if response.status_code != 200:
            print(f"  WARNING: Could not download GTEx file for {tissue} "
                  f"(HTTP {response.status_code})")
            return None

        with open(local, "wb") as fh:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                fh.write(chunk)
        print(f"  Cached: {local}")
        return local

    except requests.RequestException as exc:
        print(f"  WARNING: GTEx download failed for {tissue}: {exc}")
        return None


def _load_gtex_table(tissue: str, cache_dir: Path) -> pd.DataFrame | None:
    """Load the GTEx significant eQTL table for one tissue."""
    local = _download_gtex_file(tissue, cache_dir)
    if local is None:
        return None
    try:
        return pd.read_csv(local, sep="\t", compression="gzip")
    except Exception as exc:
        print(f"  WARNING: Could not read GTEx file for {tissue}: {exc}")
        return None


def _build_gtex_index(cache_dir: Path) -> dict[str, pd.DataFrame]:
    """Download and load GTEx tables for all kidney tissues."""
    index = {}
    for tissue in GTEX_KIDNEY_TISSUES:
        if tissue not in GTEX_BULK_URLS:
            print(f"  WARNING: No bulk URL configured for tissue '{tissue}'. Skipping.")
            continue
        df = _load_gtex_table(tissue, cache_dir)
        if df is not None:
            index[tissue] = df
    return index


def query_gtex_eqtl(
    variant_id_gtex: str,
    rsid: str,
    gtex_index: dict[str, pd.DataFrame],
) -> dict:
    """
    Look up significant eQTLs for a variant in pre-loaded GTEx tables.
    Uses slope as effect size (equivalent to NES in the API).
    """
    result = {
        "found": False,
        "tissues": [],
        "max_nes": 0.0,
        "min_pvalue": 1.0,
        "gene_symbols": [],
    }

    for tissue, df in gtex_index.items():
        hits = df[df["variant_id"] == variant_id_gtex]
        if hits.empty:
            continue

        result["found"] = True
        for _, hit in hits.iterrows():
            slope = float(hit.get("slope", 0.0))
            pval  = float(hit.get("pval_nominal", 1.0))
            gene  = str(hit.get("gene_name", hit.get("gene_id", "")))

            if abs(slope) > abs(result["max_nes"]):
                result["max_nes"] = slope
            if pval < result["min_pvalue"]:
                result["min_pvalue"] = pval
            if tissue not in result["tissues"]:
                result["tissues"].append(tissue)
            if gene and gene not in result["gene_symbols"]:
                result["gene_symbols"].append(gene)

    return result


# ---------------------------------------------------------------------------
# Layer 2b — RegulomeDB regulatory annotation
# ---------------------------------------------------------------------------

def query_regulomedb(chrom: str, pos: int) -> dict:
    """
    Query RegulomeDB v2 for a single variant position.
    Uses 0-based half-open coordinates (subtract 1 from 1-based VCF pos).
    """
    result = {
        "found": False,
        "rank": "",
        "probability": 0.0,
        "has_chip": False,
        "has_dnase": False,
        "has_qtl": False,
        "has_footprint": False,
    }

    chrom  = chrom if chrom.startswith("chr") else f"chr{chrom}"
    region = f"{chrom}:{pos - 1}-{pos}"

    try:
        response = requests.get(
            REGULOME_API,
            params={"regions": region, "genome": "GRCh38", "format": "json"},
            timeout=15,
        )
        time.sleep(API_PAUSE)

        if response.status_code != 200:
            return result

        data           = response.json()
        regulome_score = data.get("regulome_score", {})
        features       = data.get("features", {})

        if not regulome_score:
            return result

        result["found"]       = True
        result["rank"]        = str(regulome_score.get("ranking", ""))
        result["probability"] = float(regulome_score.get("probability", 0.0))

        result["has_chip"]      = bool(features.get("ChIP", False))
        result["has_dnase"]     = bool(features.get("Chromatin_accessibility", False))
        result["has_qtl"]       = bool(features.get("QTL", False))
        result["has_footprint"] = bool(features.get("Footprint", False))

    except requests.RequestException as exc:
        print(f"  RegulomeDB warning for {chrom}:{pos}: {exc}")

    return result


# ---------------------------------------------------------------------------
# Public entry point for main.py
# ---------------------------------------------------------------------------

def run_layer2(vcf: pd.DataFrame, cache_dir: Path = GTEX_CACHE_DIR) -> pd.DataFrame:
    """
    Run GTEx and RegulomeDB lookups for every variant.
    GTEx uses pre-downloaded bulk files; RegulomeDB uses the live API.
    Returns one row per variant with all Layer 2 columns.
    """
    print("  Loading GTEx bulk eQTL files (downloads on first run)...")
    gtex_index = _build_gtex_index(cache_dir)

    if not gtex_index:
        print("  WARNING: No GTEx data loaded. GTEx scores will be empty.")

    rows = []
    for _, row in tqdm(vcf.iterrows(), total=len(vcf), desc="Layer 2: GTEx + RegulomeDB"):
        rsid    = str(row["variant_id"])
        chrom   = str(row["CHROM"])
        pos     = int(row["POS"])
        ref     = str(row["REF"])
        alt     = str(row["ALT"])
        gtex_id = gtex_variant_id(chrom, pos, ref, alt)

        gtex     = query_gtex_eqtl(gtex_id, rsid, gtex_index)
        regulome = query_regulomedb(chrom, pos)

        rows.append({
            "variant_id":             rsid,
            "gtex_eqtl_found":        gtex["found"],
            "gtex_tissues":           ", ".join(gtex["tissues"]),
            "gtex_max_nes":           gtex["max_nes"],
            "gtex_min_pvalue":        gtex["min_pvalue"],
            "gtex_genes":             ", ".join(gtex["gene_symbols"]),
            "regulome_found":         regulome["found"],
            "regulome_rank":          regulome["rank"],
            "regulome_probability":   regulome["probability"],
            "regulome_has_chip":      regulome["has_chip"],
            "regulome_has_dnase":     regulome["has_dnase"],
            "regulome_has_qtl":       regulome["has_qtl"],
            "regulome_has_footprint": regulome["has_footprint"],
        })

    return pd.DataFrame(rows)