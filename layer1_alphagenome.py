"""
layer2_crossref.py
------------------
Layer 2: Cross-reference lookups.

  2a — GTEx v8 single-tissue eQTL query (kidney cortex and medulla)
  2b — RegulomeDB v2 regulatory annotation

Returns one row per variant with all Layer 2 annotation columns.
"""

import time

import pandas as pd
import requests
from tqdm import tqdm

from config import API_PAUSE, GTEX_API, GTEX_KIDNEY_TISSUES, REGULOME_API
from utils import gtex_variant_id


# ---------------------------------------------------------------------------
# Layer 2a — GTEx eQTL lookup
# ---------------------------------------------------------------------------

def query_gtex_eqtl(variant_id_gtex: str, rsid: str) -> dict:
    """
    Query GTEx API v2 for significant single-tissue eQTLs in kidney tissues.

    GTEx variant ID format:  chr1_153925321_G_A_b38
    Falls back to rsID lookup if the formatted ID returns nothing.

    Returns a dict with keys:
        found         bool
        tissues       list of tissue names with a significant eQTL
        max_nes       float  (normalised effect size, largest absolute value)
        min_pvalue    float
        gene_symbols  list
    """
    result = {
        "found": False,
        "tissues": [],
        "max_nes": 0.0,
        "min_pvalue": 1.0,
        "gene_symbols": [],
    }

    for tissue in GTEX_KIDNEY_TISSUES:
        for id_param, id_value in [("variantId", variant_id_gtex), ("snpId", rsid)]:
            try:
                response = requests.get(
                    f"{GTEX_API}/association/singleTissueEqtl",
                    params={
                        id_param: id_value,
                        "tissueSiteDetailId": tissue,
                        "datasetId": "gtex_v8",
                        "itemsPerPage": 50,
                    },
                    timeout=15,
                )
                time.sleep(API_PAUSE)

                if response.status_code != 200:
                    continue

                data = response.json()
                hits = data.get("data", data.get("singleTissueEqtl", []))

                if not hits:
                    continue

                result["found"] = True
                for hit in hits:
                    nes  = float(hit.get("nes", hit.get("effectSize", 0.0)))
                    pval = float(hit.get("pValue", 1.0))
                    gene = hit.get("geneSymbol", "")

                    if abs(nes) > abs(result["max_nes"]):
                        result["max_nes"] = nes
                    if pval < result["min_pvalue"]:
                        result["min_pvalue"] = pval
                    if tissue not in result["tissues"]:
                        result["tissues"].append(tissue)
                    if gene and gene not in result["gene_symbols"]:
                        result["gene_symbols"].append(gene)

                break  # found with this ID format; skip rsID fallback

            except requests.RequestException as exc:
                print(f"  GTEx API warning for {rsid} / {tissue}: {exc}")
                continue

    return result


# ---------------------------------------------------------------------------
# Layer 2b — RegulomeDB regulatory annotation
# ---------------------------------------------------------------------------

def query_regulomedb(chrom: str, pos: int) -> dict:
    """
    Query RegulomeDB v2 for a single variant position.

    Endpoint:
        https://regulomedb.org/regulome-search/
        ?regions=chr1:153925320-153925321&genome=GRCh38&format=json

    RegulomeDB uses 0-based half-open coordinates, so we subtract 1
    from the 1-based VCF position for the start coordinate.

    Returns a dict with keys:
        found         bool
        rank          str   e.g. '1a', '2b', '7'  (1a = highest evidence)
        probability   float (0-1; higher = more likely functional)
        has_chip      bool  (TF ChIP-seq peak overlaps variant)
        has_dnase     bool  (DNase-seq peak overlaps variant)
        has_qtl       bool  (overlaps a QTL)
        has_footprint bool
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

        data    = response.json()
        regions = data.get("regions", [])
        if not regions:
            return result

        hit = regions[0]
        result["found"]       = True
        result["rank"]        = str(hit.get("regulome_score", {}).get("ranking", ""))
        result["probability"] = float(
            hit.get("regulome_score", {}).get("probability", 0.0)
        )

        annotation = hit.get("annotation", {})
        result["has_chip"]      = bool(annotation.get("ChIP", False))
        result["has_dnase"]     = bool(annotation.get("DNase", False))
        result["has_qtl"]       = bool(annotation.get("QTL", False))
        result["has_footprint"] = bool(annotation.get("Footprint", False))

    except requests.RequestException as exc:
        print(f"  RegulomeDB warning for {chrom}:{pos}: {exc}")

    return result


# ---------------------------------------------------------------------------
# Public entry point for main.py
# ---------------------------------------------------------------------------

def run_layer2(vcf: pd.DataFrame) -> pd.DataFrame:
    """
    Run GTEx and RegulomeDB lookups for every variant.
    Returns one row per variant with all Layer 2 columns.
    """
    rows = []
    for _, row in tqdm(vcf.iterrows(), total=len(vcf), desc="Layer 2: GTEx + RegulomeDB"):
        rsid    = str(row["variant_id"])
        chrom   = str(row["CHROM"])
        pos     = int(row["POS"])
        ref     = str(row["REF"])
        alt     = str(row["ALT"])
        gtex_id = gtex_variant_id(chrom, pos, ref, alt)

        gtex     = query_gtex_eqtl(gtex_id, rsid)
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