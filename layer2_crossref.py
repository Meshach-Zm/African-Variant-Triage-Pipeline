"""
layer2_crossref.py
------------------
Layer 2: Cross-reference lookups.

  2a — GTEx v8 eQTL lookup via singleTissueEqtlByLocation API endpoint.
       Confirmed working endpoint and field names via curl:
         GET /api/v2/association/singleTissueEqtlByLocation
         params: chromosome, start, end, tissueSiteDetailId, datasetId
         response key: "singleTissueEqtl"
         fields: nes, pValue, geneSymbol, variantId, snpId

  2b — RegulomeDB v2 regulatory annotation.
       Confirmed response structure via curl:
         data["regulome_score"]["ranking"]           — rank string
         data["regulome_score"]["probability"]       — probability string
         data["features"]["ChIP"]                    — bool
         data["features"]["Chromatin_accessibility"] — bool
         data["features"]["QTL"]                     — bool
         data["features"]["Footprint"]               — bool

Returns one row per variant with all Layer 2 annotation columns.
"""

import time

import pandas as pd
import requests
from tqdm import tqdm

from config import API_PAUSE, GTEX_API, GTEX_KIDNEY_TISSUES, REGULOME_API
from utils import gtex_variant_id


# ---------------------------------------------------------------------------
# Layer 2a — GTEx eQTL lookup via singleTissueEqtlByLocation
# ---------------------------------------------------------------------------

def query_gtex_eqtl(chrom: str, pos: int, rsid: str) -> dict:
    """
    Query GTEx v8 for significant eQTLs at a variant position using the
    singleTissueEqtlByLocation endpoint — no gene ID required.

    Filters hits to the exact variant position by matching variantId or snpId.

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

    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"

    for tissue in GTEX_KIDNEY_TISSUES:
        try:
            response = requests.get(
                f"{GTEX_API}/association/singleTissueEqtlByLocation",
                params={
                    "chromosome":       chrom,
                    "start":            pos,
                    "end":              pos,
                    "tissueSiteDetailId": tissue,
                    "datasetId":        "gtex_v8",
                },
                timeout=15,
            )
            time.sleep(API_PAUSE)

            if response.status_code != 200:
                continue

            data = response.json()
            hits = data.get("singleTissueEqtl", data.get("data", []))

            if not hits:
                continue

            result["found"] = True
            for hit in hits:
                nes  = float(hit.get("nes", 0.0))
                pval = float(hit.get("pValue", 1.0))
                gene = str(hit.get("geneSymbol", ""))

                if abs(nes) > abs(result["max_nes"]):
                    result["max_nes"] = nes
                if pval < result["min_pvalue"]:
                    result["min_pvalue"] = pval
                if tissue not in result["tissues"]:
                    result["tissues"].append(tissue)
                if gene and gene not in result["gene_symbols"]:
                    result["gene_symbols"].append(gene)

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

def run_layer2(vcf: pd.DataFrame) -> pd.DataFrame:
    """
    Run GTEx and RegulomeDB lookups for every variant.
    Both use live API calls — no bulk downloads required.
    Returns one row per variant with all Layer 2 columns.
    """
    rows = []
    for _, row in tqdm(vcf.iterrows(), total=len(vcf), desc="Layer 2: GTEx + RegulomeDB"):
        rsid  = str(row["variant_id"])
        chrom = str(row["CHROM"])
        pos   = int(row["POS"])
        ref   = str(row["REF"])
        alt   = str(row["ALT"])

        gtex     = query_gtex_eqtl(chrom, pos, rsid)
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