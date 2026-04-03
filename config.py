"""
config.py
---------
All constants and thresholds for the African Variant Triage Pipeline.
"""

# ---------------------------------------------------------------------------
# Kidney tissue ontology terms (UBERON / CL identifiers)
# ---------------------------------------------------------------------------
# UBERON:0002113 = kidney     UBERON:0001225 = cortex of kidney
# CL:0000653     = podocyte
KIDNEY_ONTOLOGY_TERMS = [
    "UBERON:0002113",
    "UBERON:0001225",
    "CL:0000653",
]

# GTEx tissue IDs that correspond to the above kidney ontology terms.
GTEX_KIDNEY_TISSUES = [
    "Kidney_Cortex",
    "Kidney_Medulla",
]

# ---------------------------------------------------------------------------
# AlphaGenome quantile score thresholds
# ---------------------------------------------------------------------------
# Calibrated against common variants (MAF > 0.01) by DeepMind.
THRESHOLD_HIGH     = 0.90   # top 10%
THRESHOLD_MODERATE = 0.70   # top 30%

# Sequence window around each variant.
SEQUENCE_LENGTH_KEY = "SEQUENCE_LENGTH_1MB"

# ---------------------------------------------------------------------------
# External API endpoints
# ---------------------------------------------------------------------------
# GTEx API v2 base URL (v1 has been discontinued).
GTEX_API = "https://gtexportal.org/api/v2"

# RegulomeDB search endpoint (returns JSON).
REGULOME_API = "https://regulomedb.org/regulome-search/"

# Seconds to pause between API calls to respect rate limits.
API_PAUSE = 0.5

# ---------------------------------------------------------------------------
# Evidence convergence scoring thresholds
# ---------------------------------------------------------------------------
# A variant is flagged as actionable when convergence score >= this value.
ACTIONABLE_THRESHOLD = 4

# ---------------------------------------------------------------------------
# Example / fallback variants
# ---------------------------------------------------------------------------
# African-cohort VUS with hypertension association.
EXAMPLE_VCF = """variant_id\tCHROM\tPOS\tREF\tALT
rs138493856\tchr1\t153925321\tG\tA
rs73885319\tchr22\t36201698\tA\tC
rs186928913\tchr12\t111803962\tC\tT
rs141358937\tchr3\t58394738\tA\tT
"""