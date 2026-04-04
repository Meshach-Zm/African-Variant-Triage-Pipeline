# African Variant Triage Pipeline

> **VCF input → 3-Layer Evidence Convergence → Prioritised Wet-Lab Action List**

A multi-layer computational pipeline that converts African-cohort variants of uncertain significance (VUS) from statistical associations into multi-evidence functional hypotheses, ready for wet-lab prioritisation.

---

## At a Glance

| Input                                                | Process                                                                                 | Output                                                                                                           |
| ---------------------------------------------------- | --------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| Tab-separated VCF (variant_id, CHROM, POS, REF, ALT) | 3-layer evidence convergence: deep learning + empirical databases + population genetics | Ranked variant shortlist with convergence scores, evidence summaries, and specific wet-lab assay recommendations |

---

## Table of Contents

- [Overview](#overview)
- [Why This Pipeline vs CADD / VEP](#why-this-pipeline-vs-cadd--vep)
- [Pipeline Architecture](#pipeline-architecture)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Running on Google Colab](#running-on-google-colab)
- [Quick Start](#quick-start)
- [Input Format](#input-format)
- [Usage](#usage)
- [Output Files](#output-files)
- [Evidence Convergence Scoring](#evidence-convergence-scoring)
- [Validation](#validation)
- [Scalability](#scalability)
- [Configuration](#configuration)
- [Module Reference](#module-reference)
- [External APIs](#external-apis)
- [Troubleshooting](#troubleshooting)

---

## Overview

Large-scale African genomic cohort studies (e.g. AWI-Gen, H3Africa, APCDR) routinely identify variants associated with complex traits such as hypertension and chronic kidney disease. Many of these variants are classified as VUS — statistically significant in the discovery cohort but lacking functional annotation to explain *how* they act biologically.

This pipeline addresses the **GWAS-to-Function Gap** by integrating three independent evidence layers:

1. **AlphaGenome** (DeepMind) — deep-learning prediction of regulatory impact from DNA sequence
2. **GTEx + RegulomeDB** — cross-reference against measured kidney eQTL and regulatory element databases
3. **Population genetics** — Hardy-Weinberg equilibrium test and case/control allele frequency enrichment

A variant is flagged as a high-confidence functional candidate only when multiple independent lines of evidence converge on the same biological conclusion. This reduces false positives from any single model or database and produces a ranked, assay-ready shortlist — moving researchers from *a list of numbers* to *a lab plan*.

---

## Why This Pipeline vs CADD / VEP

General-purpose predictors like CADD and Ensembl VEP are designed for broad genomic annotation. This pipeline is purpose-built for a specific problem:

| Feature                | CADD / VEP          | This Pipeline                                                           |
| ---------------------- | ------------------- | ----------------------------------------------------------------------- |
| Tissue specificity     | Generic             | Kidney-specific ontology filtering (UBERON / CL)                        |
| African cohort focus   | Population-agnostic | Designed for African-cohort VUS with population-specific AF enrichment  |
| Evidence triangulation | Single score        | 3 independent layers — reduces FDR by requiring convergence             |
| Wet-lab output         | Annotation only     | Specific assay recommendations per variant (EMSA, luciferase, ATAC-seq) |
| Regulatory mechanism   | Limited             | TF ChIP, DNase, QTL, footprint, splice junction, CAGE, ATAC tracks      |

> Unlike general-purpose predictors like CADD, this pipeline integrates tissue-specific regulatory impact and population-specific allele frequency enrichment, reducing noise in non-coding regions and prioritising variants most likely to have a functional consequence in kidney tissue.

---

## Pipeline Architecture

```
VCF input
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  Layer 1 — AlphaGenome regulatory prediction                │
│  • Scores REF vs ALT across all recommended scorers         │
│  • Filters to kidney ontology terms (UBERON / CL)           │
│  • Classifies: HIGH / MODERATE / NEUTRAL per quantile       │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  Layer 2 — Cross-reference                                  │
│  2a. GTEx v8 single-tissue eQTL (Kidney Cortex + Medulla)  │
│  2b. RegulomeDB v2 regulatory annotation                    │
│      (rank, TF ChIP, DNase, QTL, footprint)                 │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  Layer 3 — Population genetics                              │
│  • Hardy-Weinberg equilibrium test (χ², 1 d.f.)            │
│  • Case vs control allele frequency enrichment (Δ > 2%)     │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  Evidence convergence scoring                               │
│  • Points awarded per independent evidence signal           │
│  • Empirical evidence (GTEx) and high-confidence DL         │
│    predictions (AlphaGenome) weighted more heavily than     │
│    population trends, ensuring biological plausibility      │
│    over statistical correlation                             │
│  • Actionable threshold: score >= 4 / 9                     │
│  • Wet-lab assay suggestions per variant                    │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
        CSV reports + convergence plot + (optional) track plot
```

### Convergence Score Chart

Variants scoring ≥ 4 are flagged actionable (red bars).

![Convergence Scores](convergence_scores.png)

### Deep-dive Track Plot

REF vs ALT predicted tracks across RNA-seq, ATAC, and CAGE for a single variant (generated when `--deep_dive_variant` is set).

![Track Plot](tracks_rs138493856.png)

---

## Project Structure

```
african_variant_pipeline/
│
├── config.py               # All constants, thresholds, and API endpoints
├── utils.py                # load_vcf, gtex_variant_id helpers
├── layer1_alphagenome.py   # AlphaGenome scoring, filtering, classification
├── layer2_crossref.py      # GTEx eQTL + RegulomeDB annotation
├── layer3_population.py    # Hardy-Weinberg test + AF enrichment
├── convergence.py          # Evidence scoring, report assembly, action report
├── visualisation.py        # Convergence bar chart + deep-dive track plot
├── main.py                 # CLI entry point — calls each layer in order
├── requirements.txt        # Python dependencies
└── README.md
```

---

## Installation

### 1. Clone this repository

```bash
git clone https://github.com/Meshach-Zm/African-Variant-Triage-Pipeline.git
cd African-Variant-Triage-Pipeline
```

### 2. Install AlphaGenome from source

> ⚠️ AlphaGenome is **not on PyPI**. This step must be done before installing `requirements.txt` or you will get `ModuleNotFoundError: No module named 'alphagenome'`.

```bash
git clone https://github.com/google-deepmind/alphagenome.git
pip install ./alphagenome
```

### 3. Install remaining dependencies

```bash
pip install -r requirements.txt
```

### 4. Obtain an AlphaGenome API key

Register at [deepmind.google.com/science/alphagenome](https://deepmind.google.com/science/alphagenome) to receive an API key.

### 5. Set up your environment variables

Copy the example env file and add your key:

```bash
cp .env.example .env
```

Then open `.env` and fill in your key:

```
ALPHAGENOME_API_KEY=your_api_key_here
```

The pipeline loads this automatically at startup via `python-dotenv`. The `.env` file is listed in `.gitignore` and will never be committed. Never paste your API key directly into the command line or source code.

---

## Running on Google Colab

Because Colab's GitHub file browser only shows `.ipynb` files, the recommended way to use this pipeline on Colab is to clone the repo directly into your session.

**Step 1 — Clone or update the repo**

```python
import os, sys

repo      = "African-Variant-Triage-Pipeline"
repo_url  = "https://github.com/Meshach-Zm/African-Variant-Triage-Pipeline.git"
repo_path = f"/content/{repo}"

if os.path.exists(repo_path):
    print("Repo already cloned — pulling latest changes...")
    !git -C {repo_path} pull origin main
else:
    print("Fresh clone...")
    !git clone {repo_url} {repo_path}

os.chdir(repo_path)
print(f"Working directory: {os.getcwd()}")

if repo_path not in sys.path:
    sys.path.insert(0, repo_path)
```

**Step 2 — Install AlphaGenome from source**

```python
alphagenome_path = "/content/alphagenome"
if not os.path.exists(alphagenome_path):
    !git clone https://github.com/google-deepmind/alphagenome.git {alphagenome_path}
    !pip install -q {alphagenome_path}
else:
    print("AlphaGenome already installed — skipping.")
```

**Step 3 — Install remaining dependencies**

```python
!pip install -q -r requirements.txt
print("Setup complete.")
```

**Step 4 — Set your API key using Colab Secrets**

Add `ALPHAGENOME_API_KEY` under **Secrets** in the Colab left sidebar (the 🔑 icon), then load it in your notebook:

```python
import os
from google.colab import userdata
os.environ["ALPHAGENOME_API_KEY"] = userdata.get("ALPHAGENOME_API_KEY")
```

**Step 5 — Run the pipeline**

```python
!python main.py --output_dir results/
```

> **Note:** Colab sessions are ephemeral. Save outputs to Google Drive to persist them across sessions:
>
> ```python
> from google.colab import drive
> drive.mount('/content/drive')
> !python main.py --output_dir /content/drive/MyDrive/variant_results/
> ```

---

## Quick Start

With `ALPHAGENOME_API_KEY` set in your `.env` file, no `--api_key` flag is needed:

```bash
# Run on built-in example variants
python main.py

# Run on your own variants
python main.py --vcf my_variants.vcf --output_dir results/

# With a deep-dive track plot for one variant
python main.py --vcf my_variants.vcf --deep_dive_variant rs138493856 --output_dir results/

# Re-run Layers 2 and 3 only — skips AlphaGenome API call, loads cached scores
python main.py --vcf my_variants.vcf --output_dir results/ --skip_alphagenome
```

You can still pass the key explicitly if preferred — it takes precedence over the env var:

```bash
python main.py --api_key YOUR_KEY_HERE --vcf my_variants.vcf
```

---

## Input Format

The pipeline accepts a **tab-separated file** with the following columns:

| Column       | Type    | Description                                        |
| ------------ | ------- | -------------------------------------------------- |
| `variant_id` | string  | rsID or any unique identifier (e.g. `rs138493856`) |
| `CHROM`      | string  | Chromosome with `chr` prefix (e.g. `chr1`)         |
| `POS`        | integer | 1-based position in hg38                           |
| `REF`        | string  | Reference allele                                   |
| `ALT`        | string  | Alternate allele                                   |

### Optional columns for Layer 3 (population genetics)

| Column       | Type    | Description                         |
| ------------ | ------- | ----------------------------------- |
| `N_REF_REF`  | integer | Homozygous reference count in cases |
| `N_REF_ALT`  | integer | Heterozygous count in cases         |
| `N_ALT_ALT`  | integer | Homozygous alternate count in cases |
| `CASE_AF`    | float   | Allele frequency in cases           |
| `CONTROL_AF` | float   | Allele frequency in controls        |

If the optional columns are absent, Layer 3 still runs but marks results as `insufficient_data` rather than raising an error.

> **African cohort note:** For Layer 3 to be most informative, `CASE_AF` and `CONTROL_AF` should reflect African-specific allele frequencies. Comparing against a global reference such as gnomAD will dilute African-specific signals. Where possible, use within-cohort frequencies from the discovery dataset (e.g. AWI-Gen, H3Africa) or an African-specific reference panel such as the African Genome Variation Project (AGVP).

### Example VCF

```
variant_id      CHROM   POS         REF  ALT
rs138493856     chr1    153925321   G    A
rs73885319      chr22   36201698    A    C
rs186928913     chr12   111803962   C    T
rs141358937     chr3    58394738    A    T
```

---

## Usage

```
python main.py [-h] [--api_key API_KEY] [--vcf VCF] [--output_dir OUTPUT_DIR]
               [--deep_dive_variant DEEP_DIVE_VARIANT] [--skip_alphagenome]
```

| Argument              | Required | Default          | Description                                        |
| --------------------- | -------- | ---------------- | -------------------------------------------------- |
| `--api_key`           | No       | env var          | AlphaGenome API key (or set ALPHAGENOME_API_KEY)   |
| `--vcf`               | No       | built-in example | Path to tab-separated VCF file                     |
| `--output_dir`        | No       | `.`              | Directory for all output files                     |
| `--deep_dive_variant` | No       | —                | `variant_id` for full REF vs ALT track plot        |
| `--skip_alphagenome`  | No       | False            | Skip Layer 1; load cached `kidney_scores_full.csv` |

---

## Output Files

All files are written to `--output_dir` (default: current directory).

| File                         | Description                                                  |
| ---------------------------- | ------------------------------------------------------------ |
| `kidney_scores_full.csv`     | Full AlphaGenome scores filtered to kidney ontology terms    |
| `layer2_cross_reference.csv` | GTEx eQTL and RegulomeDB results per variant                 |
| `layer3_population.csv`      | HWE p-values and allele frequency comparisons                |
| `final_report.csv`           | Merged three-layer results with convergence scores           |
| `convergence_scores.png`     | Bar chart of convergence score per variant                   |
| `tracks_<rsid>.png`          | REF vs ALT track plot (only if `--deep_dive_variant` is set) |

### Example outputs

**Convergence scores** — one bar per variant, red = actionable (score ≥ 4):

![Convergence Scores](convergence_scores.png)

**Deep-dive track plot** — REF (muted) vs ALT (vivid) overlaid for RNA-seq, ATAC, and CAGE:

![Track Plot](tracks_rs138493856.png)

### Example final report (4 built-in example variants)

| variant_id          | output_type      | biosample_name   | quantile_score | ag_impact | convergence_score | action_recommended |
| ------------------- | ---------------- | ---------------- | -------------- | --------- | ----------------- | ------------------ |
| chr1:153925321:G>A  | RNA_SEQ          | kidney           | 1.000          | HIGH      | 4                 | True               |
| chr22:36201698:A>C  | SPLICE_JUNCTIONS | cortex of kidney | 1.000          | HIGH      | 4                 | True               |
| chr12:111803962:C>T | RNA_SEQ          | kidney           | -1.000         | HIGH      | 4                 | True               |
| chr3:58394738:A>T   | RNA_SEQ          | kidney           | -0.990         | HIGH      | 4                 | True               |

All 4 example variants scored HIGH by AlphaGenome in kidney-specific tissue tracks, with RegulomeDB TF ChIP, DNase, and QTL evidence contributing to an actionable convergence score of 4/9. GTEx and Layer 3 population scores will populate when real cohort data with genotype counts is supplied.

The pipeline also prints a plain-text **action report** to stdout listing actionable variants, their evidence summaries, and suggested wet-lab assays.

---

## Evidence Convergence Scoring

Points are awarded independently across all three layers. The maximum possible score is **9**.

> This scoring system weights empirical evidence (GTEx, +2) and high-confidence deep learning predictions (AlphaGenome HIGH, +2) more heavily than population trends (+1 each), ensuring biological plausibility over statistical correlation.

| Evidence signal                              | Points |
| -------------------------------------------- | ------ |
| AlphaGenome quantile score ≥ 0.90 (HIGH)     | +2     |
| AlphaGenome quantile score ≥ 0.70 (MODERATE) | +1     |
| GTEx eQTL found in kidney tissue             | +2     |
| RegulomeDB rank 1x or 2x                     | +1     |
| RegulomeDB TF ChIP-seq + DNase-seq overlap   | +1     |
| RegulomeDB QTL overlap                       | +1     |
| HWE deviation in cases (p < 0.001)           | +1     |
| AF enriched in cases vs controls (Δ > 2%)    | +1     |

A variant is flagged as **actionable** when its convergence score is **≥ 4**. This threshold is configurable via `ACTIONABLE_THRESHOLD` in `config.py`.

### Suggested wet-lab assays

| Evidence flag triggered       | Suggested assay                     |
| ----------------------------- | ----------------------------------- |
| RegulomeDB TF ChIP overlap    | EMSA or CUT&RUN to test TF binding  |
| RegulomeDB DNase or footprint | ATAC-seq on isogenic cell line pair |
| GTEx eQTL found               | Luciferase reporter assay           |
| No specific flags             | Luciferase reporter or ATAC-seq     |

---

## Validation

### Sensitivity — built-in example variants

The pipeline was run on 4 African-cohort VUS with known hypertension associations. All 4 variants were correctly identified as actionable (convergence score ≥ 4), driven by AlphaGenome HIGH predictions in kidney-specific tissue tracks and RegulomeDB regulatory evidence. This confirms the pipeline's sensitivity to regulatory variants in the target tissue.

### Positive control recommendation

To validate the full scoring range on your own dataset, include at least one variant with a known kidney eQTL in GTEx v8 (e.g. variants near UMOD, SLC12A3, or CUBN) as a positive control. A true positive should score ≥ 6 (AlphaGenome HIGH + GTEx eQTL + RegulomeDB evidence), confirming that the GTEx +2 points are being awarded correctly.

### Full cohort validation

Layer 3 (population genetics) requires a VCF with genotype count columns (`N_REF_REF`, `N_REF_ALT`, `N_ALT_ALT`) or allele frequency columns (`CASE_AF`, `CONTROL_AF`). These are available from controlled-access African cohort datasets via the H3Africa Data Access Committee (DAC) at the European Genome-phenome Archive (EGA).

---

## Scalability

### API rate limiting

All external API calls include a configurable pause (`API_PAUSE = 0.5s` in `config.py`). For large variant sets this adds up:

| Variants | Estimated Layer 2 runtime |
| -------- | ------------------------- |
| 10       | ~10 seconds               |
| 100      | ~2 minutes                |
| 1,000    | ~17 minutes               |
| 10,000   | ~3 hours                  |

### Caching with `--skip_alphagenome`

Layer 1 (AlphaGenome) is the most time-consuming step. Once run, scores are cached to `kidney_scores_full.csv`. Use `--skip_alphagenome` on all subsequent runs to reload from cache rather than re-querying the API — this is the primary mechanism for iterative analysis on large variant sets.

### Recommended workflow for large cohorts

```bash
# Step 1 — Run Layer 1 once, cache the scores
python main.py --vcf cohort.vcf --output_dir results/

# Step 2 — Iterate on Layers 2 and 3 without re-querying AlphaGenome
python main.py --vcf cohort.vcf --output_dir results/ --skip_alphagenome
```

---

## Configuration

All thresholds and constants are centralised in `config.py`. Key settings:

| Constant                | Default                                          | Description                                               |
| ----------------------- | ------------------------------------------------ | --------------------------------------------------------- |
| `KIDNEY_ONTOLOGY_TERMS` | `UBERON:0002113`, `UBERON:0001225`, `CL:0000653` | Ontology terms used to filter AlphaGenome tracks          |
| `GTEX_KIDNEY_TISSUES`   | `Kidney_Cortex`, `Kidney_Medulla`                | GTEx tissue IDs queried for eQTLs                         |
| `THRESHOLD_HIGH`        | `0.90`                                           | AlphaGenome quantile threshold for HIGH impact            |
| `THRESHOLD_MODERATE`    | `0.70`                                           | AlphaGenome quantile threshold for MODERATE impact        |
| `ACTIONABLE_THRESHOLD`  | `4`                                              | Minimum convergence score to flag a variant as actionable |
| `SEQUENCE_LENGTH_KEY`   | `SEQUENCE_LENGTH_1MB`                            | Sequence window passed to AlphaGenome                     |
| `API_PAUSE`             | `0.5`                                            | Seconds between API calls (rate limiting)                 |

To adapt the pipeline to a different tissue, update `KIDNEY_ONTOLOGY_TERMS` and `GTEX_KIDNEY_TISSUES` in `config.py`. No other files need to change.

---

## Module Reference

### `config.py`
Centralised constants. Edit this file to change thresholds, tissues, or API endpoints.

### `utils.py`
- `load_vcf(vcf_path)` — loads and validates the input VCF; falls back to example variants if `None`.
- `gtex_variant_id(chrom, pos, ref, alt)` — formats a GTEx-compatible variant ID string.

### `layer1_alphagenome.py`
- `run_layer1(api_key, vcf)` — initialises the model, scores all variants, filters to kidney tissues, classifies impact. Returns `(dna_model, df_kidney)`.
- `summarise_alphagenome(df)` — collapses to one row per variant (highest absolute quantile score).

### `layer2_crossref.py`
- `run_layer2(vcf)` — runs GTEx and RegulomeDB lookups for every variant. Returns one row per variant.
- `query_gtex_eqtl(chrom, pos, rsid)` — queries GTEx v8 via `singleTissueEqtlByLocation` endpoint; confirmed field names: `nes`, `pValue`, `geneSymbol`.
- `query_regulomedb(chrom, pos)` — queries RegulomeDB v2; returns rank, probability, and feature flags from confirmed response structure.

### `layer3_population.py`
- `run_layer3(vcf)` — runs HWE test and AF enrichment check per variant. Gracefully handles missing optional columns.
- `hardy_weinberg_test(n_ref_ref, n_ref_alt, n_alt_alt)` — chi-squared HWE test, returns p-value.

### `convergence.py`
- `build_final_report(ag_summary, l2_df, l3_df)` — merges all layers and computes convergence scores.
- `compute_convergence_score(row)` — awards points per evidence signal; returns `(score, reasons)`.
- `print_action_report(report)` — prints ranked triage summary to stdout with assay suggestions.

### `visualisation.py`
- `plot_convergence(report, output_path)` — bar chart of convergence scores; red bars indicate actionable variants.
- `plot_variant_tracks(dna_model, variant, sequence_length, output_path)` — REF vs ALT overlay for RNA-seq, ATAC, and CAGE tracks.

### `main.py`
CLI entry point. Parses arguments, calls each layer in order, saves outputs, and triggers optional deep-dive.

---

## External APIs

| Service     | Endpoint                                  | Notes                                                                   |
| ----------- | ----------------------------------------- | ----------------------------------------------------------------------- |
| AlphaGenome | DeepMind API (key required)               | Scores regulatory impact from DNA sequence                              |
| GTEx        | `https://gtexportal.org/api/v2`           | Uses `singleTissueEqtlByLocation` — confirmed working endpoint          |
| RegulomeDB  | `https://regulomedb.org/regulome-search/` | GRCh38 coordinates; 0-based half-open regions; confirmed JSON structure |

All external calls include a `0.5 s` pause between requests (`API_PAUSE` in `config.py`) to respect rate limits. Timeouts are set to 15 seconds per request.

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'alphagenome'`**
AlphaGenome is not on PyPI and cannot be installed via `pip install -r requirements.txt` alone. You must clone and install it manually first — see [Installation](#installation) step 2 or the [Colab](#running-on-google-colab) instructions.

**`FileNotFoundError: --skip_alphagenome set but kidney_scores_full.csv not found`**
Run the pipeline once without `--skip_alphagenome` to generate the cached scores file.

**`WARNING: No kidney-tissue tracks returned. Using full unfiltered scores.`**
AlphaGenome did not return tracks matching the configured ontology terms. The pipeline falls back to the full score set. Check that `KIDNEY_ONTOLOGY_TERMS` in `config.py` matches the ontology identifiers supported by your AlphaGenome API version.

**GTEx or RegulomeDB API warnings in the log**
Both APIs are queried live during each run. Transient network errors are caught and logged as warnings; the affected variant will have empty Layer 2 fields but the pipeline will continue.

**`VCF is missing required columns: {...}`**
Ensure your input file is tab-separated and contains all five required columns: `variant_id`, `CHROM`, `POS`, `REF`, `ALT`. Lines beginning with `#` are treated as comments and skipped.

**AlphaGenome install fails**
Make sure you clone the repository and install with `pip install ./alphagenome` (local path install), not `pip install alphagenome`.

**Layer 3 shows all `n/a` values**
The optional population columns (`N_REF_REF`, `N_REF_ALT`, `N_ALT_ALT`, `CASE_AF`, `CONTROL_AF`) are not present in the VCF. Layer 3 runs but cannot compute HWE or AF enrichment without them. Add these columns from your cohort genotype data to activate Layer 3 scoring.