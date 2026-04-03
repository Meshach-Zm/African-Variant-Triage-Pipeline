"""
main.py
-------
African Variant Triage Pipeline — entry point.

Converts African-cohort VUS from statistical associations into
multi-evidence functional hypotheses, ready for wet-lab prioritisation.

Usage
-----
    python main.py --api_key YOUR_KEY_HERE
    python main.py --api_key YOUR_KEY_HERE --vcf my_variants.vcf
    python main.py --api_key YOUR_KEY_HERE \\
        --vcf my_variants.vcf \\
        --deep_dive_variant rs138493856 \\
        --output_dir results/

Installation
------------
    git clone https://github.com/google-deepmind/alphagenome.git
    pip install ./alphagenome
    pip install pandas matplotlib tqdm requests scipy
"""

import argparse
import os
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

# Load variables from .env if present (no-op if the file doesn't exist)
load_dotenv()

from alphagenome.data import genome
from alphagenome.models import dna_client

from config import KIDNEY_ONTOLOGY_TERMS, SEQUENCE_LENGTH_KEY
from convergence import build_final_report, print_action_report
from layer1_alphagenome import run_layer1, summarise_alphagenome
from layer2_crossref import run_layer2
from layer3_population import run_layer3
from utils import load_vcf
from visualisation import plot_convergence, plot_variant_tracks


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="African variant triage pipeline — Layers 1 through 3."
    )
    parser.add_argument(
        "--api_key",
        default=os.getenv("ALPHAGENOME_API_KEY"),
        help="AlphaGenome API key. Can also be set via ALPHAGENOME_API_KEY "
             "in a .env file or environment variable "
             "(obtain at deepmind.google.com/science/alphagenome).",
    )
    parser.add_argument(
        "--vcf", default=None,
        help="Tab-separated VCF file. "
             "Required columns: variant_id, CHROM, POS, REF, ALT. "
             "Optional population columns: "
             "N_REF_REF, N_REF_ALT, N_ALT_ALT, CASE_AF, CONTROL_AF.",
    )
    parser.add_argument(
        "--output_dir", default=".",
        help="Directory for CSV and plot outputs (default: current directory).",
    )
    parser.add_argument(
        "--deep_dive_variant", default=None,
        help="variant_id for a full REF vs ALT track plot "
             "(RNA-seq, ATAC, CAGE). "
             "Example: --deep_dive_variant rs138493856",
    )
    parser.add_argument(
        "--skip_alphagenome", action="store_true",
        help="Skip Layer 1 and load kidney_scores_full.csv from --output_dir. "
             "Useful for re-running Layers 2 and 3 without re-querying the API.",
    )
    args = parser.parse_args()

    if not args.api_key:
        parser.error(
            "AlphaGenome API key is required. "
            "Pass it via --api_key or set ALPHAGENOME_API_KEY in a .env file."
        )
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    vcf = load_vcf(args.vcf)
    print(f"Loaded {len(vcf)} variants.")

    # ---- Layer 1 ----
    if args.skip_alphagenome:
        cached = output_dir / "kidney_scores_full.csv"
        if not cached.exists():
            raise FileNotFoundError(
                f"--skip_alphagenome set but {cached} not found."
            )
        print(f"Skipping Layer 1. Loading cached scores from {cached}.")
        df_kidney  = pd.read_csv(cached)
        dna_model  = None
    else:
        print("\n--- Layer 1: AlphaGenome ---")
        dna_model, df_kidney = run_layer1(args.api_key, vcf)
        df_kidney.to_csv(output_dir / "kidney_scores_full.csv", index=False)
        print("  Full kidney scores saved.")

    ag_summary = summarise_alphagenome(df_kidney)

    # ---- Layer 2 ----
    print("\n--- Layer 2: GTEx + RegulomeDB ---")
    l2_df = run_layer2(vcf)
    l2_df.to_csv(output_dir / "layer2_cross_reference.csv", index=False)
    print("  Cross-reference results saved.")

    # ---- Layer 3 ----
    print("\n--- Layer 3: Population genetics ---")
    l3_df = run_layer3(vcf)
    l3_df.to_csv(output_dir / "layer3_population.csv", index=False)
    print("  Population genetics results saved.")

    # ---- Convergence ----
    report = build_final_report(ag_summary, l2_df, l3_df)
    report.to_csv(output_dir / "final_report.csv", index=False)
    print(f"\nFinal report saved: {output_dir / 'final_report.csv'}")

    print_action_report(report)
    plot_convergence(report, output_path=str(output_dir / "convergence_scores.png"))

    # ---- Optional deep-dive ----
    if args.deep_dive_variant is not None and not args.skip_alphagenome:
        row = vcf[vcf["variant_id"] == args.deep_dive_variant]
        if row.empty:
            print(
                f"WARNING: --deep_dive_variant '{args.deep_dive_variant}' "
                "not found in VCF. Skipping track plot."
            )
        else:
            row    = row.iloc[0]
            target = genome.Variant(
                chromosome=str(row["CHROM"]),
                position=int(row["POS"]),
                reference_bases=str(row["REF"]),
                alternate_bases=str(row["ALT"]),
                name=str(row["variant_id"]),
            )
            sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[SEQUENCE_LENGTH_KEY]
            plot_variant_tracks(
                dna_model=dna_client.create(args.api_key),
                variant=target,
                sequence_length=sequence_length,
                output_path=str(
                    output_dir / f"tracks_{args.deep_dive_variant}.png"
                ),
            )

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()