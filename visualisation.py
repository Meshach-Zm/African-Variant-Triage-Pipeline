"""
visualisation.py
----------------
All plotting functions for the African Variant Triage Pipeline.

  plot_convergence      — bar chart of per-variant convergence scores
  plot_variant_tracks   — deep-dive REF vs ALT track prediction for one variant
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
from matplotlib.patches import Patch

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

from config import ACTIONABLE_THRESHOLD, KIDNEY_ONTOLOGY_TERMS


# ---------------------------------------------------------------------------
# Convergence bar chart
# ---------------------------------------------------------------------------

def plot_convergence(report: pd.DataFrame, output_path: str):
    """
    Bar chart of convergence score per variant.
    Red bars meet or exceed the actionable threshold.
    """
    colours     = {True: "#c0392b", False: "#7f8c8d"}
    bar_colours = [colours[v] for v in report["action_recommended"]]

    fig, ax = plt.subplots(figsize=(11, 5))
    x = range(len(report))
    ax.bar(x, report["convergence_score"], color=bar_colours, width=0.6)
    ax.axhline(
        ACTIONABLE_THRESHOLD,
        linestyle="--",
        linewidth=1.0,
        color="#2c3e50",
        label=f"Actionable threshold (score >= {ACTIONABLE_THRESHOLD})",
    )

    ax.set_xticks(list(x))
    ax.set_xticklabels(report["variant_id"], rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Convergence score (max 9)")
    ax.set_xlabel("Variant")
    ax.set_title(
        "Multi-layer Evidence Convergence — African Cohort Variants\n"
        f"Red = actionable (>= {ACTIONABLE_THRESHOLD} independent evidence points)"
    )
    ax.set_ylim(0, 10)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    ax.legend(handles=[
        Patch(facecolor="#c0392b", label=f"Actionable (score >= {ACTIONABLE_THRESHOLD})"),
        Patch(facecolor="#7f8c8d", label="Not yet actionable"),
        plt.Line2D([0], [0], linestyle="--", color="#2c3e50", label="Threshold"),
    ], fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    print(f"Convergence plot saved: {output_path}")
    return fig


# ---------------------------------------------------------------------------
# Deep-dive track plot
# ---------------------------------------------------------------------------

def plot_variant_tracks(
    dna_model,
    variant: genome.Variant,
    sequence_length: int,
    output_path: str,
    ontology_terms: list[str] | None = None,
):
    """
    Deep-dive: predict REF vs ALT tracks for one variant.
    Plots RNA-seq, ATAC, and CAGE overlaid.

    Pass --deep_dive_variant <rsid> from the CLI to trigger this.
    """
    if ontology_terms is None:
        ontology_terms = KIDNEY_ONTOLOGY_TERMS

    interval = variant.reference_interval.resize(sequence_length)

    output = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=[
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.ATAC,
            dna_client.OutputType.CAGE,
        ],
        ontology_terms=ontology_terms,
    )

    ref, alt   = output.reference, output.alternate
    components = []

    if ref.rna_seq is not None:
        components.append(plot_components.OverlaidTracks(
            tdata={"REF": ref.rna_seq, "ALT": alt.rna_seq},
            colors={"REF": "dimgrey", "ALT": "#c0392b"},
            ylabel_template="RNA-SEQ: {biosample_name} ({strand})\n{name}",
        ))

    if ref.atac is not None:
        components.append(plot_components.OverlaidTracks(
            tdata={"REF": ref.atac, "ALT": alt.atac},
            colors={"REF": "steelblue", "ALT": "#e67e22"},
            ylabel_template="ATAC: {biosample_name} ({strand})\n{name}",
        ))

    if ref.cage is not None:
        components.append(plot_components.OverlaidTracks(
            tdata={"REF": ref.cage, "ALT": alt.cage},
            colors={"REF": "darkgreen", "ALT": "#8e44ad"},
            ylabel_template="CAGE: {biosample_name} ({strand})\n{name}",
        ))

    if not components:
        print(f"  No plottable tracks returned for {variant.name}.")
        return None

    plot_width = min(43008, interval.width)
    fig = plot_components.plot(
        components=components,
        interval=interval.resize(plot_width),
        annotations=[plot_components.VariantAnnotation([variant])],
    )
    fig.suptitle(
        f"Track prediction: {variant.name}  "
        f"({variant.chromosome}:{variant.position}  "
        f"{variant.reference_bases}>{variant.alternate_bases})",
        fontsize=11,
    )
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Track plot saved: {output_path}")
    return fig