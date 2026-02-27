"""
Generate high-resolution mutation figures and statistics for cBioPortal data.
"""

import io
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import Tuple, Optional

# High-resolution settings for publication-ready figures
DPI = 300
FIG_SIZE = (12, 8)
PALETTE = "muted"


def _setup_style():
    """Apply publication-quality style."""
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
    })


def mutation_type_distribution(df: pd.DataFrame, top_n: int = 15) -> Tuple[plt.Figure, pd.DataFrame]:
    """Bar chart of mutation type distribution with statistics table."""
    _setup_style()
    if df.empty or "mutationType" not in df.columns:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No mutation data available", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    counts = df["mutationType"].value_counts().head(top_n)
    stats_df = counts.reset_index()
    stats_df.columns = ["Mutation Type", "Count"]
    stats_df["Percentage"] = (stats_df["Count"] / stats_df["Count"].sum() * 100).round(2)

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    bars = ax.barh(range(len(counts)), counts.values, color=sns.color_palette(PALETTE, len(counts)))
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index, fontsize=10)
    ax.set_xlabel("Number of Mutations")
    ax.set_title("Mutation Type Distribution", fontweight="bold")
    ax.invert_yaxis()
    plt.tight_layout()
    return fig, stats_df


def gene_mutation_frequency(df: pd.DataFrame, top_n: int = 20, title: Optional[str] = None) -> Tuple[plt.Figure, pd.DataFrame]:
    """Bar chart of most frequently mutated genes with statistics."""
    _setup_style()
    gene_col = "geneSymbol" if "geneSymbol" in df.columns else "hugoGeneSymbol"
    if df.empty or gene_col not in df.columns:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No gene data available", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    gene_counts = df[gene_col].value_counts().head(top_n)
    stats_df = gene_counts.reset_index()
    stats_df.columns = ["Gene", "Mutation Count"]
    stats_df["Rank"] = range(1, len(stats_df) + 1)
    total = gene_counts.sum()
    stats_df["Percentage"] = (stats_df["Mutation Count"] / total * 100).round(2)

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    colors = sns.color_palette("Blues_r", len(gene_counts))
    ax.barh(range(len(gene_counts)), gene_counts.values, color=colors)
    ax.set_yticks(range(len(gene_counts)))
    ax.set_yticklabels(gene_counts.index)
    ax.set_xlabel("Number of Mutations")
    ax.set_title(title or "Top Mutated Genes", fontweight="bold")
    ax.invert_yaxis()
    plt.tight_layout()
    return fig, stats_df


def hereditary_genes_analysis(df: pd.DataFrame, top_n: int = 20, **kwargs) -> Tuple[plt.Figure, pd.DataFrame]:
    """Mutation frequency for hereditary breast cancer genes (BRCA1, BRCA2, PALB2, etc.)."""
    from cbioportal_client import HEREDITARY_BREAST_CANCER_GENES
    id_to_symbol = {v: k for k, v in HEREDITARY_BREAST_CANCER_GENES.items()}
    hbc_ids = set(HEREDITARY_BREAST_CANCER_GENES.values())
    hbc_symbols = set(HEREDITARY_BREAST_CANCER_GENES.keys())
    if df.empty:
        _setup_style()
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No mutation data", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()
    if "entrezGeneId" in df.columns:
        sub = df[df["entrezGeneId"].astype(float).fillna(-1).astype(int).isin(hbc_ids)].copy()
        if not sub.empty:
            sub["hugoGeneSymbol"] = sub["entrezGeneId"].apply(
                lambda x: id_to_symbol.get(int(float(x)) if pd.notna(x) else -1, str(int(x)) if pd.notna(x) else "")
            )
    elif "geneSymbol" in df.columns:
        sub = df[df["geneSymbol"].astype(str).str.upper().isin(hbc_symbols)]
    elif "hugoGeneSymbol" in df.columns:
        sub = df[df["hugoGeneSymbol"].astype(str).str.upper().isin(hbc_symbols)]
    else:
        sub = pd.DataFrame()
    if sub.empty:
        _setup_style()
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No hereditary breast cancer gene mutations found", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()
    return gene_mutation_frequency(sub, top_n=top_n, title="Hereditary Breast Cancer Genes")


def variant_classification_distribution(df: pd.DataFrame, top_n: int = 12) -> Tuple[plt.Figure, pd.DataFrame]:
    """Pie/bar chart of variant classification."""
    _setup_style()
    col = "variantType" if "variantType" in df.columns else "mutationType"
    if col not in df.columns:
        col = [c for c in df.columns if "variant" in c.lower() or "type" in c.lower()]
        col = col[0] if col else None
    if not col or df.empty:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No variant classification data", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    counts = df[col].value_counts().head(top_n)
    stats_df = counts.reset_index()
    stats_df.columns = ["Classification", "Count"]
    stats_df["Percentage"] = (stats_df["Count"] / stats_df["Count"].sum() * 100).round(2)

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    colors = sns.color_palette(PALETTE, len(counts))
    wedges, texts, autotexts = ax.pie(counts.values, labels=counts.index, autopct="%1.1f%%", colors=colors)
    for t in texts:
        t.set_fontsize(9)
    ax.set_title("Variant Classification Distribution", fontweight="bold")
    plt.tight_layout()
    return fig, stats_df


def sample_mutation_burden(df: pd.DataFrame, top_n: int = 25) -> Tuple[plt.Figure, pd.DataFrame]:
    """Mutation burden per sample."""
    _setup_style()
    sample_col = "sampleId" if "sampleId" in df.columns else "uniqueSampleKey"
    if df.empty or sample_col not in df.columns:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No sample data available", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    sample_counts = df[sample_col].value_counts().head(top_n)
    stats_df = sample_counts.reset_index()
    stats_df.columns = ["Sample ID", "Mutation Count"]
    stats_df["Rank"] = range(1, len(stats_df) + 1)

    fig, ax = plt.subplots(figsize=(12, max(6, len(sample_counts) * 0.25)))
    ax.barh(range(len(sample_counts)), sample_counts.values, color=sns.color_palette("viridis", len(sample_counts)))
    ax.set_yticks(range(len(sample_counts)))
    ax.set_yticklabels(sample_counts.index, fontsize=9)
    ax.set_xlabel("Number of Mutations")
    ax.set_title("Mutation Burden by Sample", fontweight="bold")
    ax.invert_yaxis()
    plt.tight_layout()
    return fig, stats_df


def oncoprint_style_matrix(df: pd.DataFrame, top_genes: int = 20, top_samples: int = 30) -> Tuple[plt.Figure, pd.DataFrame]:
    """Heatmap-style mutation matrix (gene x sample)."""
    _setup_style()
    gene_col = "geneSymbol" if "geneSymbol" in df.columns else "hugoGeneSymbol"
    sample_col = "sampleId" if "sampleId" in df.columns else "uniqueSampleKey"
    if df.empty or gene_col not in df.columns or sample_col not in df.columns:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "Insufficient data for matrix", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    top_genes_list = df[gene_col].value_counts().head(top_genes).index.tolist()
    top_samples_list = df[sample_col].value_counts().head(top_samples).index.tolist()
    sub = df[df[gene_col].isin(top_genes_list) & df[sample_col].isin(top_samples_list)]

    if sub.empty:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No data for matrix", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    pivot = sub.pivot_table(
        index=gene_col, columns=sample_col, aggfunc="size", fill_value=0
    ).clip(upper=1)

    stats_df = pivot.sum(axis=1).sort_values(ascending=False).reset_index()
    stats_df.columns = ["Gene", "Samples with Mutation"]
    stats_df["Mutation Rate (%)"] = (stats_df["Samples with Mutation"] / len(top_samples_list) * 100).round(1)

    fig, ax = plt.subplots(figsize=(14, max(8, len(top_genes_list) * 0.4)))
    sns.heatmap(pivot, cmap="YlOrRd", cbar_kws={"label": "Mutation"}, ax=ax, linewidths=0.5)
    ax.set_title("Mutation Matrix (Gene × Sample)", fontweight="bold")
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Gene")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.tight_layout()
    return fig, stats_df


def summary_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """Generate overall statistics table."""
    if df.empty:
        return pd.DataFrame({"Metric": ["No data"], "Value": ["—"]})

    stats = []
    gene_col = "geneSymbol" if "geneSymbol" in df.columns else "hugoGeneSymbol"
    sample_col = "sampleId" if "sampleId" in df.columns else "uniqueSampleKey"

    stats.append(("Total Mutations", len(df)))
    if sample_col in df.columns:
        stats.append(("Unique Samples", df[sample_col].nunique()))
    if gene_col in df.columns:
        stats.append(("Unique Genes", df[gene_col].nunique()))
    if "mutationType" in df.columns:
        stats.append(("Mutation Types", df["mutationType"].nunique()))

    return pd.DataFrame(stats, columns=["Metric", "Value"])


ANALYSIS_TYPES = {
    "Mutation Type Distribution": mutation_type_distribution,
    "Top Mutated Genes": gene_mutation_frequency,
    "Hereditary Breast Cancer Genes": hereditary_genes_analysis,
    "Variant Classification": variant_classification_distribution,
    "Sample Mutation Burden": sample_mutation_burden,
    "Mutation Matrix (Gene × Sample)": oncoprint_style_matrix,
}
