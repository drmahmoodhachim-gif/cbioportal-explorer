"""
Generate high-resolution mutation figures and statistics for cBioPortal data.
"""

import io
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import Tuple, Optional

from cbioportal_client import HEREDITARY_BREAST_CANCER_GENES, _add_gene_symbols

# High-resolution settings for publication-ready figures
DPI = 300
FIG_SIZE = (12, 8)
PALETTE = "muted"


def _setup_style():
    """Apply publication-quality style."""
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("seaborn-whitegrid")
    except OSError:
        plt.style.use("ggplot")
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
    if gene_col not in df.columns and "entrezGeneId" in df.columns:
        gene_col = "entrezGeneId"
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
    ax.set_title((title or "Top Mutated Genes"), fontweight="bold")
    ax.invert_yaxis()
    plt.tight_layout()
    return fig, stats_df


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


def oncoprint_style_matrix(df: pd.DataFrame, top_genes: int = 20, top_samples: int = 30, top_n: Optional[int] = None, **kwargs) -> Tuple[plt.Figure, pd.DataFrame]:
    """Heatmap-style mutation matrix (gene x sample)."""
    if top_n is not None:
        top_genes, top_samples = top_n, min(top_n * 2, 60)
    _setup_style()
    gene_col = "geneSymbol" if "geneSymbol" in df.columns else "hugoGeneSymbol"
    if gene_col not in df.columns and "entrezGeneId" in df.columns:
        gene_col = "entrezGeneId"
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


def hereditary_genes_analysis(df: pd.DataFrame, top_n: int = 20) -> Tuple[plt.Figure, pd.DataFrame]:
    """Focus on hereditary breast cancer genes (BRCA1, BRCA2, PALB2, etc.)."""
    df = _add_gene_symbols(df.copy())
    gene_col = "geneSymbol" if "geneSymbol" in df.columns else "hugoGeneSymbol"
    if gene_col not in df.columns and "entrezGeneId" in df.columns:
        gene_col = "entrezGeneId"
    if gene_col not in df.columns:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No gene data available", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()
    target_entrez = set(HEREDITARY_BREAST_CANCER_GENES.values())
    if gene_col == "entrezGeneId":
        def _in_target(x):
            try:
                return pd.notna(x) and int(float(x)) in target_entrez
            except (ValueError, TypeError):
                return False
        sub = df[df[gene_col].apply(_in_target)]
    else:
        target_symbols = set(HEREDITARY_BREAST_CANCER_GENES.keys())
        sub = df[df[gene_col].isin(target_symbols)]
    if sub.empty:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No hereditary gene mutations in this dataset", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()
    return gene_mutation_frequency(sub, top_n=top_n, title="Hereditary Breast Cancer Genes")


def summary_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """Generate overall statistics table."""
    if df.empty:
        return pd.DataFrame({"Metric": ["No data"], "Value": ["—"]})

    stats = []
    gene_col = "geneSymbol" if "geneSymbol" in df.columns else "hugoGeneSymbol"
    if gene_col not in df.columns and "entrezGeneId" in df.columns:
        gene_col = "entrezGeneId"
    sample_col = "sampleId" if "sampleId" in df.columns else "uniqueSampleKey"

    stats.append(("Total Mutations", len(df)))
    if sample_col in df.columns:
        stats.append(("Unique Samples", df[sample_col].nunique()))
    if gene_col in df.columns:
        stats.append(("Unique Genes", df[gene_col].nunique()))
    if "mutationType" in df.columns:
        stats.append(("Mutation Types", df["mutationType"].nunique()))

    return pd.DataFrame(stats, columns=["Metric", "Value"])


def gene_across_studies_bar(counts_df: pd.DataFrame, gene_symbol: str, top_n: int = 30) -> Tuple[plt.Figure, pd.DataFrame]:
    """Bar plot of gene mutation count across studies."""
    _setup_style()
    if counts_df.empty or "studyName" not in counts_df.columns or "count" not in counts_df.columns:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()
    df = counts_df.nlargest(top_n, "count").sort_values("count", ascending=True)
    fig, ax = plt.subplots(figsize=(12, max(6, len(df) * 0.35)))
    colors = sns.color_palette("viridis", len(df))[::-1]
    ax.barh(range(len(df)), df["count"].values, color=colors)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["studyName"].str[:50] + ("..." if df["studyName"].str.len().gt(50).any() else ""), fontsize=9)
    ax.set_xlabel("Number of Mutations")
    ax.set_title(f"{gene_symbol} Mutations Across Studies", fontweight="bold")
    ax.invert_yaxis()
    plt.tight_layout()
    return fig, df[["studyName", "studyId", "count"]].rename(columns={"count": "Mutation Count"})


def lollipop_mutations(df: pd.DataFrame, gene_symbol: str) -> Tuple[plt.Figure, pd.DataFrame]:
    """Lollipop plot: protein position vs mutation count/frequency."""
    _setup_style()
    if df.empty:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No mutation data", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    def _get_pos(row):
        if "proteinPosStart" in row.index and pd.notna(row.get("proteinPosStart")):
            try:
                return int(float(row["proteinPosStart"]))
            except (ValueError, TypeError):
                pass
        if "proteinPosition" in row.index and pd.notna(row.get("proteinPosition")):
            try:
                return int(float(row["proteinPosition"]))
            except (ValueError, TypeError):
                pass
        pc = str(row.get("proteinChange", "") or row.get("aminoAcidChange", ""))
        m = re.search(r"(\d+)", pc)
        return int(m.group(1)) if m else None

    df = df.copy()
    df["_pos"] = df.apply(_get_pos, axis=1)
    df = df.dropna(subset=["_pos"])
    if df.empty:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No protein position data for lollipop", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    pos_counts = df.groupby("_pos").size().reset_index(name="count")
    pos_counts = pos_counts.sort_values("_pos")
    x = pos_counts["_pos"].values
    y = pos_counts["count"].values

    fig, ax = plt.subplots(figsize=(14, 6))
    for xi, yi in zip(x, y):
        ax.plot([xi, xi], [0, yi], "C0-", linewidth=1.5)
    ax.scatter(x, y, s=50, c="C0", zorder=5, edgecolors="white", linewidths=1)
    ax.set_xlabel("Protein Position (amino acid)")
    ax.set_ylabel("Mutation Count")
    ax.set_title(f"{gene_symbol} Mutation Lollipop (all variations)", fontweight="bold")
    ax.set_ylim(bottom=0)
    plt.xticks(rotation=0)
    plt.tight_layout()

    stats_df = pos_counts.rename(columns={"_pos": "Protein Position", "count": "Mutation Count"})
    return fig, stats_df


def deg_downstream_chart(deg_df: pd.DataFrame, gene_symbol: str) -> Tuple[plt.Figure, pd.DataFrame]:
    """Bar chart of differential expression (median diff, -log10 p) for downstream genes."""
    _setup_style()
    if deg_df.empty or "Gene" not in deg_df.columns or "p_value" not in deg_df.columns:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No DEG data available", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    df = deg_df.copy()
    df["neglog10p"] = -np.log10(df["p_value"].replace(0, 1e-20))
    df["label"] = df["Gene"] + " (" + df["Comparison"] + ")"
    df = df.sort_values("neglog10p", ascending=True).tail(25)

    fig, ax = plt.subplots(figsize=(10, max(6, len(df) * 0.3)))
    dir_col = df["Direction"] if "Direction" in df.columns else pd.Series(["up"] * len(df))
    colors = ["#2ecc71" if str(d).lower() == "up" else "#e74c3c" for d in dir_col]
    ax.barh(range(len(df)), df["neglog10p"], color=colors)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["label"], fontsize=9)
    ax.set_xlabel("-log10(p-value)")
    ax.axvline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.7, label="p=0.05")
    ax.set_title(f"DEG downstream of {gene_symbol}: Wild vs LoF/GoF", fontweight="bold")
    ax.legend(loc="lower right")
    ax.invert_yaxis()
    plt.tight_layout()
    return fig, deg_df


def subtype_enrichment_chart(
    enrichment_results: list,
    gene_symbol: str,
) -> Tuple[plt.Figure, pd.DataFrame]:
    """Bar chart of mutation rate by molecular subtype, with chi-squared p-values."""
    _setup_style()
    if not enrichment_results:
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.text(0.5, 0.5, "No subtype data available", ha="center", va="center", fontsize=14)
        return fig, pd.DataFrame()

    all_stats = []
    n_studies = len(enrichment_results)
    fig, axes = plt.subplots(n_studies, 1, figsize=(10, 4 * n_studies), sharex=False)
    if n_studies == 1:
        axes = [axes]
    for ax, res in zip(axes, enrichment_results):
        tbl = res["subtype_df"].copy()
        tbl["studyId"] = res["studyId"]
        tbl["studyName"] = res["studyName"]
        tbl["p_value"] = res["p_value"]
        all_stats.append(tbl)
        p_str = f"p = {res['p_value']:.4f}" if not pd.isna(res["p_value"]) else "p = N/A"
        sig = " **" if res["p_value"] < 0.05 else ""
        ax.barh(tbl["Subtype"], tbl["Mutation rate (%)"], color=sns.color_palette(PALETTE, len(tbl)))
        ax.set_xlabel("Mutation rate (%)")
        ax.set_title(f"{res['studyName'][:55]} - {p_str}{sig}", fontsize=11)
        ax.invert_yaxis()
        ax.set_xlim(0, 100)
    plt.tight_layout()
    combined = pd.concat(all_stats, ignore_index=True)
    return fig, combined


def _format_survival_type(col_name: str) -> str:
    """Map OS_MONTHS, DFS_MONTHS etc. to display name."""
    m = {"OS_MONTHS": "Overall Survival", "DFS_MONTHS": "Disease-Free Survival",
         "PFS_MONTHS": "Progression-Free Survival", "RFS_MONTHS": "Recurrence-Free Survival",
         "DSS_MONTHS": "Disease-Specific Survival"}
    return m.get(col_name, col_name.replace("_", " ").title())


def _write_survival_interpretation(
    gene_symbol: str,
    survival_type: str,
    stats_df: pd.DataFrame,
    logrank_p: float,
    median_by_group: dict,
) -> str:
    """Generate plain-language interpretation of survival analysis."""
    lines = []
    lines.append(f"**Interpretation:** Kaplan-Meier curves compare {survival_type} among patients with *{gene_symbol}* **Wild Type** (no mutation), **Gain of Function** (missense/in-frame), and **Loss of Function** (nonsense/frameshift/splice) mutations.")
    lines.append("")
    if median_by_group:
        lines.append("**Median survival (months):**")
        for grp, med in median_by_group.items():
            med_str = f"{med:.1f}" if pd.notna(med) and med < 1e10 else "NR"
            lines.append(f"- {grp}: {med_str}")
        lines.append("")
    lines.append(f"**Log-rank test (across groups):** p = {logrank_p:.4f}" + (" (significant)" if logrank_p < 0.05 else " (not significant)"))
    lines.append("")
    if logrank_p < 0.05:
        lines.append("Groups differ significantly. Loss of Function mutations often associate with poorer outcomes; Gain of Function may show different biology. Consider clinical context and sample sizes.")
    else:
        lines.append("No statistically significant difference between groups in this cohort. Small numbers in mutation groups may limit power. Results are exploratory.")
    lines.append("")
    lines.append("*This analysis is for research only. Not for clinical decision-making.*")
    return "\n".join(lines)


def survival_plot_gof_lof(
    surv_df: pd.DataFrame,
    gene_symbol: str,
    survival_type: str = "Overall Survival",
) -> Tuple[plt.Figure, pd.DataFrame, str]:
    """Kaplan-Meier survival plot: GoF vs LoF vs Wild Type. Returns (fig, stats_df, interpretation)."""
    _setup_style()
    empty_fig, empty_ax = plt.subplots(figsize=FIG_SIZE)
    if surv_df.empty or "group" not in surv_df.columns or "time" not in surv_df.columns or "event" not in surv_df.columns:
        empty_ax.text(0.5, 0.5, "No survival data available", ha="center", va="center", fontsize=14)
        return empty_fig, pd.DataFrame(), ""

    try:
        from lifelines import KaplanMeierFitter
        from lifelines.statistics import multivariate_logrank_test
    except ImportError:
        empty_ax.text(0.5, 0.5, "Install lifelines: pip install lifelines", ha="center", va="center", fontsize=14)
        return empty_fig, pd.DataFrame(), ""

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {"Loss of Function": "#e74c3c", "Gain of Function": "#3498db", "Wild Type": "#2ecc71"}
    order = ["Wild Type", "Gain of Function", "Loss of Function"]
    kmf = KaplanMeierFitter()
    results = []
    median_by_group = {}

    for grp in order:
        if grp not in surv_df["group"].values:
            continue
        d = surv_df[surv_df["group"] == grp]
        T, E = d["time"], d["event"]
        kmf.fit(T, E, label=f"{grp} (n={len(d)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors.get(grp, "gray"))
        results.append({"Group": grp, "N": len(d), "Events": int(E.sum())})
        try:
            med = kmf.median_survival_time_
            median_by_group[grp] = med
        except Exception:
            median_by_group[grp] = None

    # Log-rank test across all groups
    try:
        res = multivariate_logrank_test(surv_df["time"], surv_df["group"], surv_df["event"])
        logrank_p = res.p_value
    except Exception:
        logrank_p = float("nan")

    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_title(f"{gene_symbol}: {survival_type} by Mutation Type", fontweight="bold")
    ax.legend(loc="lower left")
    ax.set_ylim(0, 1.05)
    plt.tight_layout()

    stats_df = pd.DataFrame(results) if results else pd.DataFrame()
    # Add median to stats_df if available
    if median_by_group and not stats_df.empty:
        def _fmt_med(g):
            v = median_by_group.get(g)
            if pd.isna(v) or v is None or v >= 1e10:
                return "NR"
            return f"{v:.1f}"
        stats_df["Median (months)"] = stats_df["Group"].map(_fmt_med)
    interpretation = _write_survival_interpretation(gene_symbol, survival_type, stats_df, logrank_p, median_by_group)
    return fig, stats_df, interpretation


ANALYSIS_TYPES = {
    "Mutation Type Distribution": mutation_type_distribution,
    "Top Mutated Genes": gene_mutation_frequency,
    "Hereditary Breast Cancer Genes": hereditary_genes_analysis,
    "Variant Classification": variant_classification_distribution,
    "Sample Mutation Burden": sample_mutation_burden,
    "Mutation Matrix (Gene × Sample)": oncoprint_style_matrix,
}
