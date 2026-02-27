"""
KEGG pathway enrichment from DEG gene list.
Uses KEGG REST API (rest.kegg.jp) for gene-pathway mappings.
"""

import re
import time
from typing import Callable, List, Optional, Tuple

import pandas as pd
import requests

KEGG_API = "https://rest.kegg.jp"
_RATE_LIMIT = 0.4  # ~2.5 req/sec to stay under KEGG limit


def _kegg_get(path: str) -> str:
    """GET from KEGG REST API with rate limiting."""
    time.sleep(_RATE_LIMIT)
    r = requests.get(f"{KEGG_API}{path}", timeout=30)
    r.raise_for_status()
    return r.text


_PATHWAY_GENES: Optional[dict] = None
_PATHWAY_NAMES: Optional[dict] = None


def _load_kegg_data() -> Tuple[dict, dict]:
    """Load and cache KEGG pathway data."""
    global _PATHWAY_GENES, _PATHWAY_NAMES
    if _PATHWAY_GENES is not None and _PATHWAY_NAMES is not None:
        return _PATHWAY_GENES, _PATHWAY_NAMES
    # pathway -> set of entrez IDs
    pathway_genes = {}
    text = _kegg_get("/link/hsa/pathway")
    for line in text.strip().split("\n"):
        if not line or "\t" not in line:
            continue
        path_part, gene_part = line.split("\t", 1)
        m = re.match(r"path:(hsa\d+)", path_part.strip())
        g = re.match(r"hsa:(\d+)", gene_part.strip())
        if m and g:
            pid = m.group(1)
            entrez = int(g.group(1))
            pathway_genes.setdefault(pid, set()).add(entrez)
    # pathway ID -> name
    pathway_names = {}
    text = _kegg_get("/list/pathway/hsa")
    for line in text.strip().split("\n"):
        if not line or "\t" not in line:
            continue
        pid, name = line.split("\t", 1)
        pathway_names[pid] = name.strip()
    _PATHWAY_GENES = pathway_genes
    _PATHWAY_NAMES = pathway_names
    return pathway_genes, pathway_names


def fetch_kegg_enrichment(
    deg_genes: List[str],
    get_entrez_fn: Callable[[str], Optional[int]],
    background_size: int = 20000,
) -> Tuple[pd.DataFrame, str]:
    """
    Hypergeometric enrichment of DEG genes in KEGG pathways.
    deg_genes: list of gene symbols (Hugo)
    get_entrez_fn: function(symbol) -> entrez_id or None
    Returns (enrichment_df, error_msg). Columns: Pathway_ID, Pathway_Name, p_value, FDR, n_genes, n_DEG, n_overlap, Genes
    """
    if not deg_genes:
        return pd.DataFrame(), "No genes provided"
    try:
        pathway_genes, pathway_names = _load_kegg_data()
    except Exception as e:
        return pd.DataFrame(), f"KEGG API error: {e}"
    deg_entrez = set()
    entrez_to_symbol = {}
    for g in deg_genes:
        g = g.strip()
        e = get_entrez_fn(g)
        if e:
            deg_entrez.add(e)
            entrez_to_symbol[e] = g
    if not deg_entrez:
        return pd.DataFrame(), "Could not resolve DEG genes to Entrez IDs"
    K = len(deg_entrez)
    N = background_size
    from scipy.stats import hypergeom
    results = []
    for pid, p_genes in pathway_genes.items():
        n = len(p_genes)
        overlap = deg_entrez & p_genes
        k = len(overlap)
        if k < 2:
            continue
        try:
            p_val = hypergeom.sf(k - 1, N, n, K)
        except Exception:
            continue
        name = pathway_names.get(pid, pid)
        overlap_symbols = [entrez_to_symbol[e] for e in overlap if e in entrez_to_symbol]
        results.append({
            "Pathway_ID": pid,
            "Pathway_Name": name,
            "p_value": p_val,
            "n_pathway": n,
            "n_DEG": K,
            "n_overlap": k,
            "Genes": ", ".join(sorted(overlap_symbols)[:15]) + ("..." if len(overlap_symbols) > 15 else ""),
        })
    if not results:
        return pd.DataFrame(), "No pathway overlaps"
    df = pd.DataFrame(results)
    df = df[df["p_value"] <= 0.05]
    if df.empty:
        return pd.DataFrame(), "No pathways significant at p ≤ 0.05"
    # FDR (Benjamini-Hochberg)
    pvals = df["p_value"].values
    n = len(pvals)
    ranks = pvals.argsort().argsort() + 1
    fdr = (pvals * n / ranks).clip(max=1)
    df["FDR"] = fdr
    df = df.sort_values("p_value")
    df = df[["Pathway_ID", "Pathway_Name", "p_value", "FDR", "n_overlap", "n_pathway", "n_DEG", "Genes"]]
    return df, ""
