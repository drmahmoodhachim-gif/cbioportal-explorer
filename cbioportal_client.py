"""
cBioPortal API Client for fetching cancer genomics data.
API docs: https://www.cbioportal.org/api/swagger-ui/index.html
"""

import requests
from typing import Optional
import pandas as pd

API_BASE = "https://www.cbioportal.org/api"


def _get(endpoint: str, params: Optional[dict] = None) -> list:
    """GET request to cBioPortal API."""
    url = f"{API_BASE}{endpoint}"
    r = requests.get(url, params=params or {}, timeout=60)
    r.raise_for_status()
    return r.json() if r.text else []


def _post(endpoint: str, json_data: dict) -> list:
    """POST request to cBioPortal API."""
    url = f"{API_BASE}{endpoint}"
    r = requests.post(url, json=json_data, timeout=120)
    r.raise_for_status()
    return r.json() if r.text else []


def get_all_studies(projection: str = "SUMMARY") -> pd.DataFrame:
    """Fetch all cancer studies."""
    data = _get("/studies", params={"projection": projection})
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)


def get_study(study_id: str) -> dict:
    """Get details for a specific study."""
    try:
        r = requests.get(f"{API_BASE}/studies/{study_id}", timeout=30)
        r.raise_for_status()
        d = r.json()
        return d if isinstance(d, dict) else (d[0] if d else {})
    except Exception:
        return {}


def get_molecular_profiles(study_id: str) -> pd.DataFrame:
    """Get molecular profiles (datasets) for a study."""
    data = _get(f"/studies/{study_id}/molecular-profiles", params={"projection": "DETAILED"})
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)


def get_sample_lists(study_id: str) -> pd.DataFrame:
    """Get sample lists for a study."""
    data = _get(f"/studies/{study_id}/sample-lists")
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)


def get_samples(study_id: str) -> pd.DataFrame:
    """Get samples for a study."""
    data = _get(f"/studies/{study_id}/samples", params={"projection": "SUMMARY"})
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)


def get_mutations(
    molecular_profile_id: str,
    sample_ids: Optional[list] = None,
    entrez_gene_ids: Optional[list] = None,
    study_id: Optional[str] = None,
    max_samples: int = 400,
) -> pd.DataFrame:
    """Fetch mutations for a molecular profile. Pass study_id to auto-fetch samples."""
    if not sample_ids and study_id:
        samples = get_samples(study_id)
        if not samples.empty:
            sample_ids = samples["sampleId"].tolist()
    if not sample_ids and not entrez_gene_ids:
        return pd.DataFrame()
    if sample_ids and len(sample_ids) > max_samples:
        sample_ids = sample_ids[:max_samples]
    body = {"sampleIds": sample_ids or []}
    if entrez_gene_ids:
        body["entrezGeneIds"] = entrez_gene_ids
    try:
        data = _post(f"/molecular-profiles/{molecular_profile_id}/mutations/fetch", body)
        if not data:
            return pd.DataFrame()
        return pd.DataFrame(data)
    except Exception as e:
        raise RuntimeError(f"Failed to fetch mutations: {e}") from e


def get_genes_by_entrez(entrez_ids: list) -> pd.DataFrame:
    """Fetch gene info by Entrez IDs. Returns hugoGeneSymbol."""
    if not entrez_ids:
        return pd.DataFrame()
    gene_ids = [str(int(x)) for x in entrez_ids if pd.notna(x)]
    if not gene_ids:
        return pd.DataFrame()
    try:
        data = requests.post(
            f"{API_BASE}/genes/fetch",
            json=gene_ids[:500],
            params={"geneIdType": "ENTREZ_GENE_ID"},
            timeout=30,
        ).json()
        return pd.DataFrame(data) if data else pd.DataFrame()
    except Exception:
        return pd.DataFrame()


def add_gene_symbols(mutations_df: pd.DataFrame) -> pd.DataFrame:
    """Add geneSymbol column by looking up entrezGeneId via cBioPortal API."""
    if mutations_df.empty or "entrezGeneId" not in mutations_df.columns:
        return mutations_df
    entrez_ids = mutations_df["entrezGeneId"].dropna().unique().astype(int).tolist()
    if not entrez_ids:
        return mutations_df
    genes_df = get_genes_by_entrez(entrez_ids)
    if genes_df.empty or "hugoGeneSymbol" not in genes_df.columns:
        mutations_df = mutations_df.copy()
        mutations_df["geneSymbol"] = mutations_df["entrezGeneId"].astype(str)
        return mutations_df
    gene_map = genes_df.set_index("entrezGeneId")["hugoGeneSymbol"].to_dict()
    mutations_df = mutations_df.copy()
    mutations_df["geneSymbol"] = mutations_df["entrezGeneId"].map(
        lambda x: gene_map.get(int(x), str(x))
    )
    return mutations_df


def get_cancer_types() -> pd.DataFrame:
    """Fetch cancer types for filtering."""
    data = _get("/cancer-types")
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)
