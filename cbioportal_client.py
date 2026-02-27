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
    return _get(f"/studies/{study_id}")[0] if _get(f"/studies/{study_id}") else {}


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
) -> pd.DataFrame:
    """
    Fetch mutations for a molecular profile.
    If no sample_ids or entrez_gene_ids provided, fetches all (may be slow for large studies).
    """
    if not sample_ids and not entrez_gene_ids:
        # Fetch all mutations - API may require at least one filter
        # Try without filters - some endpoints accept empty
        try:
            data = _get(f"/molecular-profiles/{molecular_profile_id}/mutations")
            if not data:
                return pd.DataFrame()
            return pd.DataFrame(data)
        except Exception:
            return pd.DataFrame()

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


def fetch_mutations_by_study(
    study_id: str,
    molecular_profile_id: str,
    sample_ids: Optional[list] = None,
) -> pd.DataFrame:
    """Fetch mutations via molecular profile endpoint (returns data; /mutations/fetch often returns empty)."""
    samples = get_samples(study_id)
    if samples.empty:
        return pd.DataFrame()

    sample_ids_to_use = sample_ids or samples["sampleId"].tolist()
    if len(sample_ids_to_use) > 500:
        sample_ids_to_use = sample_ids_to_use[:500]

    # Primary: molecular profile endpoint (reliable)
    try:
        df = get_mutations(molecular_profile_id, sample_ids=sample_ids_to_use)
        if not df.empty:
            return df
    except Exception:
        pass

    # Fallback: generic /mutations/fetch (often returns empty)
    body = {
        "sampleIdentifiers": [
            {"studyId": study_id, "sampleId": sid} for sid in sample_ids_to_use
        ],
    }
    try:
        data = _post("/mutations/fetch", body)
        if data:
            return pd.DataFrame(data)
    except Exception:
        pass
    return pd.DataFrame()


def get_cancer_types() -> pd.DataFrame:
    """Fetch cancer types for filtering."""
    data = _get("/cancer-types")
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)
