"""
cBioPortal API Client for fetching cancer genomics data.
API docs: https://www.cbioportal.org/api/swagger-ui/index.html
"""

import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Tuple
import pandas as pd

API_BASE = "https://www.cbioportal.org/api"

HEREDITARY_BREAST_CANCER_GENES = {
    "BRCA1": 672, "BRCA2": 675, "PALB2": 79728, "TP53": 7157, "PTEN": 5728,
    "CDH1": 999, "STK11": 6794, "ATM": 472, "CHEK2": 11200, "BARD1": 580,
    "RAD51C": 5889, "RAD51D": 5892, "NF1": 4763,
}


def _add_gene_symbols(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "hugoGeneSymbol" in df.columns or "geneSymbol" in df.columns:
        return df
    if "entrezGeneId" not in df.columns:
        return df
    id_to_symbol = {v: k for k, v in HEREDITARY_BREAST_CANCER_GENES.items()}

    def _to_symbol(x):
        try:
            return id_to_symbol.get(int(float(x)), str(int(float(x))))
        except (ValueError, TypeError):
            return str(x) if pd.notna(x) else ""

    out = df.copy()
    out["hugoGeneSymbol"] = df["entrezGeneId"].apply(_to_symbol)
    return out


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
    """Fetch mutations via molecular profile endpoint (reliable). /mutations/fetch often returns empty."""
    samples = get_samples(study_id)
    if samples.empty:
        return pd.DataFrame()
    sample_ids_to_use = sample_ids or samples["sampleId"].tolist()
    if len(sample_ids_to_use) > 500:
        sample_ids_to_use = sample_ids_to_use[:500]
    try:
        df = get_mutations(molecular_profile_id, sample_ids=sample_ids_to_use)
        if not df.empty:
            return _add_gene_symbols(df)
    except Exception:
        pass
    body = {"sampleIdentifiers": [{"studyId": study_id, "sampleId": sid} for sid in sample_ids_to_use]}
    try:
        data = _post("/mutations/fetch", body)
        if data:
            return _add_gene_symbols(pd.DataFrame(data))
    except Exception:
        pass
    return pd.DataFrame()


def get_entrez_id(gene_symbol: str) -> Optional[int]:
    """Resolve Hugo gene symbol to Entrez Gene ID via cBioPortal API."""
    if not gene_symbol or not str(gene_symbol).strip():
        return None
    sym = str(gene_symbol).strip().upper()
    try:
        data = _get("/genes", params={"keyword": sym, "pageSize": 10})
        for g in data or []:
            if g.get("hugoGeneSymbol", "").upper() == sym:
                return int(g.get("entrezGeneId", 0)) or None
        if data and len(data) > 0:
            return int(data[0].get("entrezGeneId", 0)) or None
    except Exception:
        pass
    return None


def fetch_gene_across_studies(
    gene_symbol: str,
    studies_df: pd.DataFrame,
    max_studies: int = 50,
    samples_per_study: int = 200,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fetch mutations for a gene across breast cancer studies.
    Returns (mutations_df with studyId, study_counts_df for bar plot).
    """
    entrez = get_entrez_id(gene_symbol)
    if not entrez:
        return pd.DataFrame(), pd.DataFrame()

    all_muts = []
    study_counts = []

    for _, row in studies_df.head(max_studies).iterrows():
        sid = row.get("studyId")
        if not sid:
            continue
        profiles = get_molecular_profiles(sid)
        if profiles.empty or "molecularAlterationType" not in profiles.columns:
            continue
        mut_profiles = profiles[
            profiles["molecularAlterationType"].str.upper().str.contains("MUTATION", na=False)
        ]
        if mut_profiles.empty:
            mut_profiles = profiles
        profile_id = mut_profiles.iloc[0].get("molecularProfileId")
        if not profile_id:
            continue
        samples = get_samples(sid)
        if samples.empty:
            continue
        sample_ids = samples["sampleId"].head(samples_per_study).tolist()
        try:
            df = get_mutations(profile_id, sample_ids=sample_ids, entrez_gene_ids=[entrez])
            if not df.empty:
                df = _add_gene_symbols(df)
                df["studyId"] = sid
                df["studyName"] = row.get("name", sid)
                all_muts.append(df)
                study_counts.append({"studyId": sid, "studyName": row.get("name", sid), "count": len(df)})
        except Exception:
            continue

    if not all_muts:
        return pd.DataFrame(), pd.DataFrame()
    muts_df = pd.concat(all_muts, ignore_index=True)
    counts_df = pd.DataFrame(study_counts)
    return muts_df, counts_df


_SUBTYPE_ATTRS = [
    "SUBTYPE", "CLAUDIN_SUBTYPE", "PAM50_MRNA", "BREAST_SUBTYPE",
    "IHC_SUBTYPE", "INTCLUST", "THREEGENE",
]


def fetch_subtype_enrichment(
    muts_df: pd.DataFrame,
    studies_df: pd.DataFrame,
    top_n_studies: int = 5,
) -> list:
    """
    For studies with mutations, fetch subtype data and compute chi-squared enrichment.
    Returns list of {studyId, studyName, subtype_stats_df, p_value, has_subtype}.
    """
    if muts_df.empty or "studyId" not in muts_df.columns or "sampleId" not in muts_df.columns:
        return []
    mutated = muts_df.groupby("studyId")["sampleId"].apply(set).to_dict()
    top_studies = muts_df["studyId"].value_counts().head(top_n_studies)
    results = []
    for sid in top_studies.index:
        study_row = studies_df[studies_df["studyId"] == sid].iloc[0] if sid in studies_df["studyId"].values else None
        study_name = study_row["name"] if study_row is not None else sid
        try:
            samples = get_samples(sid)
            if samples.empty or "patientId" not in samples.columns:
                continue
            clinical = get_clinical_data(sid, clinical_data_type="PATIENT")
            if clinical.empty:
                continue
            piv = clinical.pivot_table(index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first")
            subtype_col = None
            for col in _SUBTYPE_ATTRS:
                if col in piv.columns:
                    subtype_col = col
                    break
            if subtype_col is None:
                continue
            sub_df = piv[[subtype_col]].reset_index()
            sub_df = sub_df.rename(columns={subtype_col: "subtype"})
            samples = samples.merge(sub_df, on="patientId", how="left")
            samples["subtype"] = samples["subtype"].fillna("Unknown").astype(str).str.strip()
            samples = samples[samples["subtype"] != "Unknown"]
            if samples.empty or samples["subtype"].nunique() < 2:
                continue
            mutated_sids = mutated.get(sid, set())
            samples["mutated"] = samples["sampleId"].isin(mutated_sids).astype(int)
            tbl = samples.groupby("subtype").agg({"mutated": ["sum", "count"]})
            tbl.columns = ["mutated", "total"]
            tbl["not_mutated"] = tbl["total"] - tbl["mutated"]
            tbl = tbl[tbl["total"] >= 5]
            if tbl.shape[0] < 2 or tbl["mutated"].sum() < 2:
                continue
            try:
                from scipy.stats import chi2_contingency
                cont = tbl[["mutated", "not_mutated"]].values
                chi2, p_val, _, _ = chi2_contingency(cont)
            except Exception:
                p_val = float("nan")
            tbl["rate_pct"] = (tbl["mutated"] / tbl["total"] * 100).round(1)
            tbl = tbl.reset_index()
            tbl.columns = ["Subtype", "N mutated", "N total", "N not mutated", "Mutation rate (%)"]
            results.append({
                "studyId": sid,
                "studyName": study_name,
                "subtype_df": tbl,
                "p_value": p_val,
                "has_subtype": True,
            })
        except Exception:
            continue
    return results


# Survival attribute pairs: (time_col, event_col) - both must be patient-level
_SURVIVAL_ATTR_PAIRS = [
    ("OS_MONTHS", "OS_STATUS"),
    ("DFS_MONTHS", "DFS_STATUS"),
    ("PFS_MONTHS", "PFS_STATUS"),
    ("RFS_MONTHS", "RFS_STATUS"),
    ("DSS_MONTHS", "DSS_STATUS"),
]


def study_has_survival_data(study_id: str) -> bool:
    """Check if a study has patient-level survival attributes (OS, DFS, PFS, RFS, or DSS)."""
    try:
        data = _get(f"/studies/{study_id}/clinical-attributes", params={"projection": "SUMMARY"})
        if not data:
            return False
        attr_ids = {a.get("clinicalAttributeId") for a in data if a.get("patientAttribute")}
        for time_id, event_id in _SURVIVAL_ATTR_PAIRS:
            if time_id in attr_ids and event_id in attr_ids:
                return True
        return False
    except Exception:
        return False


def filter_studies_with_survival(studies_df: pd.DataFrame) -> pd.DataFrame:
    """Return only studies that have patient-level survival data."""
    if studies_df.empty:
        return studies_df
    study_ids = studies_df["studyId"].tolist()
    has_survival = {}
    with ThreadPoolExecutor(max_workers=10) as ex:
        futures = {ex.submit(study_has_survival_data, sid): sid for sid in study_ids}
        for f in as_completed(futures):
            sid = futures[f]
            try:
                has_survival[sid] = f.result()
            except Exception:
                has_survival[sid] = False
    mask = studies_df["studyId"].map(lambda x: has_survival.get(x, False))
    return studies_df[mask].reset_index(drop=True)


def get_clinical_data(study_id: str, clinical_data_type: str = "SAMPLE") -> pd.DataFrame:
    """Fetch clinical data for a study.
    clinical_data_type: 'PATIENT' for survival (OS_MONTHS, OS_STATUS, etc.), 'SAMPLE' for sample-level.
    """
    try:
        params = {"projection": "SUMMARY", "pageSize": 100000}
        if clinical_data_type and clinical_data_type.upper() in ("PATIENT", "SAMPLE"):
            params["clinicalDataType"] = clinical_data_type.upper()
        data = _get(f"/studies/{study_id}/clinical-data", params=params)
        if not data:
            return pd.DataFrame()
        df = pd.DataFrame(data)
        return df
    except Exception:
        return pd.DataFrame()


def fetch_survival_data_for_gene(
    study_id: str,
    molecular_profile_id: str,
    gene_symbol: str,
    sample_ids: Optional[list] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, str, str]:
    """
    Fetch mutations + clinical data for survival analysis by gene.
    Returns (survival_df, group_counts, time_col_name, error_msg).
    group = 'Loss of Function' | 'Gain of Function' | 'Wild Type'
    """
    entrez = get_entrez_id(gene_symbol)
    if not entrez:
        return pd.DataFrame(), pd.DataFrame(), "", f"Gene {gene_symbol} not found"

    samples = get_samples(study_id)
    if samples.empty:
        return pd.DataFrame(), pd.DataFrame(), "", "No samples"
    sample_ids_to_use = sample_ids or samples["sampleId"].tolist()
    if len(sample_ids_to_use) > 500:
        sample_ids_to_use = sample_ids_to_use[:500]

    muts_df = pd.DataFrame()
    try:
        muts_df = get_mutations(molecular_profile_id, sample_ids=sample_ids_to_use, entrez_gene_ids=[entrez])
        if not muts_df.empty:
            muts_df = _add_gene_symbols(muts_df)
    except Exception:
        pass

    # Survival attributes (OS_MONTHS, OS_STATUS, etc.) are patient-level; use clinicalDataType=PATIENT
    clinical = get_clinical_data(study_id, clinical_data_type="PATIENT")
    if clinical.empty:
        return pd.DataFrame(), pd.DataFrame(), "", "No clinical data"

    attr_ids = set(clinical["clinicalAttributeId"].unique()) if "clinicalAttributeId" in clinical.columns else set()
    time_col = None
    event_col = None
    for t in ["OS_MONTHS", "DFS_MONTHS", "PFS_MONTHS", "RFS_MONTHS", "DAYS_TO_LAST_FOLLOWUP", "DSS_MONTHS"]:
        if t in attr_ids:
            time_col = t
            break
    for e in ["OS_STATUS", "DFS_STATUS", "PFS_STATUS", "RFS_STATUS", "DSS_STATUS"]:
        if e in attr_ids:
            event_col = e
            break
    if not time_col or not event_col:
        return pd.DataFrame(), pd.DataFrame(), "", "No survival data (looked for OS/DFS/PFS/RFS)"

    idx = "patientId" if "patientId" in clinical.columns else "sampleId"
    piv = clinical.pivot_table(
        index=idx,
        columns="clinicalAttributeId",
        values="value",
        aggfunc="first",
    ).reset_index()
    id_col = idx
    if time_col not in piv.columns or event_col not in piv.columns:
        return pd.DataFrame(), pd.DataFrame(), "", f"Missing {time_col} or {event_col}"

    sample_to_patient = samples.set_index("sampleId")["patientId"].to_dict() if "patientId" in samples.columns else {}
    if not sample_to_patient:
        sample_to_patient = {s: s for s in sample_ids_to_use}

    LOF_TYPES = {
        "frame_shift_del", "frame_shift_ins", "nonsense_mutation", "splice_site",
        "nonstop_mutation", "de_novo_start_outofframe", "start_codon_snp",
    }
    GOF_TYPES = {"missense_mutation", "in_frame_del", "in_frame_ins", "missense_variant"}

    def _classify(mut_type):
        mt = str(mut_type).lower().replace(" ", "_") if pd.notna(mut_type) else ""
        if mt in LOF_TYPES:
            return "Loss of Function"
        if mt in GOF_TYPES:
            return "Gain of Function"
        return "Other"

    patient_groups = {}
    if not muts_df.empty and "sampleId" in muts_df.columns:
        mut_col = "mutationType" if "mutationType" in muts_df.columns else "variantType"
        type_col = mut_col if mut_col in muts_df.columns else None
        for _, row in muts_df.iterrows():
            pid = sample_to_patient.get(row["sampleId"], row["sampleId"])
            grp = _classify(row.get(type_col, "")) if type_col else "Gain of Function"
            if pid not in patient_groups or grp == "Loss of Function":
                patient_groups[pid] = grp
        for pid in patient_groups:
            if patient_groups[pid] == "Other":
                patient_groups[pid] = "Gain of Function"
    for sid in sample_ids_to_use:
        pid = sample_to_patient.get(sid, sid)
        if pid not in patient_groups:
            patient_groups[pid] = "Wild Type"

    piv["group"] = piv[id_col].map(patient_groups).fillna("Wild Type")
    piv = piv[piv["group"].isin(["Loss of Function", "Gain of Function", "Wild Type"])]
    if piv.empty or piv["group"].nunique() < 2:
        return pd.DataFrame(), pd.DataFrame(), "", "Need at least 2 groups for survival comparison"

    piv[time_col] = pd.to_numeric(piv[time_col], errors="coerce")
    piv = piv.dropna(subset=[time_col])
    piv["event"] = piv[event_col].astype(str).str.upper().str.contains("DECEASED|RECURRED|PROGRESSED|1").astype(int)

    surv_df = piv[[id_col, "group", time_col, "event"]].copy()
    surv_df = surv_df.rename(columns={time_col: "time"})
    # Add molecular subtype if available - try common breast cancer subtype attributes
    subtype_cols = [
        "SUBTYPE",           # TCGA PanCan (BRCA_LumA, BRCA_Basal, etc.)
        "CLAUDIN_SUBTYPE",   # METABRIC (LumA, LumB, Her2, Basal, Claudin-low)
        "PAM50_MRNA", "BREAST_SUBTYPE", "IHC_SUBTYPE",
        "INTCLUST",          # METABRIC integrative clusters
        "THREEGENE",         # 3-gene classifier
    ]
    for col in subtype_cols:
        if col in piv.columns:
            surv_df["subtype"] = piv[col].fillna("Unknown").astype(str).str.strip()
            break
    else:
        surv_df["subtype"] = "All"
    if "DAYS" in time_col:
        surv_df["time"] = surv_df["time"] / 30.44
    counts = surv_df["group"].value_counts().reset_index()
    counts.columns = ["Group", "N"]
    return surv_df, counts, time_col, ""

def get_cancer_types() -> pd.DataFrame:
    """Fetch cancer types for filtering."""
    data = _get("/cancer-types")
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)
