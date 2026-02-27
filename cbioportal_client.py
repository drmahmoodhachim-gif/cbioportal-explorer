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
# Cache for API-resolved Entrez ID -> symbol (so all genes appear as symbols, not numbers)
_ENTREZ_TO_SYMBOL_CACHE: dict = {v: k for k, v in HEREDITARY_BREAST_CANCER_GENES.items()}


def _fetch_one_gene(entrez_id: int) -> None:
    """Fetch one gene symbol via API; update cache. Used for parallel batch."""
    if entrez_id in _ENTREZ_TO_SYMBOL_CACHE:
        return
    try:
        r = requests.get(f"{API_BASE}/genes/{entrez_id}", timeout=5)
        if r.ok and r.text:
            data = r.json()
            sym = data.get("hugoGeneSymbol") or data.get("geneSymbol")
            if sym:
                _ENTREZ_TO_SYMBOL_CACHE[entrez_id] = str(sym)
                return
    except Exception:
        pass
    _ENTREZ_TO_SYMBOL_CACHE[entrez_id] = str(entrez_id)


def _add_gene_symbols(df: pd.DataFrame, max_api_lookups: int = 75) -> pd.DataFrame:
    """Add hugoGeneSymbol. Resolve unknown IDs via API in parallel; cap lookups to avoid slowness."""
    if df.empty or "hugoGeneSymbol" in df.columns or "geneSymbol" in df.columns:
        return df
    if "entrezGeneId" not in df.columns:
        return df
    # Prioritize most frequent genes (top mutated) for API resolution
    id_counts = df["entrezGeneId"].dropna().astype(str).value_counts()
    ids_to_fetch = []
    for uid in id_counts.index:
        try:
            eid = int(float(uid))
            if eid not in _ENTREZ_TO_SYMBOL_CACHE and len(ids_to_fetch) < max_api_lookups:
                ids_to_fetch.append(eid)
        except (ValueError, TypeError):
            pass
    if ids_to_fetch:
        with ThreadPoolExecutor(max_workers=10) as ex:
            list(ex.map(_fetch_one_gene, ids_to_fetch))
    # Use cache (IDs not resolved stay as numbers)
    def _to_symbol(x):
        try:
            return _ENTREZ_TO_SYMBOL_CACHE.get(int(float(x)), str(int(float(x))))
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
    sample_list_id: Optional[str] = None,
) -> pd.DataFrame:
    """
    Fetch mutations for a molecular profile.
    Prefer sample_list_id (queries ALL samples in list) over sample_ids (limited).
    """
    # Use GET with sampleListId to query all samples - matches cBioPortal's full query
    if sample_list_id and entrez_gene_ids:
        try:
            all_data = []
            page = 0
            page_size = 10000
            while True:
                data = _get(
                    f"/molecular-profiles/{molecular_profile_id}/mutations",
                    params={
                        "sampleListId": sample_list_id,
                        "entrezGeneId": entrez_gene_ids[0],
                        "projection": "SUMMARY",
                        "pageSize": page_size,
                        "pageNumber": page,
                    },
                )
                if not data:
                    break
                all_data.extend(data)
                if len(data) < page_size:
                    break
                page += 1
            if all_data:
                return pd.DataFrame(all_data)
        except Exception:
            pass

    # Fallback: POST fetch with sample IDs (limited to provided list)
    if not sample_ids and sample_list_id:
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
    max_samples: int = 500,
) -> pd.DataFrame:
    """Fetch mutations via molecular profile endpoint (reliable). /mutations/fetch often returns empty."""
    samples = get_samples(study_id)
    if samples.empty:
        return pd.DataFrame()
    sample_ids_to_use = sample_ids or samples["sampleId"].tolist()
    if len(sample_ids_to_use) > max_samples:
        sample_ids_to_use = sample_ids_to_use[:max_samples]
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


def _pick_sample_list_for_mutations(sample_lists_df: pd.DataFrame, study_id: str) -> Optional[str]:
    """Pick sample list that covers all/cases with mutation data. Prefer sequenced > cnaseq > all."""
    if sample_lists_df.empty or "sampleListId" not in sample_lists_df.columns:
        return None
    sl = sample_lists_df
    # Prefer: sequenced (mutation data), cnaseq (mutation+cna), all
    for pat in ["sequenced", "cnaseq", "mutation", "all"]:
        m = sl["sampleListId"].str.lower().str.contains(pat, na=False)
        if m.any():
            cand = sl[m].iloc[0]["sampleListId"]
            return str(cand)
    return sl.iloc[0]["sampleListId"] if len(sl) > 0 else None


def _pick_sample_list_for_expression(sample_lists_df: pd.DataFrame) -> Optional[str]:
    """Pick sample list likely to have expression data. Prefer rna/mrna/3way_complete > sequenced."""
    if sample_lists_df.empty or "sampleListId" not in sample_lists_df.columns:
        return None
    sl = sample_lists_df
    for pat in ["rna_seq", "mrna", "rna", "3way_complete", "complete", "sequenced", "cnaseq"]:
        m = sl["sampleListId"].str.lower().str.contains(pat, na=False)
        if m.any():
            return str(sl[m].iloc[0]["sampleListId"])
    return str(sl.iloc[0]["sampleListId"]) if len(sl) > 0 else None


def fetch_gene_across_studies(
    gene_symbol: str,
    studies_df: pd.DataFrame,
    max_studies: int = 50,
    samples_per_study: int = 2000,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fetch mutations for a gene across breast cancer studies.
    Uses sample lists to query ALL samples per study (matches cBioPortal).
    Fallback: sample_ids if no suitable list (limited to samples_per_study).
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
        sample_list_id = None
        try:
            sample_lists = get_sample_lists(sid)
            sample_list_id = _pick_sample_list_for_mutations(sample_lists, sid)
        except Exception:
            pass
        try:
            if sample_list_id:
                df = get_mutations(
                    profile_id,
                    sample_list_id=sample_list_id,
                    entrez_gene_ids=[entrez],
                )
            else:
                samples = get_samples(sid)
                if samples.empty:
                    continue
                sample_ids = samples["sampleId"].head(samples_per_study).tolist()
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


def study_has_expression_data(study_id: str) -> bool:
    """Check if a study has mRNA/expression molecular profile."""
    try:
        profiles = get_molecular_profiles(study_id)
        if profiles.empty or "molecularAlterationType" not in profiles.columns:
            return False
        expr = profiles["molecularAlterationType"].astype(str).str.upper().str.contains("MRNA|EXPRESSION", na=False)
        if expr.any():
            return True
        if "datatype" in profiles.columns:
            dt = profiles["datatype"].astype(str).str.upper().str.contains("Z-SCORE|CONTINUOUS", na=False)
            if dt.any():
                return True
        return False
    except Exception:
        return False


def filter_studies_with_expression(studies_df: pd.DataFrame) -> pd.DataFrame:
    """Return only studies that have mRNA/expression data."""
    if studies_df.empty:
        return studies_df
    study_ids = studies_df["studyId"].tolist()
    has_expr = {}
    with ThreadPoolExecutor(max_workers=10) as ex:
        futures = {ex.submit(study_has_expression_data, sid): sid for sid in study_ids}
        for f in as_completed(futures):
            sid = futures[f]
            try:
                has_expr[sid] = f.result()
            except Exception:
                has_expr[sid] = False
    mask = studies_df["studyId"].map(lambda x: has_expr.get(x, False))
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
    max_samples: int = 500,
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
    if len(sample_ids_to_use) > max_samples:
        sample_ids_to_use = sample_ids_to_use[:max_samples]

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

# Known downstream/pathway genes for DEG (gene-specific; from pathway DBs, literature)
DOWNSTREAM_GENES = {
    "BRCA1": ["RAD51", "BRCA2", "ATM", "CHEK2", "TP53", "CCND1", "MYC", "E2F1", "BARD1", "PALB2", "CDKN1A", "FANCD2", "MRE11", "RAD50", "NBN"],
    "BRCA2": ["RAD51", "BRCA1", "PALB2", "ATM", "TP53", "CDKN1A", "FANCD2", "FANCA", "RAD51C", "RAD51D"],
    "TP53": ["CDKN1A", "MDM2", "BAX", "BBC3", "CCND1", "CDK4", "E2F1", "BRCA1"],
    "PTEN": ["AKT1", "PIK3CA", "MTOR", "GSK3B", "FOXO1", "CCND1", "CDKN1A"],
    "PIK3CA": ["AKT1", "MTOR", "GSK3B", "PTEN", "FOXO1", "CCND1", "RPS6KB1"],
    "ERBB2": ["GRB2", "SOS1", "MAP2K1", "MAPK1", "PIK3CA", "AKT1", "MTOR", "MYC", "CCND1", "ERBB3", "EGFR"],
    "ESR1": ["CCND1", "MYC", "GREB1", "TFF1", "PGR", "BCL2", "CDKN1A", "FOXA1", "GATA3"],
    "CDH1": ["CTNNA1", "CTNND1", "JUP", "CDH2", "EGFR", "WNT3A"],
    "GATA3": ["ESR1", "FOXA1", "KRT18", "KRT19", "TFF1", "PGR"],
    "MAP3K1": ["MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "JUN", "ELK1"],
    "NF1": ["KRAS", "BRAF", "MAP2K1", "MAPK1", "RASA1", "RAF1"],
    "AKT1": ["MTOR", "GSK3B", "FOXO1", "BAD", "TSC2", "MDM2", "CCND1"],
    "ATM": ["TP53", "CHEK2", "BRCA1", "BRCA2", "RAD50", "MRE11", "NBN", "CDKN1A"],
    "CHEK2": ["TP53", "CDC25A", "BRCA1", "CDKN1A", "E2F1"],
    "PALB2": ["BRCA1", "BRCA2", "RAD51", "RAD51C", "FANCD2"],
    "CDKN1A": ["CDK2", "CDK4", "CDK6", "CCND1", "E2F1", "TP53"],
    "MYC": ["CCND1", "CDKN1A", "BCL2", "E2F1", "CDK4"],
}
# Fallback for genes without curated list: shared pathway / cell-cycle genes
DEFAULT_DOWNSTREAM = ["TP53", "CDKN1A", "CCND1", "MYC", "E2F1", "AKT1", "MTOR", "BRCA1", "BRCA2", "PTEN", "EGFR", "MAPK1"]

# Broad gene set for standalone DEG analysis (union of pathway genes + cancer-relevant)
DEG_GENES_FULL = list(dict.fromkeys(
    [g for genes in DOWNSTREAM_GENES.values() for g in genes]
    + DEFAULT_DOWNSTREAM
    + ["GRB2", "SOS1", "RAF1", "BRAF", "KRAS", "NRAS", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2",
       "BCL2", "BCL2L1", "BAD", "BAX", "BBC3", "MDM2", "CDK4", "CDK6", "CDK2", "RB1",
       "FANCA", "FANCD2", "FANCC", "ERBB3", "EGFR", "GREB1", "TFF1", "PGR", "FOXA1", "GATA3",
       "KRT18", "KRT19", "CTNNA1", "CTNND1", "JUP", "CDH2", "WNT3A", "GSK3B", "FOXO1",
       "RPS6KB1", "TSC2", "RASA1", "ELK1", "JUN", "E2F1", "CCND1", "CDKN1A", "CDKN2A",
       "RAD51", "RAD51C", "RAD51D", "MRE11", "RAD50", "NBN", "BARD1"]
))


def get_molecular_data(
    molecular_profile_id: str,
    sample_list_id: str,
    entrez_gene_ids: list,
) -> pd.DataFrame:
    """Fetch molecular (e.g. expression) data for genes."""
    try:
        body = {"sampleListId": sample_list_id, "entrezGeneIds": entrez_gene_ids}
        data = _post(f"/molecular-profiles/{molecular_profile_id}/molecular-data/fetch", body)
        if not data:
            return pd.DataFrame()
        return pd.DataFrame(data)
    except Exception:
        return pd.DataFrame()


def fetch_deg_downstream(
    study_id: str,
    mut_profile_id: str,
    gene_symbol: str,
    sample_list_id: Optional[str] = None,
) -> Tuple[pd.DataFrame, str]:
    """
    Differential expression of downstream genes: Wild vs LoF vs GoF.
    Returns (deg_df, error_msg). deg_df has Gene, Comparison, Median_diff, p_value, Direction.
    Uses sample list with expression data (rna_seq/mrna) when sample_list_id not provided.
    """
    entrez = get_entrez_id(gene_symbol)
    if not entrez:
        return pd.DataFrame(), f"Gene {gene_symbol} not found"

    if not sample_list_id:
        sl_df = get_sample_lists(study_id)
        sample_list_id = _pick_sample_list_for_expression(sl_df)
        if not sample_list_id:
            return pd.DataFrame(), "No suitable sample list (need expression data)"

    profiles = get_molecular_profiles(study_id)
    if profiles.empty:
        return pd.DataFrame(), "No molecular profiles"
    expr_profiles = profiles[
        profiles["molecularAlterationType"].astype(str).str.upper().str.contains("MRNA|EXPRESSION", na=False)
    ]
    if expr_profiles.empty and "datatype" in profiles.columns:
        expr_profiles = profiles[profiles["datatype"].astype(str).str.upper().str.contains("Z-SCORE|CONTINUOUS", na=False)]
    if expr_profiles.empty:
        return pd.DataFrame(), "No expression profile in study"

    expr_profile_id = expr_profiles.iloc[0]["molecularProfileId"]
    samples = get_samples(study_id)
    if samples.empty or "patientId" not in samples.columns:
        return pd.DataFrame(), "No samples"
    sample_to_patient = samples.set_index("sampleId")["patientId"].to_dict()

    # Get mutation groups (same logic as survival)
    muts_df = pd.DataFrame()
    try:
        muts_df = get_mutations(mut_profile_id, sample_list_id=sample_list_id, entrez_gene_ids=[entrez])
        if not muts_df.empty:
            muts_df = _add_gene_symbols(muts_df)
    except Exception:
        pass

    LOF_TYPES = {"frame_shift_del", "frame_shift_ins", "nonsense_mutation", "splice_site", "nonstop_mutation", "de_novo_start_outofframe", "start_codon_snp"}
    GOF_TYPES = {"missense_mutation", "in_frame_del", "in_frame_ins", "missense_variant"}

    def _classify(mt):
        mt = str(mt).lower().replace(" ", "_") if pd.notna(mt) else ""
        if mt in LOF_TYPES:
            return "Loss of Function"
        if mt in GOF_TYPES:
            return "Gain of Function"
        return "Gain of Function"

    sample_groups = {}
    all_sids = samples["sampleId"].tolist()
    if not muts_df.empty and "sampleId" in muts_df.columns:
        mut_col = "mutationType" if "mutationType" in muts_df.columns else "variantType"
        for _, row in muts_df.iterrows():
            sid = row["sampleId"]
            grp = _classify(row.get(mut_col, ""))
            if sid not in sample_groups or grp == "Loss of Function":
                sample_groups[sid] = grp
        for sid in list(sample_groups):
            if sample_groups[sid] == "Other":
                sample_groups[sid] = "Gain of Function"
    for sid in all_sids:
        if sid not in sample_groups:
            sample_groups[sid] = "Wild Type"

    downstream = DOWNSTREAM_GENES.get(gene_symbol.upper(), DEFAULT_DOWNSTREAM)
    downstream = [g for g in downstream if g.upper() != gene_symbol.upper()]
    down_entrez = [get_entrez_id(g) for g in downstream]
    down_entrez = [(e, s) for e, s in zip(down_entrez, downstream) if e]
    if not down_entrez:
        return pd.DataFrame(), "No downstream genes resolved"

    expr_df = get_molecular_data(expr_profile_id, sample_list_id, [e for e, _ in down_entrez])
    if expr_df.empty:
        return pd.DataFrame(), "No expression data"
    # API returns sampleId, value; gene id as entrezGeneId or geneId
    id_col = "entrezGeneId" if "entrezGeneId" in expr_df.columns else "geneId"
    if id_col not in expr_df.columns or "sampleId" not in expr_df.columns or "value" not in expr_df.columns:
        return pd.DataFrame(), "No expression data (missing sampleId/value/geneId)"
    expr_df["entrezGeneId"] = expr_df[id_col]

    expr_df["group"] = expr_df["sampleId"].map(sample_groups)
    expr_df = expr_df.dropna(subset=["group", "value"])
    expr_df["value"] = pd.to_numeric(expr_df["value"], errors="coerce")
    expr_df = expr_df.dropna(subset=["value"])

    entrez_to_symbol = {e: s for e, s in down_entrez}
    results = []
    try:
        from scipy.stats import mannwhitneyu
    except ImportError:
        return pd.DataFrame(), "Install scipy for DEG"

    for eg, sym in down_entrez:
        sub = expr_df[expr_df["entrezGeneId"] == eg]
        if sub.empty or sub["group"].nunique() < 2:
            continue
        for grp in ["Loss of Function", "Gain of Function"]:
            a = sub[sub["group"] == grp]["value"]
            b = sub[sub["group"] == "Wild Type"]["value"]
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                u, p = mannwhitneyu(a, b, alternative="two-sided")
                fc = a.median() - b.median()
                direction = "up" if fc > 0 else "down"
                results.append({
                    "Gene": sym,
                    "Comparison": f"{grp} vs Wild",
                    "Median_diff": round(fc, 3),
                    "p_value": round(p, 4),
                    "Direction": direction,
                    "N_mut": len(a),
                    "N_wild": len(b),
                })
            except Exception:
                continue

    if not results:
        return pd.DataFrame(), "No significant comparisons"
    deg_df = pd.DataFrame(results)
    deg_df = deg_df.sort_values("p_value")
    return deg_df, ""


def fetch_deg_full(
    study_id: str,
    mut_profile_id: str,
    gene_symbol: str,
    sample_list_id: Optional[str] = None,
    p_threshold: float = 0.05,
) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    """
    Full DEG analysis: classify samples into Wild, LoF, GoF; test all comparisons.
    Comparisons: LoF vs Wild, GoF vs Wild, LoF vs GoF.
    Returns (deg_df, group_counts_df, error_msg).
    deg_df: Gene, Comparison, Median_diff, p_value, Direction, N_A, N_B.
    """
    entrez = get_entrez_id(gene_symbol)
    if not entrez:
        return pd.DataFrame(), pd.DataFrame(), f"Gene {gene_symbol} not found"

    if not sample_list_id:
        sl_df = get_sample_lists(study_id)
        sample_list_id = _pick_sample_list_for_expression(sl_df)
        if not sample_list_id:
            return pd.DataFrame(), pd.DataFrame(), "No suitable sample list (need expression data)"

    profiles = get_molecular_profiles(study_id)
    if profiles.empty:
        return pd.DataFrame(), pd.DataFrame(), "No molecular profiles"
    expr_profiles = profiles[
        profiles["molecularAlterationType"].astype(str).str.upper().str.contains("MRNA|EXPRESSION", na=False)
    ]
    if expr_profiles.empty and "datatype" in profiles.columns:
        expr_profiles = profiles[profiles["datatype"].astype(str).str.upper().str.contains("Z-SCORE|CONTINUOUS", na=False)]
    if expr_profiles.empty:
        return pd.DataFrame(), pd.DataFrame(), "No expression profile in study"

    expr_profile_id = expr_profiles.iloc[0]["molecularProfileId"]
    samples = get_samples(study_id)
    if samples.empty:
        return pd.DataFrame(), pd.DataFrame(), "No samples"

    # Mutation classification
    muts_df = pd.DataFrame()
    try:
        muts_df = get_mutations(mut_profile_id, sample_list_id=sample_list_id, entrez_gene_ids=[entrez])
        if not muts_df.empty:
            muts_df = _add_gene_symbols(muts_df)
    except Exception:
        pass

    LOF_TYPES = {"frame_shift_del", "frame_shift_ins", "nonsense_mutation", "splice_site",
                 "nonstop_mutation", "de_novo_start_outofframe", "start_codon_snp"}
    GOF_TYPES = {"missense_mutation", "in_frame_del", "in_frame_ins", "missense_variant"}

    def _classify(mt):
        mt = str(mt).lower().replace(" ", "_") if pd.notna(mt) else ""
        if mt in LOF_TYPES:
            return "Loss of Function"
        if mt in GOF_TYPES:
            return "Gain of Function"
        return "Gain of Function"

    sample_groups = {}
    all_sids = samples["sampleId"].tolist()
    if not muts_df.empty and "sampleId" in muts_df.columns:
        mut_col = "mutationType" if "mutationType" in muts_df.columns else "variantType"
        for _, row in muts_df.iterrows():
            sid = row["sampleId"]
            grp = _classify(row.get(mut_col, ""))
            if sid not in sample_groups or grp == "Loss of Function":
                sample_groups[sid] = grp
        for sid in list(sample_groups):
            if sample_groups[sid] == "Other":
                sample_groups[sid] = "Gain of Function"
    for sid in all_sids:
        if sid not in sample_groups:
            sample_groups[sid] = "Wild Type"

    # Group counts
    from collections import Counter
    cnt = Counter(sample_groups.values())
    group_counts = pd.DataFrame([{"Group": k, "N": v} for k, v in sorted(cnt.items())])

    # Gene set: DEG_GENES_FULL minus query gene
    genes_to_test = [g for g in DEG_GENES_FULL if g.upper() != gene_symbol.upper()]
    down_entrez = [(get_entrez_id(g), g) for g in genes_to_test]
    down_entrez = [(e, s) for e, s in down_entrez if e]
    if not down_entrez:
        return pd.DataFrame(), group_counts, "No genes resolved for DEG"

    expr_df = get_molecular_data(expr_profile_id, sample_list_id, [e for e, _ in down_entrez])
    if expr_df.empty:
        return pd.DataFrame(), group_counts, "No expression data"
    id_col = "entrezGeneId" if "entrezGeneId" in expr_df.columns else "geneId"
    if id_col not in expr_df.columns or "sampleId" not in expr_df.columns or "value" not in expr_df.columns:
        return pd.DataFrame(), group_counts, "No expression data (missing columns)"
    expr_df["entrezGeneId"] = expr_df[id_col]
    expr_df["group"] = expr_df["sampleId"].map(sample_groups)
    expr_df = expr_df.dropna(subset=["group", "value"])
    expr_df["value"] = pd.to_numeric(expr_df["value"], errors="coerce")
    expr_df = expr_df.dropna(subset=["value"])

    try:
        from scipy.stats import mannwhitneyu
    except ImportError:
        return pd.DataFrame(), group_counts, "Install scipy for DEG"

    comparisons = [
        ("Loss of Function", "Wild Type", "LoF vs Wild"),
        ("Gain of Function", "Wild Type", "GoF vs Wild"),
        ("Loss of Function", "Gain of Function", "LoF vs GoF"),
    ]
    results = []
    for eg, sym in down_entrez:
        sub = expr_df[expr_df["entrezGeneId"] == eg]
        if sub.empty:
            continue
        for grp_a, grp_b, label in comparisons:
            a = sub[sub["group"] == grp_a]["value"]
            b = sub[sub["group"] == grp_b]["value"]
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                u, p = mannwhitneyu(a, b, alternative="two-sided")
                fc = a.median() - b.median()
                direction = "up" if fc > 0 else "down"
                results.append({
                    "Gene": sym,
                    "Comparison": label,
                    "Median_diff": round(fc, 3),
                    "p_value": round(p, 4),
                    "Direction": direction,
                    "N_A": len(a),
                    "N_B": len(b),
                })
            except Exception:
                continue

    if not results:
        return pd.DataFrame(), group_counts, "No comparisons (need ≥3 samples per group)"
    deg_df = pd.DataFrame(results)
    deg_df = deg_df[deg_df["p_value"] <= p_threshold]
    deg_df = deg_df.sort_values("p_value")
    return deg_df, group_counts, ""


def get_cancer_types() -> pd.DataFrame:
    """Fetch cancer types for filtering."""
    data = _get("/cancer-types")
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)
