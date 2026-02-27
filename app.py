"""
cBioPortal Explorer - Interactive Streamlit app for cancer genomics data.
Explore studies, datasets, run mutations, and export figures + tables.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import io
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from cbioportal_client import (
    get_all_studies,
    get_molecular_profiles,
    get_samples,
    get_mutations,
    fetch_mutations_by_study,
    fetch_gene_across_studies,
    fetch_survival_data_for_gene,
    fetch_subtype_enrichment,
    fetch_deg_downstream,
    filter_studies_with_survival,
)
from visualizations import (
    ANALYSIS_TYPES,
    mutation_type_distribution,
    gene_mutation_frequency,
    variant_classification_distribution,
    sample_mutation_burden,
    oncoprint_style_matrix,
    summary_statistics,
    gene_across_studies_bar,
    lollipop_mutations,
    subtype_enrichment_chart,
    survival_plot_gof_lof,
    deg_downstream_chart,
)

# Owner info
AUTHOR_NAME = "Mahmood Al Mashhadani"
AUTHOR_TITLE = "Assistant Professor - Molecular Medicine"
AUTHOR_AFFILIATION = "MBRU - College of Medicine"
AUTHOR_EMAIL = "mahmood.almashhadani@dubaihealth.ae"

st.set_page_config(
    page_title="cBioPortal Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    .banner-top { background: linear-gradient(135deg, #0f4c81 0%, #1E3A5F 50%, #2d5a87 100%); color: white; padding: 1.25rem 1.5rem; border-radius: 8px; margin-bottom: 1rem; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
    .banner-top h1 { margin: 0; font-size: 1.8rem; font-weight: 700; }
    .banner-top p { margin: 0.25rem 0 0 0; font-size: 0.95rem; opacity: 0.95; }
    .welcome-box { background: linear-gradient(135deg, #f0f7ff 0%, #e8f4fd 100%); padding: 1rem 1.5rem; border-radius: 8px; border-left: 4px solid #1E3A5F; margin-bottom: 1rem; }
    .data-source-banner { background: #fff8e6; border: 1px solid #f0d675; border-radius: 8px; padding: 1rem 1.25rem; margin-bottom: 1.25rem; font-size: 0.9rem; }
    .data-source-banner a { color: #1E3A5F; font-weight: 600; }
    .footer-legal { font-size: 0.75rem; color: #6b7280; margin-top: 2rem; padding: 1rem; border-top: 1px solid #e5e7eb; }
</style>
""", unsafe_allow_html=True)

# Main banner
st.markdown("""
<div class="banner-top">
    <h1>🧬 cBioPortal Explorer – Breast Cancer</h1>
    <p>Mutation analysis for breast cancer genomics with hereditary gene panels</p>
</div>
""", unsafe_allow_html=True)

# Welcome / Author banner
st.markdown(f"""
<div class="welcome-box">
    <strong>Welcome to the tool.</strong><br>
    Developed by <strong>{AUTHOR_NAME}</strong> — {AUTHOR_TITLE}<br>
    <small>{AUTHOR_AFFILIATION} • <a href="mailto:{AUTHOR_EMAIL}">{AUTHOR_EMAIL}</a></small>
</div>
""", unsafe_allow_html=True)

# cBioPortal data source attribution (per cBioPortal requirements for legal protection)
st.markdown("""
<div class="data-source-banner">
    <strong>📊 Data Source:</strong> All genomic data is retrieved from 
    <a href="https://www.cbioportal.org" target="_blank">cBioPortal for Cancer Genomics</a> via their public API.
    This tool does not host or store data; it fetches data on demand from cBioPortal.<br><br>
    <strong>Citation (required when using cBioPortal data):</strong><br>
    • Cerami E, et al. The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data. <em>Cancer Discovery</em>. 2012;2(5):401-404.<br>
    • Gao J, et al. Integrative Analysis of Complex Cancer Genomics and Clinical Profiles Using the cBioPortal. <em>Science Signaling</em>. 2013;6(269):pl1.<br>
    <small>Full citation info: <a href="https://docs.cbioportal.org/user-guide/faq/#how-do-i-cite-the-cbioportal" target="_blank">cBioPortal Citation Guide</a></small>
</div>
""", unsafe_allow_html=True)

st.sidebar.header("📋 Mode")
mode = st.sidebar.radio("Choose analysis mode", ["Study Analysis", "Gene Search Across Studies", "Survival Plotter (GoF vs LoF vs Wild)"], key="mode")
st.sidebar.caption("Data from [cBioPortal](https://www.cbioportal.org)")

@st.cache_data(ttl=3600)
def load_breast_cancer_studies():
    try:
        df = get_all_studies(projection="SUMMARY")
        if df.empty:
            return df
        keywords = ("breast", "brca", "mammary")
        mask = df.apply(
            lambda r: any(kw in str(r.get("name","")).lower() or kw in str(r.get("studyId","")).lower() or kw in str(r.get("cancerTypeId","")).lower() for kw in keywords),
            axis=1,
        )
        return df[mask].sort_values("studyId").reset_index(drop=True)
    except Exception as e:
        st.error(f"Failed to load studies: {e}")
        return pd.DataFrame()


@st.cache_data(ttl=3600)
def load_breast_cancer_studies_with_survival():
    """Breast cancer studies that have patient-level survival data (OS/DFS/PFS/RFS)."""
    df = load_breast_cancer_studies()
    if df.empty:
        return df
    return filter_studies_with_survival(df)


studies_df = load_breast_cancer_studies()
if studies_df.empty:
    st.warning("Could not load breast cancer studies from cBioPortal.")
    st.stop()

if mode == "Gene Search Across Studies":
    st.subheader("🔍 Gene Search Across All Breast Cancer Studies")
    st.caption("Queries all sequenced samples per study (matches cBioPortal scope). Breast cancer studies only.")
    gene_input = st.text_input("Enter gene symbol (e.g., BRCA1, TP53, PIK3CA)", value="BRCA1", key="gene_search")
    max_studies_search = st.sidebar.slider("Max studies to query", 5, 80, 30, key="max_studies_search")
    if st.button("🔍 Search Gene", type="primary", key="btn_gene_search"):
        with st.spinner(f"Searching {gene_input} across {max_studies_search} studies..."):
            try:
                muts_df, counts_df = fetch_gene_across_studies(gene_input, studies_df, max_studies=max_studies_search)
            except Exception as e:
                st.error(f"Error: {e}")
                muts_df, counts_df = pd.DataFrame(), pd.DataFrame()
        if muts_df.empty and counts_df.empty:
            st.warning(f"No mutations found for **{gene_input}** in breast cancer studies. Check gene symbol.")
        else:
            st.success(f"Found **{len(muts_df)}** mutations for **{gene_input}** across **{len(counts_df)}** studies.")
            col_bar, col_lollipop = st.columns(2)
            with col_bar:
                fig_bar, stats_bar = gene_across_studies_bar(counts_df, gene_input)
                st.pyplot(fig_bar)
                plt.close(fig_bar)
                st.dataframe(stats_bar, use_container_width=True, hide_index=True)
            with col_lollipop:
                fig_lolli, stats_lolli = lollipop_mutations(muts_df, gene_input)
                st.pyplot(fig_lolli)
                plt.close(fig_lolli)
                if not stats_lolli.empty:
                    st.dataframe(stats_lolli, use_container_width=True, hide_index=True)
            st.subheader("Molecular subtype association")
            with st.spinner("Fetching subtype data for top studies..."):
                subtype_results = fetch_subtype_enrichment(muts_df, studies_df, top_n_studies=5)
            if subtype_results:
                fig_sub, stats_sub = subtype_enrichment_chart(subtype_results, gene_input)
                st.pyplot(fig_sub)
                plt.close(fig_sub)
                st.markdown("**Chi-squared test:** Mutation rate (%) by subtype. p < 0.05 = significant association.")
                for res in subtype_results:
                    pv = res["p_value"]
                p_str = f"{pv:.4f}" if not pd.isna(pv) else "N/A"
                sig = " (significant)" if not pd.isna(pv) and pv < 0.05 else ""
                with st.expander(f"{res['studyName'][:65]} - p = {p_str}{sig}", expanded=not pd.isna(pv) and pv < 0.05):
                        st.dataframe(res["subtype_df"], use_container_width=True, hide_index=True)
            else:
                st.info("No molecular subtype data available for studies with mutations.")
            st.subheader("All variations")
            disp_cols = [c for c in muts_df.columns if c in ["geneSymbol", "hugoGeneSymbol", "proteinChange", "aminoAcidChange", "mutationType", "variantType", "studyId", "studyName"]]
            if disp_cols:
                st.dataframe(muts_df[disp_cols].drop_duplicates(), use_container_width=True, hide_index=True)
            st.download_button("Download all mutations (CSV)", muts_df.to_csv(index=False), file_name=f"{gene_input}_across_studies.csv", mime="text/csv", key="dl_gene_muts")
elif mode == "Survival Plotter (GoF vs LoF vs Wild)":
    st.subheader("📈 Survival by Mutation Type (GoF vs LoF vs Wild Type)")
    st.markdown("Compare **Gain of Function** (missense, in-frame), **Loss of Function** (nonsense, frameshift, splice), and **Wild Type** (no mutation) for a gene. Only studies with survival data (OS/DFS/PFS/RFS) are shown.")
    with st.spinner("Loading studies with survival data (may take 30-60 sec on first load)..."):
        surv_studies_df = load_breast_cancer_studies_with_survival()
    if surv_studies_df.empty:
        st.warning("No breast cancer studies with survival data found in cBioPortal.")
        st.stop()
    surv_study_options = {f"{r.get('name','')} ({r.get('studyId','')})": r.get("studyId") for _, r in surv_studies_df.iterrows()}
    surv_study_label = st.selectbox("Select study (with survival data)", options=list(surv_study_options.keys()), key="surv_study")
    surv_study_id = surv_study_options[surv_study_label]
    surv_profiles_df = get_molecular_profiles(surv_study_id)
    if not surv_profiles_df.empty and "molecularAlterationType" in surv_profiles_df.columns:
        surv_mut_profiles = surv_profiles_df[surv_profiles_df["molecularAlterationType"].str.upper().str.contains("MUTATION", na=False)]
        surv_mut_profiles = surv_mut_profiles if not surv_mut_profiles.empty else surv_profiles_df
    else:
        surv_mut_profiles = surv_profiles_df
    surv_profile_options = {r.get("name", r.get("molecularProfileId","")): r.get("molecularProfileId") for _, r in surv_mut_profiles.iterrows()} if not surv_mut_profiles.empty else {}
    surv_profile_label = st.selectbox("Select mutation profile", options=list(surv_profile_options.keys()) or [""], key="surv_profile")
    surv_profile_id = surv_profile_options.get(surv_profile_label) if surv_profile_options else None
    surv_gene_input = st.text_input("Enter gene symbol", value="BRCA1", key="surv_gene")
    if st.button("📈 Plot Survival", type="primary", key="btn_surv"):
        if not surv_profile_id:
            st.error("No mutation profile selected.")
        else:
            with st.spinner("Fetching survival data..."):
                surv_df, surv_counts, time_col, err = fetch_survival_data_for_gene(surv_study_id, surv_profile_id, surv_gene_input)
            if err:
                st.warning(err)
            else:
                st.success(f"Loaded survival data for **{surv_gene_input}** ({time_col})")
                st.dataframe(surv_counts, use_container_width=True, hide_index=True)
                surv_type_display = {"OS_MONTHS": "Overall Survival", "DFS_MONTHS": "Disease-Free Survival",
                    "PFS_MONTHS": "Progression-Free Survival", "RFS_MONTHS": "Recurrence-Free Survival",
                    "DSS_MONTHS": "Disease-Specific Survival"}.get(time_col, time_col)

                # Overall survival (all patients)
                st.subheader("Overall (all patients)")
                fig_surv, stats_surv, interp = survival_plot_gof_lof(surv_df, surv_gene_input, surv_type_display)
                st.pyplot(fig_surv)
                plt.close(fig_surv)
                if not stats_surv.empty:
                    st.dataframe(stats_surv, use_container_width=True, hide_index=True)
                if interp:
                    st.markdown("---")
                    st.markdown(interp)

                # Per-subtype survival (if subtype data exists)
                st.subheader("Per molecular subtype")
                has_subtypes = (
                    "subtype" in surv_df.columns
                    and surv_df["subtype"].nunique() > 1
                    and (surv_df["subtype"].astype(str) != "All").any()
                )
                if has_subtypes:
                    subtypes = [
                        s for s in surv_df["subtype"].dropna().unique()
                        if str(s) not in ("All", "Unknown", "nan")
                    ]
                    min_per_subtype = 15
                    shown = 0
                    for sub in sorted(subtypes):
                        sub_df = surv_df[surv_df["subtype"] == sub]
                        if len(sub_df) < min_per_subtype or sub_df["group"].nunique() < 2:
                            continue
                        with st.expander(f"**{sub}** (n={len(sub_df)})", expanded=True):
                            fig_sub, stats_sub, interp_sub = survival_plot_gof_lof(
                                sub_df, surv_gene_input, f"{surv_type_display} - {sub}"
                            )
                            st.pyplot(fig_sub)
                            plt.close(fig_sub)
                            if not stats_sub.empty:
                                st.dataframe(stats_sub, use_container_width=True, hide_index=True)
                            if interp_sub:
                                st.markdown(interp_sub)
                        shown += 1
                    if shown == 0:
                        st.info("No subtypes with enough patients and mutation groups for separate analysis.")
                else:
                    st.info("Molecular subtype data not available for this study. Only overall survival is shown.")

                # DEG downstream of gene (Wild vs LoF vs GoF)
                st.subheader(f"DEG downstream of {surv_gene_input}")
                with st.spinner("Fetching expression data for downstream genes..."):
                    deg_df, deg_err = fetch_deg_downstream(surv_study_id, surv_profile_id, surv_gene_input)
                if deg_err:
                    st.info(deg_err)
                elif not deg_df.empty:
                    fig_deg, _ = deg_downstream_chart(deg_df, surv_gene_input)
                    st.pyplot(fig_deg)
                    plt.close(fig_deg)
                    st.dataframe(deg_df, use_container_width=True, hide_index=True)

                st.markdown("---")
                buf = io.BytesIO()
                fig_surv2, _, _ = survival_plot_gof_lof(surv_df, surv_gene_input, surv_type_display)
                fig_surv2.savefig(buf, format="png", dpi=300, bbox_inches="tight")
                plt.close(fig_surv2)
                buf.seek(0)
                st.download_button("Download survival plot (PNG)", buf, file_name=f"{surv_study_id}_{surv_gene_input}_survival.png", mime="image/png", key="dl_surv")
else:
    study_options = {
        f"{row.get('name', row.get('studyId', ''))} ({row.get('studyId', '')})": row.get("studyId")
        for _, row in studies_df.iterrows()
    }
    selected_study_label = st.sidebar.selectbox("1. Select Cancer Study", options=list(study_options.keys()), key="study")
    study_id = study_options[selected_study_label]

    @st.cache_data(ttl=1800)
    def load_profiles(sid):
        try:
            return get_molecular_profiles(sid)
        except Exception:
            return pd.DataFrame()

    profiles_df = load_profiles(study_id)
    if not profiles_df.empty and "molecularAlterationType" in profiles_df.columns:
        mutation_profiles = profiles_df[profiles_df["molecularAlterationType"].str.upper().str.contains("MUTATION", na=False)]
        if mutation_profiles.empty:
            mutation_profiles = profiles_df
    else:
        mutation_profiles = profiles_df

    if mutation_profiles.empty:
        st.warning(f"No molecular profiles found for study **{study_id}**.")
        st.stop()

    profile_options = {
        f"{row.get('name', row.get('molecularProfileId', ''))}": row.get("molecularProfileId")
        for _, row in mutation_profiles.iterrows()
    }
    selected_profile_label = st.sidebar.selectbox("2. Select Mutation Dataset", options=list(profile_options.keys()), key="profile")
    molecular_profile_id = profile_options[selected_profile_label]

    analysis_options = list(ANALYSIS_TYPES.keys())
    selected_analysis = st.sidebar.selectbox("3. Select Analysis", options=analysis_options, key="analysis")
    sample_limit = st.sidebar.slider("Max samples to analyze (for speed)", 50, 500, 200)
    top_n = st.sidebar.slider("Top N items in charts", 5, 50, 20)

    st.sidebar.divider()
    st.sidebar.markdown("**Data from [cBioPortal](https://www.cbioportal.org)** (see attribution below)")

    with st.sidebar.expander("© Legal & Data Usage", expanded=False):
        st.markdown(f"""
**Copyright © {AUTHOR_NAME}. All rights reserved.**

**This Tool:**
- This software is the exclusive property of the author.
- Unauthorized copying, modification, or distribution is prohibited.
- Contact: **{AUTHOR_EMAIL}**

**Data (cBioPortal):**
- Genomic data is sourced from cBioPortal. Users must cite cBioPortal per their [citation guidelines](https://docs.cbioportal.org/user-guide/faq/#how-do-i-cite-the-cbioportal).
- This tool does not guarantee data accuracy; refer to original studies and cBioPortal.

**Medical Disclaimer:**
This tool is for **research and educational purposes only**. It is NOT intended for clinical decision-making. Always consult qualified healthcare professionals for medical advice.
        """)

    if st.sidebar.button("▶️ Run Analysis", type="primary"):
        with st.spinner("Fetching mutation data..."):
            try:
                mutations_df = fetch_mutations_by_study(study_id, molecular_profile_id)
            except Exception as e:
                st.error(f"Error fetching mutations: {e}")
                st.stop()

        if mutations_df.empty:
            st.warning("No mutation data returned. Try another study.")
            st.stop()

        if "sampleId" in mutations_df.columns and mutations_df["sampleId"].nunique() > sample_limit:
            top_samples = mutations_df["sampleId"].value_counts().head(sample_limit).index.tolist()
            mutations_df = mutations_df[mutations_df["sampleId"].isin(top_samples)]

        st.success(f"Loaded **{len(mutations_df)}** mutations from **{mutations_df['sampleId'].nunique() if 'sampleId' in mutations_df.columns else '?'}** samples.")

        st.subheader("📊 Summary Statistics")
        summary_df = summary_statistics(mutations_df)
        st.dataframe(summary_df, use_container_width=True, hide_index=True)
        st.download_button("Download Summary (CSV)", summary_df.to_csv(index=False), file_name=f"{study_id}_summary_stats.csv", mime="text/csv", key="dl_summary")

        st.divider()
        st.subheader(f"📈 {selected_analysis}")

        analysis_func = ANALYSIS_TYPES[selected_analysis]
        try:
            fig, stats_df = analysis_func(mutations_df, top_n=top_n)
        except Exception as e:
            st.error(f"Analysis error: {e}")
            st.stop()

        col1, col2 = st.columns([2, 1])
        with col1:
            st.pyplot(fig)
            plt.close(fig)
            buf = io.BytesIO()
            fig2, _ = analysis_func(mutations_df, top_n=top_n)
            fig2.savefig(buf, format="png", dpi=300, bbox_inches="tight")
            plt.close(fig2)
            buf.seek(0)
            st.download_button("⬇️ Download Figure (PNG, 300 DPI)", buf, file_name=f"{study_id}_{selected_analysis.replace(' ', '_').replace('(', '').replace(')', '')}.png", mime="image/png", key="dl_fig")
        with col2:
            st.dataframe(stats_df, use_container_width=True, hide_index=True)
            if not stats_df.empty:
                st.download_button("Download Table (CSV)", stats_df.to_csv(index=False), file_name=f"{study_id}_{selected_analysis.replace(' ', '_')}_stats.csv", mime="text/csv", key="dl_table")

        st.divider()
        with st.expander("Download raw mutation data"):
            st.download_button("Download Full Mutation Data (CSV)", mutations_df.to_csv(index=False), file_name=f"{study_id}_mutations.csv", mime="text/csv", key="dl_raw")
    else:
        st.info("👆 Select study, dataset, and analysis, then click **Run Analysis**.")
    st.markdown("### Quick Start\n- **Study**: Pick a cancer cohort\n- **Dataset**: Mutation profile\n- **Analysis**: Figure type\n\nFigures at **300 DPI** for publication.")

st.divider()
st.markdown(f"""
<div class="footer-legal">
    <strong>© {AUTHOR_NAME}</strong> | All rights reserved. Data from <a href="https://www.cbioportal.org">cBioPortal</a> — cite per their guidelines. 
    Research use only; not for clinical decisions. Contact: {AUTHOR_EMAIL}
</div>
""", unsafe_allow_html=True)
