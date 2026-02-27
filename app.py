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
)
from visualizations import (
    ANALYSIS_TYPES,
    mutation_type_distribution,
    gene_mutation_frequency,
    variant_classification_distribution,
    sample_mutation_burden,
    oncoprint_style_matrix,
    summary_statistics,
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
    .main-header { font-size: 2.2rem; font-weight: 700; color: #1E3A5F; margin-bottom: 0.5rem; }
    .sub-header { color: #5A7A9A; font-size: 1rem; margin-bottom: 0.5rem; }
    .welcome-box { background: linear-gradient(135deg, #f0f7ff 0%, #e8f4fd 100%); padding: 1rem 1.5rem; border-radius: 8px; border-left: 4px solid #1E3A5F; margin-bottom: 1.5rem; }
    .footer-legal { font-size: 0.75rem; color: #6b7280; margin-top: 2rem; padding: 1rem; border-top: 1px solid #e5e7eb; }
</style>
""", unsafe_allow_html=True)

st.markdown('<p class="main-header">🧬 cBioPortal Explorer – Breast Cancer</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Mutation analysis for breast cancer genomics with hereditary gene panels</p>', unsafe_allow_html=True)

st.markdown(f"""
<div class="welcome-box">
    <strong>Welcome to the tool.</strong><br>
    Developed by <strong>{AUTHOR_NAME}</strong> — {AUTHOR_TITLE}<br>
    <small>{AUTHOR_AFFILIATION} • <a href="mailto:{AUTHOR_EMAIL}">{AUTHOR_EMAIL}</a></small>
</div>
""", unsafe_allow_html=True)

st.sidebar.header("📋 Data Selection")

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

studies_df = load_breast_cancer_studies()
if studies_df.empty:
    st.warning("Could not load breast cancer studies from cBioPortal.")
    st.stop()

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
st.sidebar.caption("Data from [cBioPortal](https://www.cbioportal.org)")

with st.sidebar.expander("© Legal & Property", expanded=False):
    st.markdown(f"""
    **Copyright © {AUTHOR_NAME}. All rights reserved.**

    - This software and its content are the exclusive property of the author.
    - Unauthorized copying, modification, distribution, or use is strictly prohibited.
    - For permission to use, reproduce, or cite, contact: **{AUTHOR_EMAIL}**
    - **Disclaimer:** This tool is for research and educational purposes only. It is not intended for clinical decision-making. Always consult qualified healthcare professionals for medical advice.
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
    © {AUTHOR_NAME} | All rights reserved. This tool is proprietary. No unauthorized use, reproduction, or distribution permitted. Contact: {AUTHOR_EMAIL}
</div>
""", unsafe_allow_html=True)
