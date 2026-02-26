"""
cBioPortal Explorer - Interactive Streamlit app for cancer genomics data.
Explore studies, datasets, run mutation analyses, and export high-resolution figures + tables.
"""

import io
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from cbioportal_client import (
    get_all_studies,
    get_molecular_profiles,
    get_samples,
    get_mutations,
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

st.set_page_config(
    page_title="cBioPortal Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1E3A5F;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        color: #5A7A9A;
        font-size: 1rem;
        margin-bottom: 2rem;
    }
    .stDownloadButton button {
        background: linear-gradient(90deg, #2563eb 0%, #1d4ed8 100%);
        color: white;
        font-weight: 600;
    }
    div[data-testid="stSidebar"] {
        background: linear-gradient(180deg, #f8fafc 0%, #f1f5f9 100%);
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<p class="main-header">🧬 cBioPortal Explorer</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="sub-header">Interactive cancer genomics data exploration with mutation analysis and publication-ready exports</p>',
    unsafe_allow_html=True,
)

# Sidebar: Data selection
st.sidebar.header("📋 Data Selection")

# Load studies
@st.cache_data(ttl=3600)
def load_studies():
    try:
        return get_all_studies(projection="SUMMARY")
    except Exception as e:
        st.error(f"Failed to load studies: {e}")
        return pd.DataFrame()

studies_df = load_studies()

if studies_df.empty:
    st.warning("Could not load studies from cBioPortal. Check your internet connection.")
    st.stop()

# Study dropdown - filter by cancer type if desired
cancer_filter = st.sidebar.text_input("Filter by keyword (e.g. TCGA, breast)", "")
if cancer_filter:
    mask = studies_df.apply(
        lambda r: cancer_filter.lower() in str(r.get("name", "")).lower()
        or cancer_filter.lower() in str(r.get("studyId", "")).lower()
        or cancer_filter.lower() in str(r.get("cancerTypeId", "")).lower(),
        axis=1,
    )
    studies_df = studies_df[mask]

if studies_df.empty:
    st.sidebar.info("No studies match your filter.")
    st.stop()

study_options = {
    f"{row.get('name', row.get('studyId', ''))} ({row.get('studyId', '')})": row.get("studyId")
    for _, row in studies_df.iterrows()
}
selected_study_label = st.sidebar.selectbox(
    "1. Select Cancer Study",
    options=list(study_options.keys()),
    key="study",
)
study_id = study_options[selected_study_label]

# Molecular profiles (datasets)
@st.cache_data(ttl=1800)
def load_profiles(sid):
    try:
        return get_molecular_profiles(sid)
    except Exception:
        return pd.DataFrame()

profiles_df = load_profiles(study_id)

# Filter to mutation profiles for this app
if not profiles_df.empty and "molecularAlterationType" in profiles_df.columns:
    mutation_profiles = profiles_df[
        profiles_df["molecularAlterationType"].str.upper().str.contains("MUTATION", na=False)
    ]
    if mutation_profiles.empty:
        mutation_profiles = profiles_df
else:
    mutation_profiles = profiles_df

if mutation_profiles.empty:
    st.warning(f"No molecular profiles found for study **{study_id}**. Try another study.")
    st.stop()

profile_options = {
    f"{row.get('name', row.get('molecularProfileId', ''))}": row.get("molecularProfileId")
    for _, row in mutation_profiles.iterrows()
}
selected_profile_label = st.sidebar.selectbox(
    "2. Select Mutation Dataset",
    options=list(profile_options.keys()),
    key="profile",
)
molecular_profile_id = profile_options[selected_profile_label]

# Analysis type
analysis_options = list(ANALYSIS_TYPES.keys())
selected_analysis = st.sidebar.selectbox(
    "3. Select Analysis",
    options=analysis_options,
    key="analysis",
)

# Optional: limit samples for faster loading
sample_limit = st.sidebar.slider("Max samples to analyze (for speed)", 50, 500, 200)
top_n = st.sidebar.slider("Top N items in charts", 5, 50, 20)

st.sidebar.divider()
st.sidebar.caption("Data from [cBioPortal](https://www.cbioportal.org)")

# Main content
if st.sidebar.button("▶️ Run Analysis", type="primary"):
    with st.spinner("Fetching mutation data..."):
        try:
            mutations_df = get_mutations(
                molecular_profile_id,
                study_id=study_id,
                max_samples=sample_limit,
            )
        except Exception as e:
            st.error(f"Error fetching mutations: {e}")
            st.stop()

    if mutations_df.empty:
        st.warning("No mutation data returned. The dataset may be empty or require gene filters.")
        st.stop()

    # Limit samples for large datasets
    if "sampleId" in mutations_df.columns:
        sample_ids = mutations_df["sampleId"].unique().tolist()
        if len(sample_ids) > sample_limit:
            sample_ids = sample_ids[:sample_limit]
            mutations_df = mutations_df[mutations_df["sampleId"].isin(sample_ids)]

    st.success(f"Loaded **{len(mutations_df)}** mutations from **{mutations_df['sampleId'].nunique() if 'sampleId' in mutations_df.columns else '?'}** samples.")

    # Summary statistics
    st.subheader("📊 Summary Statistics")
    summary_df = summary_statistics(mutations_df)
    st.dataframe(summary_df, use_container_width=True, hide_index=True)

    csv_summary = summary_df.to_csv(index=False)
    st.download_button(
        "Download Summary (CSV)",
        csv_summary,
        file_name=f"{study_id}_summary_stats.csv",
        mime="text/csv",
        key="dl_summary",
    )

    st.divider()
    st.subheader(f"📈 {selected_analysis}")

    # Run selected analysis
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

        # High-resolution download
        buf = io.BytesIO()
        fig2, _ = analysis_func(mutations_df, top_n=top_n)
        fig2.savefig(buf, format="png", dpi=300, bbox_inches="tight")
        plt.close(fig2)
        buf.seek(0)
        st.download_button(
            "⬇️ Download Figure (PNG, 300 DPI)",
            buf,
            file_name=f"{study_id}_{selected_analysis.replace(' ', '_').replace('(', '').replace(')', '')}.png",
            mime="image/png",
            key="dl_fig",
        )

    with col2:
        st.dataframe(stats_df, use_container_width=True, hide_index=True)
        if not stats_df.empty:
            csv_stats = stats_df.to_csv(index=False)
            st.download_button(
                "Download Table (CSV)",
                csv_stats,
                file_name=f"{study_id}_{selected_analysis.replace(' ', '_')}_stats.csv",
                mime="text/csv",
                key="dl_table",
            )

    # Raw data download
    st.divider()
    with st.expander("Download raw mutation data"):
        raw_csv = mutations_df.to_csv(index=False)
        st.download_button(
            "Download Full Mutation Data (CSV)",
            raw_csv,
            file_name=f"{study_id}_mutations.csv",
            mime="text/csv",
            key="dl_raw",
        )
else:
    st.info("👆 Select your study, dataset, and analysis type in the sidebar, then click **Run Analysis**.")
    st.markdown("""
    ### Quick Start
    - **Study**: Pick a cancer cohort (e.g. TCGA studies)
    - **Dataset**: Choose the mutation profile for that study
    - **Analysis**: Select the type of mutation figure/table to generate

    All figures can be downloaded at **300 DPI** for publication use.
    """)
