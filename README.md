# cBioPortal Explorer

Interactive Streamlit app for exploring cancer genomics data from [cBioPortal](https://www.cbioportal.org). Pick studies, mutation datasets, run analyses, and export high-resolution figures and statistics tables.

## Features

- **Study selection**: Browse and filter cancer studies from cBioPortal
- **Dataset selection**: Choose mutation profiles for your selected study
- **Analysis types**:
  - Mutation Type Distribution
  - Top Mutated Genes
  - Variant Classification (pie chart)
  - Sample Mutation Burden
  - Mutation Matrix (Gene × Sample heatmap)
- **Downloads**:
  - 300 DPI figures (PNG) for publication
  - Statistics tables (CSV)
  - Raw mutation data (CSV)

## Quick Start

```bash
# Create virtual environment (optional)
python -m venv venv
venv\Scripts\activate   # Windows
# source venv/bin/activate  # Linux/Mac

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

Open http://localhost:8501 in your browser.

## Project Structure

```
cbioportal-explorer/
├── app.py              # Main Streamlit app
├── cbioportal_client.py # cBioPortal API client
├── visualizations.py   # Figure and table generators
├── requirements.txt
└── README.md
```

## Live App (Share with Anyone)

**GitHub:** https://github.com/drmahmoodhachim-gif/cbioportal-explorer

**Deploy on Streamlit Cloud** (free, ~2 min):
1. Go to https://share.streamlit.io
2. Sign in with GitHub
3. New app → Repo: `drmahmoodhachim-gif/cbioportal-explorer`, Branch: `master`, File: `app.py`
4. Click **Deploy** → Get live URL like `https://cbioportal-explorer.streamlit.app`

## Clone and Run

```bash
git clone https://github.com/drmahmoodhachim-gif/cbioportal-explorer.git
cd cbioportal-explorer
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## Data Source

All data is fetched from the [cBioPortal Public API](https://www.cbioportal.org/api/swagger-ui/index.html). No API key required.
