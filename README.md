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

## Deploy as Live Page (Streamlit Cloud)

Free hosting on [Streamlit Community Cloud](https://share.streamlit.io):

1. **Push to GitHub** – Upload this repo to your GitHub account.
2. **Go to** [share.streamlit.io](https://share.streamlit.io)
3. **Sign in** with your GitHub account.
4. **New app** → choose your repo, branch (e.g. `main`), and set main file to `app.py`.
5. **Deploy** – Streamlit Cloud will build and host your app.

You’ll get a live URL like `https://your-app-name.streamlit.app`.

## GitHub

Clone and run:

```bash
git clone https://github.com/YOUR_USERNAME/cbioportal-explorer.git
cd cbioportal-explorer
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## Data Source

All data is fetched from the [cBioPortal Public API](https://www.cbioportal.org/api/swagger-ui/index.html). No API key required.
