# Deploy with Private Repo (Option 1)

This guide walks you through making the code private while keeping the app usable via Streamlit Teams or similar paid hosting.

---

## Step 1: Make the GitHub repository private

1. Go to your repo: **https://github.com/drmahmoodhachim-gif/cbioportal-explorer**
2. Click **Settings**
3. Scroll to **Danger Zone**
4. Click **Change repository visibility**
5. Select **Make private**
6. Confirm by typing the repository name

⚠️ **Note:** After this, only you and collaborators you add can see the code. The current public Streamlit deployment will **stop working** (Community Cloud requires public repos).

---

## Step 2: Set up Streamlit Teams (for private repo deployment)

Streamlit Teams supports deploying from private GitHub repositories.

1. Go to **https://share.streamlit.io**
2. Sign in with your GitHub account
3. Click **Upgrade** or **Get Streamlit Teams** (or **Streamlit for Teams**)
4. Choose a plan and complete billing
5. In the dashboard:
   - Click **New app**
   - Select your **private** repo: `drmahmoodhachim-gif/cbioportal-explorer`
   - Ensure your GitHub account has access to the repo
   - Choose branch: `master` (or `main`)
   - Main file: `app.py`
   - Click **Deploy**

Your app will deploy from the private repo. Users access the app URL; they never see the source code.

---

## Alternative: Streamlit Community Cloud + separate private repo

If you don't want to pay for Streamlit Teams:

1. **Create a new private repo** (e.g. `cbioportal-explorer-private`) with your code
2. **Keep the current public repo** as a minimal “shell” for Community Cloud:
   - Remove or replace the actual app code with a redirect or lightweight placeholder
   - Or: use a different hosting service that supports private repos

**Paid alternatives that support private repos:**

- **Streamlit Teams** – https://streamlit.io/cloud
- **Hugging Face Spaces** (private Spaces on paid plan) – https://huggingface.co/spaces
- **Render** – https://render.com (can connect private GitHub repos on paid plans)
- **Railway** – https://railway.app
- **DigitalOcean App Platform** – connect private GitHub repo

---

## Summary

| Action | Where |
|--------|-------|
| Make repo private | GitHub → Settings → Danger Zone → Change visibility |
| Deploy from private repo | Streamlit Teams / Hugging Face / Render / Railway / etc. |
| Share app | Share the deployment URL (e.g. `https://yourapp.streamlit.app`) |

---

## Quick checklist

- [ ] Make GitHub repo private
- [ ] Sign up for Streamlit Teams (or alternative)
- [ ] Connect private repo to hosting
- [ ] Deploy app
- [ ] Test app URL
- [ ] Share URL with users (they use the app; they cannot see your code)
