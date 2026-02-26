# Opens Streamlit Cloud deploy page - sign in with GitHub and click Deploy
$url = "https://share.streamlit.io/deploy?repository=drmahmoodhachim-gif/cbioportal-explorer"
Start-Process $url
Write-Host "Browser opened. Sign in with GitHub, then click Deploy."
Write-Host "Your app will be live at: https://cbioportal-explorer.streamlit.app"
