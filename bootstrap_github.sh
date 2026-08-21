#!/usr/bin/env bash
# One-shot: create private repo under martinfrasch and push this folder.
# Requires GitHub CLI authenticated: gh auth login
set -euo pipefail
REPO="${1:-nonrecip-modularity}"
git init -b main
cat > .gitignore << 'GI'
data/
venv/
__pycache__/
*.npz
*.csv
*.png
!pilot/*.png
GI
git add .
git commit -m "Nonreciprocal colloid modularity test (Hara et al. PRL 137, 068302) + NWAP analysis pipeline"
gh repo create "martinfrasch/${REPO}" --private --source=. --push
echo "Done: https://github.com/martinfrasch/${REPO}"
