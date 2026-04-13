#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# One-time setup for Protein Sequence Alignment Tool (macOS)
# Run once from Terminal:   bash setup_mac.sh
# ─────────────────────────────────────────────────────────────────────────────
set -e

echo "=== Protein Sequence Alignment Tool — macOS Setup ==="

# 1. Check Python 3
if ! command -v python3 &>/dev/null; then
  echo "ERROR: python3 not found. Install from https://www.python.org"
  exit 1
fi
echo "✓ Python $(python3 --version)"

# 2. Install Homebrew if needed (for HMMER)
if ! command -v brew &>/dev/null; then
  echo "Installing Homebrew..."
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi

# 3. Install HMMER (required by abnumber for CDR annotation)
if ! command -v hmmscan &>/dev/null; then
  echo "Installing HMMER (required for CDR annotation)..."
  brew install hmmer
else
  echo "✓ HMMER already installed"
fi

# 4. Install Python packages
echo "Installing Python packages..."
pip3 install abnumber biopython --upgrade

echo ""
echo "=== Setup complete! ==="
echo "To launch the tool, run:"
echo "   python3 protein_aligner.py"
echo ""
echo "Or double-click 'run.command' in Finder."
