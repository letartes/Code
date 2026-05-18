#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  Protein Liability Analyzer — macOS Build Script
#  Creates a .app bundle and a distributable .dmg disk image.
#
#  Run from a terminal:
#      cd path/to/Code/build
#      chmod +x build_mac.sh
#      ./build_mac.sh
#
#  Requires: Python 3.8+, pip
#  Optional: brew install create-dmg   (prettier DMG with drag-to-install UI)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

APP_NAME="ProteinLiabilityAnalyzer"
APP_DISPLAY="Protein Liability Analyzer"
VERSION="2.0"
BUNDLE_ID="com.letartescientific.proteinliabilityanalyzer"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(dirname "$SCRIPT_DIR")"          # Code/
DIST_DIR="$SRC_DIR/dist"

echo ""
echo "======================================================"
echo "  $APP_DISPLAY — macOS Build"
echo "  Letarte Scientific Consulting"
echo "======================================================"
echo ""

# ── Check Python ──────────────────────────────────────────────────────────────
if ! command -v python3 &>/dev/null; then
    echo "[ERROR] python3 not found."
    exit 1
fi
echo "Python: $(python3 --version)"

# ── Check PyInstaller and tk support ─────────────────────────────────────────
echo ""
echo "[1/4]  Checking PyInstaller..."
if command -v pyinstaller &>/dev/null; then
    echo "  PyInstaller $(pyinstaller --version) found at $(which pyinstaller)"
elif python3 -m PyInstaller --version &>/dev/null; then
    echo "  PyInstaller available via python3 -m PyInstaller"
else
    echo "  PyInstaller not found — installing..."
    python3 -m pip install --upgrade pyinstaller || \
    python3 -m pip install --upgrade pyinstaller --break-system-packages
fi

# Ensure tkinter is available (Homebrew Python needs python-tk)
if ! python3 -c "import tkinter" &>/dev/null; then
    echo ""
    echo "  tkinter not found — attempting: brew install python-tk"
    brew install python-tk@"$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')" || true
fi

# ── Run PyInstaller ───────────────────────────────────────────────────────────
echo ""
echo "[2/4]  Building .app bundle (1-3 minutes)..."
cd "$SCRIPT_DIR"
python3 -m PyInstaller --clean -y protein_liability.spec

# ── Copy .app to dist/ ────────────────────────────────────────────────────────
mkdir -p "$DIST_DIR"
echo ""
echo "[3/4]  Copying $APP_NAME.app to dist/..."
rm -rf "$DIST_DIR/$APP_NAME.app"
cp -R "$SCRIPT_DIR/dist/$APP_NAME.app" "$DIST_DIR/$APP_NAME.app"

# ── Create DMG ────────────────────────────────────────────────────────────────
echo ""
echo "[4/4]  Creating .dmg disk image..."

DMG_NAME="${APP_NAME}_v${VERSION}.dmg"
DMG_PATH="$DIST_DIR/$DMG_NAME"
STAGING="$SCRIPT_DIR/dmg_staging"

rm -f "$DMG_PATH"
rm -rf "$STAGING"
mkdir -p "$STAGING"

cp -R "$DIST_DIR/$APP_NAME.app" "$STAGING/"
ln -s /Applications "$STAGING/Applications"

if command -v create-dmg &>/dev/null; then
    create-dmg \
        --volname "$APP_DISPLAY $VERSION" \
        --window-pos 200 120 \
        --window-size 600 400 \
        --icon-size 128 \
        --icon "$APP_NAME.app" 160 185 \
        --hide-extension "$APP_NAME.app" \
        --app-drop-link 430 185 \
        --no-internet-enable \
        "$DMG_PATH" \
        "$STAGING"
else
    hdiutil create \
        -volname "$APP_DISPLAY $VERSION" \
        -srcfolder "$STAGING" \
        -ov -format UDZO \
        "$DMG_PATH"
fi

rm -rf "$STAGING"

echo ""
echo "======================================================"
echo "  Build complete!"
echo ""
echo "  App bundle : dist/$APP_NAME.app"
echo "  Disk image : dist/$DMG_NAME"
echo ""
echo "  Share the .dmg file. Users open it and drag the app"
echo "  to their Applications folder."
echo ""
echo "  First-launch note: if macOS blocks the app with"
echo "  'unidentified developer', right-click -> Open -> Open."
echo "  (Code-signing requires an Apple Developer account.)"
echo "======================================================"
echo ""
