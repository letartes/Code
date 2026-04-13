#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  Protein Liability Analyzer — macOS Build Script
#  Creates a .app bundle and a distributable .dmg disk image.
#
#  Run from a terminal:
#      cd "path/to/sequence_liabilities HOS/build"
#      chmod +x build_mac.sh
#      ./build_mac.sh
#
#  Requires: Python 3.8+ (system or Homebrew), pip
#  Optional: Homebrew  (brew install create-dmg)  for a prettier DMG
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

APP_NAME="ProteinLiabilityAnalyzer"
APP_DISPLAY="Protein Liability Analyzer"
VERSION="2.0"
BUNDLE_ID="com.letartescientific.proteinliabilityanalyzer"

# Resolve directories relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(dirname "$SCRIPT_DIR")"
DIST_DIR="$SRC_DIR/dist"

echo ""
echo "======================================================"
echo "  $APP_DISPLAY — macOS Build"
echo "  Letarte Scientific Consulting"
echo "======================================================"
echo ""

# ── Check Python ──────────────────────────────────────────────────────────────
if ! command -v python3 &>/dev/null; then
    echo "[ERROR] python3 not found. Install from https://www.python.org/downloads/"
    exit 1
fi
echo "Python: $(python3 --version)"

# ── Install / upgrade PyInstaller ─────────────────────────────────────────────
echo ""
echo "[1/4]  Installing PyInstaller..."
python3 -m pip install --upgrade pyinstaller

# ── Run PyInstaller ───────────────────────────────────────────────────────────
echo ""
echo "[2/4]  Building .app bundle (this takes 1-3 minutes)..."
cd "$SCRIPT_DIR"
python3 -m PyInstaller --clean protein_liability.spec

# ── Copy .app to dist/ ────────────────────────────────────────────────────────
mkdir -p "$DIST_DIR"
echo ""
echo "[3/4]  Copying $APP_NAME.app to dist/ ..."
rm -rf "$DIST_DIR/$APP_NAME.app"
cp -R "$SCRIPT_DIR/dist/$APP_NAME.app" "$DIST_DIR/$APP_NAME.app"

# ── Create DMG ────────────────────────────────────────────────────────────────
echo ""
echo "[4/4]  Creating distributable .dmg disk image..."

DMG_NAME="${APP_NAME}_v${VERSION}.dmg"
DMG_PATH="$DIST_DIR/$DMG_NAME"
STAGING="$SCRIPT_DIR/dmg_staging"

# Remove old artifacts
rm -f "$DMG_PATH"
rm -rf "$STAGING"
mkdir -p "$STAGING"

# Copy .app into staging area
cp -R "$DIST_DIR/$APP_NAME.app" "$STAGING/"
# Add a symlink to /Applications for drag-install
ln -s /Applications "$STAGING/Applications"

if command -v create-dmg &>/dev/null; then
    # ── Pretty DMG via create-dmg (Homebrew: brew install create-dmg) ─────────
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
    # ── Plain DMG via hdiutil (always available on macOS) ─────────────────────
    hdiutil create \
        -volname "$APP_DISPLAY $VERSION" \
        -srcfolder "$STAGING" \
        -ov -format UDZO \
        "$DMG_PATH"
fi

# Clean up staging
rm -rf "$STAGING"

echo ""
echo "======================================================"
echo "  Build complete!"
echo ""
echo "  App bundle : dist/$APP_NAME.app"
echo "  Disk image : dist/$DMG_NAME"
echo ""
echo "  To distribute: share $DMG_NAME"
echo "  Users open the DMG and drag the app to Applications."
echo ""
echo "  NOTE: If macOS shows 'unidentified developer' on first"
echo "  launch, right-click the .app → Open → Open anyway."
echo "  (Signing requires an Apple Developer account.)"
echo "======================================================"
echo ""
