# -*- mode: python ; coding: utf-8 -*-
# ─────────────────────────────────────────────────────────────────────────────
# PyInstaller spec — Protein Liability Analyzer  (Auto-HOS Edition)
# Works on Windows, macOS, and Linux.
# ─────────────────────────────────────────────────────────────────────────────
import sys
from pathlib import Path

# Resolve the source directory (one level up from this spec file)
SRC = str(Path(SPECPATH).parent)

a = Analysis(
    [str(Path(SRC) / 'protein_liability_gui_autohos.py')],
    pathex=[SRC],
    binaries=[],
    datas=[],
    # protein_liability_analyzer_v2 is a plain Python module in the same
    # directory — PyInstaller discovers it automatically through import analysis.
    hiddenimports=['protein_liability_analyzer_v2'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # Trim unused heavy packages if accidentally pulled in
        'matplotlib', 'numpy', 'pandas', 'scipy', 'PIL',
        'PyQt5', 'PyQt6', 'wx',
    ],
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data)

# ── Single-file executable ────────────────────────────────────────────────────
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='ProteinLiabilityAnalyzer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,          # no console window on Windows
    disable_windowed_traceback=False,
    argv_emulation=False,   # macOS: let Tkinter handle events
    target_arch=None,       # None = native arch; set 'x86_64' or 'arm64' to override
    codesign_identity=None,
    entitlements_file=None,
    # icon='../assets/icon.ico',   # uncomment + supply icon.ico for Windows
)

# ── macOS .app bundle ─────────────────────────────────────────────────────────
# Only generated on macOS; ignored on Windows/Linux.
app = BUNDLE(
    exe,
    name='ProteinLiabilityAnalyzer.app',
    # icon='../assets/icon.icns',  # uncomment + supply icon.icns for macOS
    bundle_identifier='com.letartescientific.proteinliabilityanalyzer',
    info_plist={
        'CFBundleDisplayName':       'Protein Liability Analyzer',
        'CFBundleShortVersionString': '2.0.0',
        'CFBundleVersion':           '2.0.0',
        'NSHighResolutionCapable':    True,
        'LSMinimumSystemVersion':     '10.13.0',
        'NSHumanReadableCopyright':  'Letarte Scientific Consulting',
    },
)
