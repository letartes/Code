# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec — Protein Liability Analyzer  (Auto-HOS Edition)
from pathlib import Path

# Source directory is one level up from this spec file (the Code/ folder)
SRC = str(Path(SPECPATH).parent)

a = Analysis(
    [str(Path(SRC) / 'protein_liability_gui_autohos.py')],
    pathex=[SRC],
    binaries=[],
    datas=[],
    hiddenimports=['protein_liability_analyzer_v2'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'matplotlib', 'numpy', 'pandas', 'scipy', 'PIL',
        'PyQt5', 'PyQt6', 'wx',
    ],
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='ProteinLiabilityAnalyzer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='ProteinLiabilityAnalyzer',
)

app = BUNDLE(
    coll,
    name='ProteinLiabilityAnalyzer.app',
    bundle_identifier='com.letartescientific.proteinliabilityanalyzer',
    info_plist={
        'CFBundleDisplayName':        'Protein Liability Analyzer',
        'CFBundleShortVersionString': '2.0.0',
        'CFBundleVersion':            '2.0.0',
        'NSHighResolutionCapable':     True,
        'LSMinimumSystemVersion':      '10.13.0',
        'NSHumanReadableCopyright':   'Letarte Scientific Consulting',
    },
)
