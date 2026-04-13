@echo off
REM ─────────────────────────────────────────────────────────────────────────────
REM  Protein Liability Analyzer — Windows Build Script
REM  Run this script from inside the  build\  folder, or double-click it.
REM  Requires Python 3.8+ on PATH.  Everything else is installed automatically.
REM ─────────────────────────────────────────────────────────────────────────────

setlocal enabledelayedexpansion

echo.
echo =====================================================
echo  Protein Liability Analyzer — Windows Build
echo  Letarte Scientific Consulting
echo =====================================================
echo.

REM ── Check Python ──────────────────────────────────────────────────────────
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found on PATH.
    echo         Download from https://www.python.org/downloads/
    pause & exit /b 1
)
python --version

REM ── Install / upgrade PyInstaller ─────────────────────────────────────────
echo.
echo [1/3]  Installing PyInstaller...
python -m pip install --upgrade pyinstaller
if errorlevel 1 ( echo [ERROR] pip install failed. & pause & exit /b 1 )

REM ── Run PyInstaller ───────────────────────────────────────────────────────
echo.
echo [2/3]  Building executable (this takes 1-3 minutes)...
cd /d "%~dp0"
python -m PyInstaller --clean protein_liability.spec
if errorlevel 1 ( echo [ERROR] PyInstaller build failed. & pause & exit /b 1 )

REM ── Copy result to dist root ──────────────────────────────────────────────
echo.
echo [3/3]  Copying ProteinLiabilityAnalyzer.exe to dist\ folder...
if not exist "..\dist\" mkdir "..\dist\"
copy /Y "dist\ProteinLiabilityAnalyzer.exe" "..\dist\ProteinLiabilityAnalyzer.exe"

echo.
echo =====================================================
echo  Build complete!
echo.
echo  Executable:  dist\ProteinLiabilityAnalyzer.exe
echo.
echo  To create a full Windows Installer (.exe setup file):
echo    1. Download and install Inno Setup from https://jrsoftware.org/isinfo.php
echo    2. Open  build\installer_windows.iss  in Inno Setup
echo    3. Press F9 (Build) — the Setup.exe will appear in  dist\
echo =====================================================
echo.
pause
