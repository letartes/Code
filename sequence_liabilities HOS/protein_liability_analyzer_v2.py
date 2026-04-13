#!/usr/bin/env python3
"""
Protein Sequence Liability Analyzer v2 — Higher Order Structure Extension
==========================================================================
Identifies PTM risks and sequence liabilities in biopharmaceutical protein
sequences.  When no PDB structure is available this tool can optionally:

  • Predict secondary structure locally using the Chou-Fasman algorithm
  • Estimate per-residue Relative Solvent Accessibility (RSA) from sequence
    using a statistical hydrophobicity-context model
  • Report the structural context (secondary structure + exposure class) for
    every detected liability site

When a PDB structure file IS provided the tool:
  • Parses HELIX / SHEET records for accurate secondary-structure assignment
  • Computes per-residue SASA using the Shrake-Rupley algorithm (all-atom,
    probe radius 1.4 Å) — purely local, no internet required
  • Normalises SASA by residue-type maximum ASA to give true RSA values

All computation runs on the local machine.
Sequences are NEVER transmitted to any external server or API.

Detected liabilities (same as v1):
  - Deamidation   : NG (high), NS (high), NA/NT/NQ/QG/QS (medium)
  - Oxidation     : M/W (high), H/C (medium), F/Y (low)
  - Isomerization : DG/DS (high), DT/DA/DN/DD (medium)
  - Glycosylation : N-X-S/T sequon (N-linked), S-P / T-P (O-linked)
  - Truncation    : DP bond, N-term E/Q pyroglutamate, K/R tryptic sites,
                    C-term K clipping, N-term Met

Usage:
  python protein_liability_analyzer_v2.py sequence.fasta
  python protein_liability_analyzer_v2.py --sequence EVQLVESGGGLVQ...
  python protein_liability_analyzer_v2.py --interactive
  python protein_liability_analyzer_v2.py sequence.fasta --hos
  python protein_liability_analyzer_v2.py sequence.fasta --pdb struct.pdb
  python protein_liability_analyzer_v2.py sequence.fasta --pdb struct.pdb --chain A
  python protein_liability_analyzer_v2.py sequence.fasta --output report.html
"""

import re
import sys
import os
import math
import argparse
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Optional


# ══════════════════════════════════════════════════════════════════════════════
# 1. LIABILITY CATALOGUE  (unchanged from v1)
# ══════════════════════════════════════════════════════════════════════════════

LIABILITIES = {
    # ── Deamidation ──────────────────────────────────────────────────────────
    "deamid_high": {
        "label": "Deamidation – High Risk",
        "short": "Deam-H",
        "category": "Deamidation",
        "risk": "high",
        "color_text": "#c0392b",
        "color_bg": "#fde8e8",
        "color_border": "#e74c3c",
        "description": (
            "NG and NS motifs. Asn–Gly is the most susceptible sequence to "
            "deamidation, forming succinimide intermediates rapidly at "
            "physiological pH. Asn–Ser is second most prone."
        ),
        "references": "Geiger & Clarke (1987); Robinson & Robinson (2001)",
        "span_length": 2,
    },
    "deamid_medium": {
        "label": "Deamidation – Medium Risk",
        "short": "Deam-M",
        "category": "Deamidation",
        "risk": "medium",
        "color_text": "#d35400",
        "color_bg": "#fef0e6",
        "color_border": "#e67e22",
        "description": (
            "NA, NT, NQ, NR, NH motifs and Gln-containing NG/NS equivalents "
            "(QG, QS). Glutamine can also deamidate, though more slowly than "
            "asparagine."
        ),
        "references": "Hao et al. (2011); Chelius et al. (2005)",
        "span_length": 2,
    },

    # ── Oxidation ─────────────────────────────────────────────────────────────
    "oxid_high": {
        "label": "Oxidation – High Risk",
        "short": "Ox-H",
        "category": "Oxidation",
        "risk": "high",
        "color_text": "#7d3c98",
        "color_bg": "#f4ecf7",
        "color_border": "#9b59b6",
        "description": (
            "Met (M) and Trp (W). Methionine sulfoxide formation is a primary "
            "oxidative liability in mAbs; tryptophan oxidation (hydroxytryptophan, "
            "kynurenine) commonly observed under photo-oxidative stress."
        ),
        "references": "Lam et al. (1997); Houde et al. (2009)",
        "span_length": 1,
    },
    "oxid_medium": {
        "label": "Oxidation – Medium Risk",
        "short": "Ox-M",
        "category": "Oxidation",
        "risk": "medium",
        "color_text": "#1f618d",
        "color_bg": "#eaf4fb",
        "color_border": "#2980b9",
        "description": (
            "His (H) and Cys (C). Histidine is oxidised to 2-oxo-histidine "
            "under metal-catalysed oxidation conditions. Free cysteines are "
            "prone to disulfide scrambling and oxidation to sulfinic/sulfonic acid."
        ),
        "references": "Uchida & Kawakishi (1993); Liu et al. (2008)",
        "span_length": 1,
    },
    "oxid_low": {
        "label": "Oxidation – Low Risk",
        "short": "Ox-L",
        "category": "Oxidation",
        "risk": "low",
        "color_text": "#117a65",
        "color_bg": "#e8f8f5",
        "color_border": "#1abc9c",
        "description": (
            "Phe (F) and Tyr (Y). These can oxidise under extreme photo- or "
            "metal-catalysed conditions but are generally lower concern."
        ),
        "references": "Stadtman & Levine (2003)",
        "span_length": 1,
    },

    # ── Isomerization ─────────────────────────────────────────────────────────
    "isom_high": {
        "label": "Isomerization – High Risk",
        "short": "Isom-H",
        "category": "Isomerization",
        "risk": "high",
        "color_text": "#1e8449",
        "color_bg": "#eafaf1",
        "color_border": "#27ae60",
        "description": (
            "DG and DS motifs. Asp–Gly undergoes succinimide-mediated "
            "isomerization to isoAsp most readily; Asp–Ser is also highly prone. "
            "These result in backbone insertion of one CH₂ and altered charge."
        ),
        "references": "Clarke (1987); Geiger & Clarke (1987)",
        "span_length": 2,
    },
    "isom_medium": {
        "label": "Isomerization – Medium Risk",
        "short": "Isom-M",
        "category": "Isomerization",
        "risk": "medium",
        "color_text": "#196f3d",
        "color_bg": "#d5f5e3",
        "color_border": "#239b56",
        "description": (
            "DT, DA, DN, DD, DH motifs. These aspartate-containing sequences "
            "can undergo slower isomerization relative to DG/DS."
        ),
        "references": "Stephenson & Clarke (1989)",
        "span_length": 2,
    },

    # ── Glycosylation ─────────────────────────────────────────────────────────
    "n_glycan": {
        "label": "N-Glycosylation Sequon (NXS/T, X≠P)",
        "short": "N-Glycan",
        "category": "Glycosylation",
        "risk": "info",
        "color_text": "#154360",
        "color_bg": "#d6eaf8",
        "color_border": "#2e86c1",
        "description": (
            "Canonical N-linked glycosylation sequon N–X–S/T where X is any "
            "amino acid except Pro. May represent intended or adventitious "
            "glycosylation. Impacts efficacy, PK, and immunogenicity."
        ),
        "references": "Imperiali & Hendrickson (1995); Apweiler et al. (1999)",
        "span_length": 3,
    },
    "o_glycan": {
        "label": "O-Glycosylation Site (S/T-P motif)",
        "short": "O-Glycan",
        "category": "Glycosylation",
        "risk": "info",
        "color_text": "#4a235a",
        "color_bg": "#f5eef8",
        "color_border": "#8e44ad",
        "description": (
            "Ser–Pro and Thr–Pro motifs are characteristic O-glycosylation "
            "sites (mucin-type). Also flags isolated Ser/Thr clusters that may "
            "be O-glycosylated depending on cellular context."
        ),
        "references": "Van den Steen et al. (1998)",
        "span_length": 2,
    },

    # ── Truncation / Clipping ─────────────────────────────────────────────────
    "trunc_high": {
        "label": "Truncation – High Risk",
        "short": "Trunc-H",
        "category": "Truncation",
        "risk": "high",
        "color_text": "#78281f",
        "color_bg": "#fdedec",
        "color_border": "#c0392b",
        "description": (
            "Asp–Pro (DP) bonds are acid-labile and cleave readily under "
            "mildly acidic purification/formulation conditions. "
            "N-terminal Gln (Q) and Glu (E) undergo spontaneous cyclisation to "
            "pyroglutamate, causing loss of one or two mass units and charge "
            "heterogeneity."
        ),
        "references": "Piszkiewicz et al. (1970); Yu et al. (2006)",
        "span_length": 2,
    },
    "trunc_medium": {
        "label": "Truncation – Medium Risk (Tryptic/C-term clipping)",
        "short": "Trunc-M",
        "category": "Truncation",
        "risk": "medium",
        "color_text": "#784212",
        "color_bg": "#fef9e7",
        "color_border": "#f39c12",
        "description": (
            "Lys (K) and Arg (R) are tryptic cleavage sites and are also prone "
            "to C-terminal lysine clipping in CHO-expressed mAbs. "
            "N-terminal Met is subject to Met aminopeptidase removal and oxidation."
        ),
        "references": "Dick et al. (2008); Harris (1995)",
        "span_length": 1,
    },
}

RISK_ORDER = {"high": 0, "medium": 1, "low": 2, "info": 3}


# ══════════════════════════════════════════════════════════════════════════════
# 2. STRUCTURE PREDICTION DATA
# ══════════════════════════════════════════════════════════════════════════════

# Chou-Fasman propensity values (Pa=helix, Pb=sheet, Pt=turn)
# Source: Chou & Fasman (1974) Biochemistry 13(2):211-222
CF_PROPENSITIES = {
    #        Pa     Pb     Pt
    "A": (1.42, 0.83, 0.66),
    "R": (0.98, 0.93, 0.95),
    "N": (0.67, 0.89, 1.56),
    "D": (1.01, 0.54, 1.46),
    "C": (0.70, 1.19, 1.19),
    "Q": (1.11, 1.10, 0.98),
    "E": (1.51, 0.37, 0.74),
    "G": (0.57, 0.75, 1.56),
    "H": (1.00, 0.87, 0.95),
    "I": (1.08, 1.60, 0.47),
    "L": (1.21, 1.30, 0.59),
    "K": (1.16, 0.74, 1.01),
    "M": (1.45, 1.05, 0.60),
    "F": (1.13, 1.38, 0.60),
    "P": (0.57, 0.55, 1.52),
    "S": (0.77, 0.75, 1.43),
    "T": (0.83, 1.19, 0.96),
    "W": (1.08, 1.37, 0.96),
    "Y": (0.69, 1.47, 1.14),
    "V": (1.06, 1.70, 0.50),
}

# Kyte-Doolittle hydrophobicity scale (used for RSA window correction)
KD_SCALE = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}

# Average fractional RSA in coil regions per residue type
# Derived from PDB-based statistics (Rost & Sander 1994; Ahmad et al. 2004)
COIL_AVG_RSA = {
    "A": 0.38, "C": 0.19, "D": 0.49, "E": 0.55, "F": 0.27,
    "G": 0.47, "H": 0.44, "I": 0.21, "K": 0.62, "L": 0.30,
    "M": 0.27, "N": 0.49, "P": 0.45, "Q": 0.52, "R": 0.58,
    "S": 0.45, "T": 0.38, "V": 0.22, "W": 0.30, "Y": 0.40,
}

# Multiplier applied to coil-based RSA estimate for each SS class
# (Residues in helices/sheets are on average more buried than coil residues)
SS_RSA_MULT = {"H": 0.78, "E": 0.72, "T": 0.95, "C": 1.00}

# Maximum ASA (Å²) per residue type — Tien et al. (2013) PLOS ONE
MAX_ASA = {
    "A": 121.0, "R": 265.0, "N": 187.0, "D": 187.0, "C": 148.0,
    "Q": 214.0, "E": 214.0, "G":  97.0, "H": 216.0, "I": 195.0,
    "L": 191.0, "K": 230.0, "M": 203.0, "F": 228.0, "P": 154.0,
    "S": 143.0, "T": 163.0, "W": 264.0, "Y": 255.0, "V": 165.0,
}

# Van der Waals radii (Å) for heavy atoms used in SASA calculation
VDW_RADII = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "H": 1.20, "P": 1.80, "SE": 1.90, "FE": 1.47,
    "ZN": 1.39, "MG": 1.73,
}

# 3-letter to 1-letter amino acid code (includes common modified residues)
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Modified amino acids
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H",
    "HIE": "H", "HIP": "H", "CSD": "C", "CSX": "C", "CYX": "C",
    "LYZ": "K", "GLH": "E", "ASH": "D", "SEP": "S", "TPO": "T",
    "PTR": "Y", "MLY": "K",
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. SEQUENCE PARSING  (same as v1)
# ══════════════════════════════════════════════════════════════════════════════

def parse_fasta(text: str) -> list:
    """Return list of (name, sequence) tuples from FASTA text."""
    sequences = []
    current_name = "Sequence"
    current_seq = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_seq:
                sequences.append((current_name, "".join(current_seq).upper()))
            current_name = line[1:].strip() or "Sequence"
            current_seq = []
        else:
            current_seq.append(re.sub(r"[^A-Za-z]", "", line))
    if current_seq:
        sequences.append((current_name, "".join(current_seq).upper()))
    return sequences


def clean_sequence(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper()


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def validate_sequence(seq: str):
    warnings = []
    invalid = sorted(set(seq) - VALID_AA)
    if invalid:
        warnings.append(f"Non-standard characters found and ignored: {', '.join(invalid)}")
    if len(seq) < 5:
        warnings.append("Sequence is very short (< 5 residues).")
    return len(invalid) == 0, warnings


# ══════════════════════════════════════════════════════════════════════════════
# 4. LIABILITY DETECTION  (same as v1)
# ══════════════════════════════════════════════════════════════════════════════

def find_liabilities(seq: str) -> list:
    findings = []
    n = len(seq)

    def add(liability_key: str, pos0: int, span: int, note: str = ""):
        findings.append({
            "type": liability_key,
            "pos0": pos0,
            "span": span,
            "pos1": pos0 + 1,
            "residues": seq[pos0: pos0 + span],
            "note": note,
            **{k: v for k, v in LIABILITIES[liability_key].items()
               if k not in ("description", "references", "span_length")},
        })

    HIGH_DEAMID = {"NG", "NS"}
    MED_DEAMID  = {"NA", "NT", "NQ", "NR", "NH", "QG", "QS"}
    for i in range(n - 1):
        pair = seq[i:i+2]
        if pair in HIGH_DEAMID:
            add("deamid_high", i, 2)
        elif pair in MED_DEAMID:
            add("deamid_medium", i, 2)

    for i, aa in enumerate(seq):
        if aa in "MW":
            note = "N-terminal Met (also truncation risk)" if aa == "M" and i == 0 else ""
            add("oxid_high", i, 1, note)
        elif aa in "HC":
            add("oxid_medium", i, 1)
        elif aa in "FY":
            add("oxid_low", i, 1)

    HIGH_ISOM = {"DG", "DS"}
    MED_ISOM  = {"DT", "DA", "DN", "DD", "DH"}
    for i in range(n - 1):
        pair = seq[i:i+2]
        if pair in HIGH_ISOM:
            add("isom_high", i, 2)
        elif pair in MED_ISOM:
            add("isom_medium", i, 2)

    for i in range(n - 2):
        if seq[i] == "N" and seq[i+1] != "P" and seq[i+2] in "ST":
            add("n_glycan", i, 3)

    for i in range(n - 1):
        if seq[i] in "ST" and seq[i+1] == "P":
            add("o_glycan", i, 2)

    for i in range(n - 1):
        if seq[i:i+2] == "DP":
            add("trunc_high", i, 2, "Acid-labile Asp–Pro bond")

    if seq and seq[0] in "QE":
        note = (
            "N-terminal Gln → pyroglutamate (–17 Da)"
            if seq[0] == "Q"
            else "N-terminal Glu → pyroglutamate (–18 Da)"
        )
        add("trunc_high", 0, 1, note)

    for i, aa in enumerate(seq):
        if aa in "KR":
            if i == n - 1 and aa == "K":
                add("trunc_medium", i, 1, "C-terminal Lys clipping (CHO expression)")
            else:
                next_aa = seq[i+1] if i + 1 < n else ""
                if next_aa != "P":
                    add("trunc_medium", i, 1, f"Tryptic site ({aa})")
                else:
                    add("trunc_medium", i, 1, f"Trypsin-resistant {aa}–Pro (low risk)")

    if seq and seq[0] == "M":
        add("trunc_medium", 0, 1, "N-terminal Met: subject to Met aminopeptidase removal")

    return findings


# ══════════════════════════════════════════════════════════════════════════════
# 5. SECONDARY STRUCTURE PREDICTION — Chou-Fasman algorithm
# ══════════════════════════════════════════════════════════════════════════════

def predict_secondary_structure(seq: str) -> list:
    """
    Predict secondary structure using the Chou-Fasman (1974) algorithm.

    Returns a list of characters, one per residue:
      'H' = alpha-helix   'E' = beta-strand   'T' = beta-turn   'C' = coil

    Algorithm steps:
      1. Assign propensity scores (Pa, Pb, Pt) to each residue.
      2. Find helix nucleation windows (6-mer with ≥4 residues having Pa≥1.0).
      3. Find sheet nucleation windows (5-mer with ≥3 residues having Pb≥1.0).
      4. Extend each nucleus while the 4-residue sliding average stays ≥1.0.
      5. Identify turns in unassigned 4-mers where Pt > Pa and Pt > Pb.
      6. Resolve H/E conflicts: whichever has higher propensity at that site wins.
    """
    n = len(seq)
    if n < 4:
        return ["C"] * n

    pa = [CF_PROPENSITIES.get(aa, (1.0, 1.0, 1.0))[0] for aa in seq]
    pb = [CF_PROPENSITIES.get(aa, (1.0, 1.0, 1.0))[1] for aa in seq]
    pt = [CF_PROPENSITIES.get(aa, (1.0, 1.0, 1.0))[2] for aa in seq]

    # --- Helix nucleation and extension --------------------------------------
    helix = [False] * n
    i = 0
    while i <= n - 6:
        window = pa[i: i + 6]
        if sum(1 for p in window if p >= 1.0) >= 4:
            left, right = i, i + 6
            # Extend right
            while right < n and sum(pa[max(0, right - 3): right + 1]) / 4 >= 1.0:
                right += 1
            # Extend left
            while left > 0 and sum(pa[left - 1: left + 3]) / 4 >= 1.0:
                left -= 1
            for j in range(left, right):
                helix[j] = True
            i = right
        else:
            i += 1

    # --- Sheet nucleation and extension --------------------------------------
    sheet = [False] * n
    i = 0
    while i <= n - 5:
        window = pb[i: i + 5]
        if sum(1 for p in window if p >= 1.0) >= 3:
            left, right = i, i + 5
            while right < n and sum(pb[max(0, right - 3): right + 1]) / 4 >= 1.0:
                right += 1
            while left > 0 and sum(pb[left - 1: left + 3]) / 4 >= 1.0:
                left -= 1
            for j in range(left, right):
                sheet[j] = True
            i = right
        else:
            i += 1

    # --- Turn detection in unassigned coil regions ---------------------------
    turn = [False] * n
    for i in range(n - 3):
        if all(not helix[j] and not sheet[j] for j in range(i, i + 4)):
            avg_pt = sum(pt[i: i + 4]) / 4
            avg_pa = sum(pa[i: i + 4]) / 4
            avg_pb = sum(pb[i: i + 4]) / 4
            if avg_pt > 1.0 and avg_pt > avg_pa and avg_pt > avg_pb:
                for j in range(i, i + 4):
                    turn[j] = True

    # --- Final assignment ----------------------------------------------------
    ss = ["C"] * n
    for i in range(n):
        if helix[i] and sheet[i]:
            ss[i] = "H" if pa[i] >= pb[i] else "E"
        elif helix[i]:
            ss[i] = "H"
        elif sheet[i]:
            ss[i] = "E"
        elif turn[i]:
            ss[i] = "T"

    return ss


# ══════════════════════════════════════════════════════════════════════════════
# 6. RSA PREDICTION FROM SEQUENCE
# ══════════════════════════════════════════════════════════════════════════════

def predict_rsa_from_sequence(seq: str, ss: list) -> list:
    """
    Estimate per-residue Relative Solvent Accessibility (RSA) from sequence.

    Method:
      1. Assign a base RSA from the residue type's average in coil regions
         (derived from PDB statistics).
      2. Multiply by an SS-class modifier:
            Helix ×0.78,  Strand ×0.72,  Turn ×0.95,  Coil ×1.00
      3. Apply a sliding-window (11-mer) hydrophobicity correction:
            hydrophobic surroundings → bury the residue (lower RSA)
            hydrophilic surroundings → expose the residue (higher RSA)
      4. Boost terminal residues slightly (they are often more flexible/exposed).

    Returns list of floats in [0.02, 0.98].
    Note: This is a PREDICTION from sequence alone.  Accuracy improves with
    a real PDB structure (use --pdb).
    """
    n = len(seq)
    WINDOW = 11
    half   = WINDOW // 2

    rsa = []
    for i, (aa, s) in enumerate(zip(seq, ss)):
        base = COIL_AVG_RSA.get(aa, 0.40)
        mult = SS_RSA_MULT.get(s, 1.00)
        rsa_i = base * mult

        # Hydrophobicity window correction
        start = max(0, i - half)
        end   = min(n, i + half + 1)
        win   = seq[start:end]
        avg_kd = sum(KD_SCALE.get(a, 0.0) for a in win) / len(win)
        correction = -0.09 * (avg_kd / 4.5)   # KD range ≈ ±4.5
        rsa_i = max(0.02, min(0.98, rsa_i + correction))
        rsa.append(rsa_i)

    # Terminal boost (N-terminal & C-terminal 5 residues)
    for j in range(min(5, n)):
        rsa[j]       = min(0.98, rsa[j]       * 1.15)
        rsa[n - 1 - j] = min(0.98, rsa[n - 1 - j] * 1.15)

    return rsa


# ══════════════════════════════════════════════════════════════════════════════
# 7. PDB PARSING
# ══════════════════════════════════════════════════════════════════════════════

def parse_pdb_atoms(pdb_text: str, chain_id: Optional[str] = None) -> list:
    """
    Parse ATOM records from PDB text.
    Returns a list of dicts: {name, chain, res_num, res_name, x, y, z, element}
    Hydrogen atoms are excluded.  If chain_id is given, only that chain is returned.
    """
    atoms = []
    seen_res = {}   # (chain, res_num) → res_name  (deduplicate alt-locs)

    for line in pdb_text.splitlines():
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        try:
            atom_name = line[12:16].strip()
            alt_loc   = line[16].strip()
            res_name  = line[17:20].strip()
            chain     = line[21].strip() if len(line) > 21 else "A"
            res_num   = int(line[22:26].strip())
            x         = float(line[30:38])
            y         = float(line[38:46])
            z         = float(line[46:54])
            # Element from cols 76-78 if present; else infer from atom name
            element   = line[76:78].strip().upper() if len(line) > 77 else ""
            if not element:
                element = re.sub(r"[^A-Z]", "", atom_name[0:2]).upper()
                element = element[0] if element else "C"
        except (ValueError, IndexError):
            continue

        # Skip alternate locations other than A or blank
        if alt_loc not in ("", "A"):
            continue
        # Skip hydrogens
        if element in ("H", "D"):
            continue
        # Skip water
        if res_name in ("HOH", "WAT", "H2O"):
            continue

        if chain_id and chain != chain_id:
            continue

        atoms.append({
            "name":     atom_name,
            "chain":    chain,
            "res_num":  res_num,
            "res_name": res_name,
            "x": x, "y": y, "z": z,
            "element":  element,
        })

    return atoms


def parse_pdb_helix_sheet(pdb_text: str):
    """
    Parse HELIX and SHEET records to build a secondary-structure lookup.
    Returns dict: (chain, res_num) → 'H' or 'E'
    """
    ss_map = {}

    for line in pdb_text.splitlines():
        if line.startswith("HELIX "):
            try:
                chain_i = line[19].strip()
                start_i = int(line[21:25].strip())
                chain_j = line[31].strip()
                end_j   = int(line[33:37].strip())
                chain   = chain_i if chain_i else chain_j
                for r in range(start_i, end_j + 1):
                    ss_map[(chain, r)] = "H"
            except (ValueError, IndexError):
                pass

        elif line.startswith("SHEET "):
            try:
                chain_i = line[21].strip()
                start_i = int(line[22:26].strip())
                chain_j = line[32].strip()
                end_j   = int(line[33:37].strip())
                chain   = chain_i if chain_i else chain_j
                for r in range(start_i, end_j + 1):
                    ss_map[(chain, r)] = "E"
            except (ValueError, IndexError):
                pass

    return ss_map


def extract_pdb_residues(atoms: list) -> list:
    """
    Return ordered list of (chain, res_num, res_name) without duplicates.
    """
    seen = set()
    residues = []
    for a in atoms:
        key = (a["chain"], a["res_num"])
        if key not in seen:
            seen.add(key)
            residues.append((a["chain"], a["res_num"], a["res_name"]))
    return residues


def map_pdb_to_sequence(pdb_residues: list, seq: str, chain_id: Optional[str] = None):
    """
    Attempt to map PDB residue numbers to 0-based sequence positions.

    Converts 3-letter residue names to 1-letter and tries a substring match
    or a global offset alignment.  Returns dict: seq_pos → (chain, res_num).
    Falls back to positional mapping (1-to-1) if sequences differ.
    """
    # Build PDB sequence from ATOM records
    pdb_aa = []
    for chain, res_num, res_name in pdb_residues:
        aa = AA3TO1.get(res_name.upper(), "X")
        pdb_aa.append((chain, res_num, aa))

    pdb_seq = "".join(a[2] for a in pdb_aa)

    # Try to find seq inside pdb_seq (or vice versa)
    mapping = {}
    best_offset = 0
    best_score  = -1

    # Try every start offset within ±50 residues
    for offset in range(-50, 51):
        score = 0
        for i, aa in enumerate(seq):
            pdb_i = i + offset
            if 0 <= pdb_i < len(pdb_seq) and pdb_seq[pdb_i] == aa:
                score += 1
        if score > best_score:
            best_score  = score
            best_offset = offset

    for i in range(len(seq)):
        pdb_i = i + best_offset
        if 0 <= pdb_i < len(pdb_aa):
            mapping[i] = (pdb_aa[pdb_i][0], pdb_aa[pdb_i][1])

    match_pct = best_score / len(seq) * 100 if seq else 0
    return mapping, match_pct


# ══════════════════════════════════════════════════════════════════════════════
# 8. SHRAKE-RUPLEY SASA CALCULATION  (pure Python, local only)
# ══════════════════════════════════════════════════════════════════════════════

def _fibonacci_sphere(n: int) -> list:
    """Generate n uniformly distributed points on the unit sphere (Fibonacci lattice)."""
    pts  = []
    gold = math.pi * (3.0 - math.sqrt(5.0))   # golden angle ≈ 2.3999 rad
    for i in range(n):
        y     = 1.0 - (i / max(n - 1, 1)) * 2.0
        r     = math.sqrt(max(0.0, 1.0 - y * y))
        theta = gold * i
        pts.append((r * math.cos(theta), y, r * math.sin(theta)))
    return pts


# Pre-compute sphere points at module load time
_SPHERE_92  = _fibonacci_sphere(92)
_SPHERE_194 = _fibonacci_sphere(194)   # higher accuracy option


def compute_sasa_per_residue(atoms: list, n_sphere_points: int = 92,
                              probe: float = 1.4,
                              progress: bool = False) -> dict:
    """
    Compute per-residue Solvent-Accessible Surface Area (Å²) using the
    Shrake-Rupley (1973) rolling-probe algorithm.

    Algorithm:
      For each atom i:
        1. Place n_sphere_points uniformly on the sphere of radius (VDW_r + probe).
        2. For each test point, check if it lies inside any neighbouring atom's
           extended sphere (VDW_r_j + probe).
        3. ASA(i) = (accessible_points / n_sphere_points) × 4π(VDW_r + probe)²

    Neighbour search is accelerated with a 3-D grid (cell size = 6 Å).
    Only heavy atoms (no H) should be in the input list.

    Returns {res_num: SASA_in_Å²}
    """
    sphere_pts = _SPHERE_194 if n_sphere_points > 92 else _SPHERE_92
    n_pts = len(sphere_pts)

    # Build spatial grid for fast neighbour lookup
    CELL = 6.0   # Å  (larger than VDW + probe + VDW + probe ≈ 5 Å)
    grid: dict = defaultdict(list)
    for idx, a in enumerate(atoms):
        key = (int(a["x"] / CELL), int(a["y"] / CELL), int(a["z"] / CELL))
        grid[key].append(idx)

    res_sasa: dict = defaultdict(float)
    n_atoms = len(atoms)

    for i, atom in enumerate(atoms):
        if progress and i % 500 == 0:
            pct = i / n_atoms * 100
            print(f"\r  Computing SASA … {pct:5.1f}%", end="", flush=True)

        r_i = VDW_RADII.get(atom["element"], 1.70) + probe
        ax, ay, az = atom["x"], atom["y"], atom["z"]

        cx = int(ax / CELL)
        cy = int(ay / CELL)
        cz = int(az / CELL)

        # Collect neighbours (atoms whose extended spheres could overlap)
        neighbours = []
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                for dz in range(-2, 3):
                    for j in grid.get((cx + dx, cy + dy, cz + dz), []):
                        if j == i:
                            continue
                        nb = atoms[j]
                        r_j = VDW_RADII.get(nb["element"], 1.70) + probe
                        ddx = ax - nb["x"]
                        ddy = ay - nb["y"]
                        ddz = az - nb["z"]
                        dist2 = ddx * ddx + ddy * ddy + ddz * ddz
                        if dist2 < (r_i + r_j) ** 2:
                            neighbours.append((nb["x"], nb["y"], nb["z"], r_j))

        # Count accessible probe positions
        accessible = 0
        for sp in sphere_pts:
            tx = ax + r_i * sp[0]
            ty = ay + r_i * sp[1]
            tz = az + r_i * sp[2]
            buried = False
            for nx, ny, nz, r_j in neighbours:
                dx2 = tx - nx
                dy2 = ty - ny
                dz2 = tz - nz
                if dx2 * dx2 + dy2 * dy2 + dz2 * dz2 < r_j * r_j:
                    buried = True
                    break
            if not buried:
                accessible += 1

        atom_sasa = (accessible / n_pts) * 4.0 * math.pi * r_i * r_i
        res_sasa[atom["res_num"]] += atom_sasa

    if progress:
        print("\r  Computing SASA … done.       ")

    return dict(res_sasa)


# ══════════════════════════════════════════════════════════════════════════════
# 9. STRUCTURAL CONTEXT — annotate findings with SS and RSA
# ══════════════════════════════════════════════════════════════════════════════

def add_structural_context(findings: list, ss: list, rsa: list,
                            source: str = "predicted (Chou-Fasman)") -> None:
    """
    Annotate each finding dict in-place with:
      ss            : secondary structure at the most-exposed residue in span
      rsa           : max RSA within the span (the most accessible residue)
      rsa_source    : 'predicted (Chou-Fasman)' or 'PDB (Shrake-Rupley SASA)'
      exposure_class: 'Exposed', 'Partial', or 'Buried'
      exposure_note : brief implication for liability reactivity
    """
    n = len(ss)
    for f in findings:
        pos  = f["pos0"]
        span = f["span"]
        end  = min(pos + span, n)

        if pos >= n:
            f["ss"] = "C"; f["rsa"] = 0.40; f["rsa_source"] = source
            f["exposure_class"] = "Unknown"; f["exposure_note"] = "—"
            continue

        # For SS: use the first residue of the liability span
        f["ss"] = ss[pos]

        # For RSA: report the maximum within the span (most exposed residue)
        span_rsa = rsa[pos:end]
        f["rsa"] = round(max(span_rsa) if span_rsa else 0.40, 3)
        f["rsa_source"] = source

        rsa_val = f["rsa"]
        if rsa_val > 0.25:
            f["exposure_class"] = "Exposed"
            f["exposure_note"]  = "Surface-exposed — elevated chemical reactivity"
        elif rsa_val > 0.10:
            f["exposure_class"] = "Partial"
            f["exposure_note"]  = "Partially accessible — moderate reactivity"
        else:
            f["exposure_class"] = "Buried"
            f["exposure_note"]  = "Buried — reduced chemical accessibility"


# ══════════════════════════════════════════════════════════════════════════════
# 10. STATISTICS
# ══════════════════════════════════════════════════════════════════════════════

def summarize(findings: list, seq_len: int) -> dict:
    by_category = defaultdict(list)
    by_risk     = defaultdict(list)
    by_type     = defaultdict(list)
    for f in findings:
        by_category[f["category"]].append(f)
        by_risk[f["risk"]].append(f)
        by_type[f["type"]].append(f)
    return {
        "total":       len(findings),
        "by_category": dict(by_category),
        "by_risk":     dict(by_risk),
        "by_type":     dict(by_type),
        "seq_len":     seq_len,
    }


def hos_stats(ss: list, rsa: list) -> dict:
    """Summary statistics for the HOS analysis."""
    n = len(ss)
    counts = {"H": 0, "E": 0, "T": 0, "C": 0}
    for s in ss:
        counts[s] = counts.get(s, 0) + 1
    avg_rsa = sum(rsa) / len(rsa) if rsa else 0.0
    exposed  = sum(1 for r in rsa if r > 0.25)
    buried   = sum(1 for r in rsa if r <= 0.10)
    return {
        "n":       n,
        "counts":  counts,
        "avg_rsa": round(avg_rsa, 3),
        "exposed": exposed,
        "buried":  buried,
    }


# ══════════════════════════════════════════════════════════════════════════════
# 11. HTML REPORT
# ══════════════════════════════════════════════════════════════════════════════

CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  background: #f0f2f5;
  color: #1a1a2e;
  font-size: 15px;
}
.page-wrap { max-width: 1200px; margin: 0 auto; padding: 24px 16px 60px; }

/* Header */
.report-header {
  background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
  color: #fff;
  border-radius: 12px;
  padding: 32px 36px;
  margin-bottom: 28px;
}
.report-header h1 { font-size: 1.9em; font-weight: 700; margin-bottom: 6px; }
.report-header .subtitle { color: #a8b2c8; font-size: 0.95em; }
.report-header .meta { margin-top: 16px; font-size: 0.85em; color: #c0c8dc;
  display:flex; gap:24px; flex-wrap:wrap; }
.report-header .meta span { display:flex; align-items:center; gap:6px; }

/* Summary cards grid */
.summary-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 14px;
  margin-bottom: 28px;
}
.summary-card {
  background: #fff;
  border-radius: 10px;
  padding: 20px 16px;
  box-shadow: 0 1px 6px rgba(0,0,0,.08);
  text-align: center;
}
.card-icon  { font-size: 1.8em; margin-bottom: 6px; }
.card-title { font-weight: 600; font-size: 0.9em; color: #555; margin-bottom: 4px; }
.card-count { font-size: 2em; font-weight: 700; color: #1a1a2e; margin-bottom: 8px; }
.card-bars  { display:flex; flex-wrap:wrap; gap:4px; justify-content:center; }

/* Risk badges */
.risk-badge {
  display: inline-block;
  padding: 2px 7px;
  border-radius: 20px;
  font-size: 0.72em;
  font-weight: 700;
  letter-spacing: .4px;
  white-space: nowrap;
}
.risk-high   { background:#fde8e8; color:#c0392b; border:1px solid #e74c3c; }
.risk-medium { background:#fef0e6; color:#d35400; border:1px solid #e67e22; }
.risk-low    { background:#e8f8f5; color:#117a65; border:1px solid #1abc9c; }
.risk-info   { background:#d6eaf8; color:#154360; border:1px solid #2e86c1; }

/* Sections */
.section {
  background: #fff;
  border-radius: 10px;
  box-shadow: 0 1px 6px rgba(0,0,0,.08);
  margin-bottom: 24px;
  overflow: hidden;
}
.section-header {
  padding: 16px 24px;
  border-bottom: 1px solid #eee;
  font-weight: 700;
  font-size: 1.05em;
  color: #1a1a2e;
  background: #fafbfc;
}
.section-body { padding: 20px 24px; }

/* Sequence display */
.seq-block {
  font-family: "Courier New", "Courier", monospace;
  font-size: 0.90em;
  line-height: 1.9;
  letter-spacing: 0.05em;
  background: #f8f9fa;
  border-radius: 8px;
  padding: 16px;
  overflow-x: auto;
  white-space: pre;
}
.seq-line    { display:block; white-space:pre; }
.ss-line     { display:block; white-space:pre; line-height:1.4; }
.rsa-line    { display:block; white-space:pre; line-height:1.4; margin-bottom:6px; }
.aa          { display:inline-block; }
.liability   { border-radius: 2px; cursor: help; }
.pos-label   { color: #9ba5b0; user-select: none; font-size: 0.88em; }
.ss-H { background:#00b4d8; color:#fff; font-size:0.7em; padding:0 1px; border-radius:1px; }
.ss-E { background:#f9c74f; color:#333; font-size:0.7em; padding:0 1px; border-radius:1px; }
.ss-T { background:#90be6d; color:#fff; font-size:0.7em; padding:0 1px; border-radius:1px; }
.ss-C { background:#dee2e6; color:#888; font-size:0.7em; padding:0 1px; border-radius:1px; }

/* HOS stat bar */
.hos-bar-wrap {
  display:flex; gap:4px; align-items:center; margin-top:10px; flex-wrap:wrap;
}
.hos-bar-seg {
  height:18px;
  border-radius:3px;
  display:inline-flex;
  align-items:center;
  justify-content:center;
  font-size:0.72em;
  font-weight:700;
  color:#fff;
  min-width:30px;
  white-space:nowrap;
  padding:0 4px;
}

/* Exposure badges */
.exp-exposed { background:#ffe0e0; color:#c0392b; border:1px solid #e74c3c;
               border-radius:10px; padding:1px 7px; font-size:0.75em; font-weight:700; }
.exp-partial { background:#fef5e7; color:#d35400; border:1px solid #f39c12;
               border-radius:10px; padding:1px 7px; font-size:0.75em; font-weight:700; }
.exp-buried  { background:#e8f4fd; color:#1565c0; border:1px solid #2980b9;
               border-radius:10px; padding:1px 7px; font-size:0.75em; font-weight:700; }

/* Findings table */
table { width: 100%; border-collapse: collapse; font-size: 0.9em; }
thead th {
  background: #f4f5f7;
  padding: 10px 12px;
  text-align: left;
  font-weight: 600;
  color: #444;
  border-bottom: 2px solid #e0e3e8;
  white-space: nowrap;
}
tbody tr { border-bottom: 1px solid #f0f0f0; }
tbody tr:hover { background: #fafbff; }
tbody td { padding: 8px 12px; vertical-align: middle; }
.pos-cell { font-family: monospace; color: #7f8c8d; font-weight: 600; }

/* Legend table */
.legend-table { font-size: 0.85em; }
.legend-table td { padding: 6px 10px; vertical-align: top; }
.legend-table tr { border-bottom: 1px solid #f0f0f0; }

/* Info box */
.info-box {
  border-left: 4px solid #2980b9;
  background: #f0f7fd;
  padding: 12px 16px;
  border-radius: 0 8px 8px 0;
  font-size: 0.88em;
  color: #2c3e50;
  margin-bottom: 18px;
  line-height: 1.6;
}

/* Footer */
.footer {
  text-align: center;
  margin-top: 40px;
  font-size: 0.82em;
  color: #aaa;
}

@media print {
  body { background: #fff; }
  .section { box-shadow: none; border: 1px solid #ddd; }
}
"""

SS_COLORS = {
    "H": ("#00b4d8", "#fff"),   # helix: teal / white text
    "E": ("#f9c74f", "#333"),   # sheet: gold / dark text
    "T": ("#90be6d", "#fff"),   # turn:  sage / white text
    "C": ("#dee2e6", "#888"),   # coil:  light grey / grey text
}

SS_LABELS = {"H": "H", "E": "E", "T": "T", "C": "·"}


def _rsa_color(rsa: float) -> str:
    """Map RSA 0→1 to a color: blue (buried) → white → red (exposed)."""
    if rsa <= 0.25:
        t = rsa / 0.25
        r = int(70  + 185 * t)
        g = int(100 + 155 * t)
        b = int(200 -  45 * t)
    else:
        t = (rsa - 0.25) / 0.75
        r = 255
        g = int(255 - 180 * t)
        b = int(155 - 155 * t)
    return f"rgb({r},{g},{b})"


def build_annotated_sequence(seq: str, findings: list,
                              ss: Optional[list] = None,
                              rsa: Optional[list] = None,
                              wrap: int = 60) -> str:
    """
    Build the HTML annotated-sequence block.
    When ss/rsa are provided, adds two extra tracks below each sequence line:
      - Secondary-structure colour track (H/E/T/C)
      - RSA heatmap track (blue=buried → red=exposed)
    """
    n = len(seq)
    has_struct = ss is not None and rsa is not None

    # Build per-position liability mapping
    pos_findings = [[] for _ in range(n)]
    for f in findings:
        for j in range(f["pos0"], min(f["pos0"] + f["span"], n)):
            pos_findings[j].append(f)

    def best_finding(flist):
        if not flist:
            return None
        return min(flist, key=lambda f: RISK_ORDER.get(f["risk"], 9))

    # Build HTML spans for each residue
    html_chars = []
    for i, aa in enumerate(seq):
        flist = pos_findings[i]
        if not flist:
            html_chars.append(f'<span class="aa">{aa}</span>')
        else:
            bf    = best_finding(flist)
            types = " | ".join(sorted({f["short"] for f in flist}))
            notes = "; ".join(filter(None, [f.get("note", "") for f in flist]))
            title = f"Pos {i+1}: {types}"
            if notes:
                title += f" — {notes}"
            style = (
                f"background:{bf['color_bg']};"
                f"color:{bf['color_text']};"
                f"border-bottom:2px solid {bf['color_border']};"
            )
            html_chars.append(
                f'<span class="aa liability" style="{style}" title="{title}">{aa}</span>'
            )

    lines_html = []
    for start in range(0, n, wrap):
        chunk    = html_chars[start: start + wrap]
        end_pos  = min(start + wrap, n)
        pos_lbl  = f'<span class="pos-label">{start+1:>5} </span>'
        end_lbl  = f' <span class="pos-label">{end_pos}</span>'

        seq_line = f'<span class="seq-line">{pos_lbl}{"".join(chunk)}{end_lbl}</span>'
        lines_html.append(seq_line)

        if has_struct:
            # Secondary-structure track
            ss_spans = []
            for i in range(start, end_pos):
                s = ss[i] if i < len(ss) else "C"
                bg, fg = SS_COLORS.get(s, ("#dee2e6", "#888"))
                lbl = SS_LABELS.get(s, "·")
                ss_spans.append(
                    f'<span class="aa ss-{s}" style="background:{bg};color:{fg};"'
                    f' title="Pos {i+1}: {s}">{lbl}</span>'
                )
            ss_line = (
                f'<span class="ss-line">'
                f'<span class="pos-label">{"":>6}</span>'
                f'{"".join(ss_spans)}'
                f'</span>'
            )
            lines_html.append(ss_line)

            # RSA heatmap track
            rsa_spans = []
            for i in range(start, end_pos):
                r = rsa[i] if i < len(rsa) else 0.4
                col  = _rsa_color(r)
                txt_col = "#333" if r > 0.4 else "#ccc"
                lbl  = "█"
                rsa_spans.append(
                    f'<span class="aa" style="background:{col};color:{col};"'
                    f' title="Pos {i+1}: RSA={r:.2f}">{lbl}</span>'
                )
            rsa_line = (
                f'<span class="rsa-line">'
                f'<span class="pos-label">{"":>6}</span>'
                f'{"".join(rsa_spans)}'
                f'</span>'
            )
            lines_html.append(rsa_line)

    return "\n".join(lines_html)


def liability_legend_html() -> str:
    rows = []
    for key, lib in LIABILITIES.items():
        risk_badge = (
            f'<span class="risk-badge risk-{lib["risk"]}">'
            f'{lib["risk"].upper()}</span>'
        )
        rows.append(
            f'<tr>'
            f'<td><span style="background:{lib["color_bg"]};color:{lib["color_text"]};'
            f'border:1px solid {lib["color_border"]};padding:2px 6px;border-radius:3px;'
            f'font-family:monospace;font-weight:600">{lib["short"]}</span></td>'
            f'<td>{lib["label"]}</td>'
            f'<td>{risk_badge}</td>'
            f'<td style="color:#555;font-size:0.85em">{lib["description"]}</td>'
            f'</tr>'
        )
    return "\n".join(rows)


def findings_table_html(findings: list, has_structure: bool = False) -> str:
    if not findings:
        return "<p><em>No liabilities found.</em></p>"

    sorted_f = sorted(findings, key=lambda f: (RISK_ORDER.get(f["risk"], 9), f["pos0"]))

    rows = []
    for f in sorted_f:
        risk_badge = (
            f'<span class="risk-badge risk-{f["risk"]}">{f["risk"].upper()}</span>'
        )
        res_span = (
            f'<span style="background:{f["color_bg"]};color:{f["color_text"]};'
            f'border:1px solid {f["color_border"]};padding:2px 5px;'
            f'border-radius:3px;font-family:monospace;font-weight:700">'
            f'{f["residues"]}</span>'
        )
        note = f["note"] if f.get("note") else "—"

        struct_cols = ""
        if has_structure:
            ss_val  = f.get("ss", "?")
            rsa_val = f.get("rsa", None)
            exp_cls = f.get("exposure_class", "—")
            exp_note= f.get("exposure_note", "—")

            # SS badge
            bg, fg = SS_COLORS.get(ss_val, ("#dee2e6", "#888"))
            ss_full = {"H": "Helix", "E": "Strand", "T": "Turn", "C": "Coil"}.get(ss_val, ss_val)
            ss_badge = (
                f'<span style="background:{bg};color:{fg};padding:2px 7px;'
                f'border-radius:3px;font-family:monospace;font-weight:700;font-size:0.85em">'
                f'{ss_val} <span style="font-weight:400;font-size:0.9em">{ss_full}</span></span>'
            )

            # RSA bar
            if rsa_val is not None:
                bar_w   = int(rsa_val * 60)
                bar_col = _rsa_color(rsa_val)
                rsa_html = (
                    f'<div style="display:flex;align-items:center;gap:6px;">'
                    f'<div style="width:60px;height:8px;background:#eee;border-radius:4px;">'
                    f'<div style="width:{bar_w}px;height:8px;background:{bar_col};border-radius:4px;"></div>'
                    f'</div>'
                    f'<span style="font-family:monospace;font-size:0.88em">{rsa_val:.2f}</span>'
                    f'</div>'
                )
            else:
                rsa_html = "—"

            # Exposure class badge
            cls_lower = exp_cls.lower()
            exp_badge = f'<span class="exp-{cls_lower}">{exp_cls}</span>'

            struct_cols = (
                f'<td>{ss_badge}</td>'
                f'<td>{rsa_html}</td>'
                f'<td style="font-size:0.82em">{exp_badge}</td>'
                f'<td style="font-size:0.80em;color:#666">{exp_note}</td>'
            )

        rows.append(
            f"<tr>"
            f"<td class='pos-cell'>{f['pos1']}</td>"
            f"<td>{res_span}</td>"
            f"<td>{f['category']}</td>"
            f"<td>{risk_badge}</td>"
            f"<td style='font-size:0.88em'>{f['label']}</td>"
            f"<td style='font-size:0.82em;color:#555'>{note}</td>"
            f"{struct_cols}"
            f"</tr>"
        )
    return "\n".join(rows)


def summary_cards_html(summary: dict) -> str:
    order = ["Deamidation", "Oxidation", "Isomerization", "Glycosylation", "Truncation"]
    cards = []
    for cat in order:
        items  = summary["by_category"].get(cat, [])
        high   = sum(1 for f in items if f["risk"] == "high")
        medium = sum(1 for f in items if f["risk"] == "medium")
        low    = sum(1 for f in items if f["risk"] == "low")
        info   = sum(1 for f in items if f["risk"] == "info")
        total  = len(items)
        bars   = ""
        if high:   bars += f'<span class="risk-badge risk-high">{high} High</span> '
        if medium: bars += f'<span class="risk-badge risk-medium">{medium} Med</span> '
        if low:    bars += f'<span class="risk-badge risk-low">{low} Low</span> '
        if info:   bars += f'<span class="risk-badge risk-info">{info} Info</span>'
        icon = {"Deamidation": "🧪", "Oxidation": "⚡", "Isomerization": "🔄",
                "Glycosylation": "🍬", "Truncation": "✂️"}.get(cat, "•")
        cards.append(
            f'<div class="summary-card">'
            f'<div class="card-icon">{icon}</div>'
            f'<div class="card-title">{cat}</div>'
            f'<div class="card-count">{total}</div>'
            f'<div class="card-bars">{bars}</div>'
            f'</div>'
        )
    return "\n".join(cards)


def hos_summary_card_html(stats: dict, source: str) -> str:
    """Render the HOS summary card showing SS composition and RSA overview."""
    n = stats["n"]
    counts = stats["counts"]
    total_assigned = n

    bar_segs = []
    for ss_char, label, color in [
        ("H", "Helix",  "#00b4d8"),
        ("E", "Strand", "#f9c74f"),
        ("T", "Turn",   "#90be6d"),
        ("C", "Coil",   "#dee2e6"),
    ]:
        cnt = counts.get(ss_char, 0)
        if cnt == 0:
            continue
        pct = cnt / total_assigned * 100 if total_assigned else 0
        txt_col = "#fff" if ss_char != "E" and ss_char != "C" else "#333"
        bar_segs.append(
            f'<span class="hos-bar-seg" '
            f'style="background:{color};color:{txt_col};width:{pct:.0f}%" '
            f'title="{label}: {cnt} res ({pct:.1f}%)">'
            f'{ss_char} {pct:.0f}%</span>'
        )

    exposed_pct  = stats["exposed"]  / n * 100 if n else 0
    buried_pct   = stats["buried"]   / n * 100 if n else 0

    return f"""
    <div class="section">
      <div class="section-header">🏗️ Higher-Order Structure Analysis</div>
      <div class="section-body">
        <div class="info-box">
          <strong>Source:</strong> {source}<br>
          Prediction is based on local sequence context only.
          For higher accuracy, provide a PDB structure file using <code>--pdb</code>.
        </div>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:16px;margin-bottom:20px">
          <div class="summary-card">
            <div class="card-icon">🌀</div>
            <div class="card-title">α-Helix</div>
            <div class="card-count">{counts.get("H",0)}</div>
            <div style="font-size:0.82em;color:#888">{counts.get("H",0)/n*100:.1f}% of residues</div>
          </div>
          <div class="summary-card">
            <div class="card-icon">🏹</div>
            <div class="card-title">β-Strand</div>
            <div class="card-count">{counts.get("E",0)}</div>
            <div style="font-size:0.82em;color:#888">{counts.get("E",0)/n*100:.1f}% of residues</div>
          </div>
          <div class="summary-card">
            <div class="card-icon">🔴</div>
            <div class="card-title">Surface-Exposed</div>
            <div class="card-count">{stats["exposed"]}</div>
            <div style="font-size:0.82em;color:#888">RSA &gt; 0.25 ({exposed_pct:.1f}%)</div>
          </div>
          <div class="summary-card">
            <div class="card-icon">🔵</div>
            <div class="card-title">Buried</div>
            <div class="card-count">{stats["buried"]}</div>
            <div style="font-size:0.82em;color:#888">RSA ≤ 0.10 ({buried_pct:.1f}%)</div>
          </div>
        </div>

        <div style="margin-bottom:14px">
          <div style="font-weight:600;margin-bottom:6px;color:#1a1a2e">
            Secondary Structure Composition
          </div>
          <div class="hos-bar-wrap" style="width:100%">
            {''.join(bar_segs)}
          </div>
        </div>

        <div style="margin-top:14px;display:flex;gap:16px;flex-wrap:wrap">
          <span><span style="display:inline-block;width:14px;height:14px;background:#00b4d8;border-radius:2px;vertical-align:middle;margin-right:4px"></span>α-Helix (H)</span>
          <span><span style="display:inline-block;width:14px;height:14px;background:#f9c74f;border-radius:2px;vertical-align:middle;margin-right:4px"></span>β-Strand (E)</span>
          <span><span style="display:inline-block;width:14px;height:14px;background:#90be6d;border-radius:2px;vertical-align:middle;margin-right:4px"></span>β-Turn (T)</span>
          <span><span style="display:inline-block;width:14px;height:14px;background:#dee2e6;border-radius:2px;vertical-align:middle;margin-right:4px"></span>Coil (C)</span>
          <span style="margin-left:20px;color:#888;font-size:0.9em">
            RSA track: <span style="color:#4682b4">■</span> Buried →
            <span style="color:#999">■</span> Partial →
            <span style="color:#c0392b">■</span> Exposed
          </span>
        </div>
      </div>
    </div>
    """


def build_html_report(sequences: list, title: str = "Protein Sequence Liability Analysis",
                       generated_by: str = "Protein Liability Analyzer v2") -> str:
    """
    Build the full HTML report.
    sequences: list of dicts with keys:
      name, seq, findings, summary, ss (optional), rsa (optional),
      hos_stats (optional), hos_source (optional)
    """
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    seq_sections = []

    for entry in sequences:
        name       = entry["name"]
        seq        = entry["seq"]
        findings   = entry["findings"]
        summary    = entry["summary"]
        ss         = entry.get("ss")
        rsa        = entry.get("rsa")
        has_struct = ss is not None and rsa is not None

        annotated  = build_annotated_sequence(seq, findings, ss=ss, rsa=rsa)
        table_rows = findings_table_html(findings, has_structure=has_struct)
        cards      = summary_cards_html(summary)
        total      = summary["total"]
        high_ct    = len(summary["by_risk"].get("high", []))

        alert_color = "#c0392b" if high_ct > 0 else "#2e86c1"
        alert_bg    = "#fde8e8" if high_ct > 0 else "#d6eaf8"
        alert_text  = (
            f"⚠️  {high_ct} high-risk liabilit{'y' if high_ct == 1 else 'ies'} detected."
            if high_ct > 0
            else "✅  No high-risk liabilities detected."
        )

        # Structure-specific columns header
        struct_th = ""
        if has_struct:
            struct_th = (
                "<th>Sec. Structure</th>"
                "<th>RSA</th>"
                "<th>Exposure</th>"
                "<th>Structural Note</th>"
            )

        # HOS section for this sequence
        hos_html = ""
        if has_struct and "hos_stats" in entry:
            hos_html = hos_summary_card_html(entry["hos_stats"], entry.get("hos_source", "predicted"))

        seq_sections.append(f"""
        {hos_html}
        <div class="section">
          <div class="section-header">
            {'📌 ' if len(sequences) > 1 else ''}{name}
            <span style="font-weight:400;font-size:0.88em;color:#888;margin-left:10px">
              {summary['seq_len']} residues · {total} liabilities flagged
            </span>
          </div>
          <div class="section-body">
            <div style="background:{alert_bg};color:{alert_color};padding:10px 16px;
                        border-radius:6px;margin-bottom:18px;font-weight:600;font-size:0.9em">
              {alert_text}
            </div>

            <div style="margin-bottom:18px">
              <div class="summary-grid">{cards}</div>
            </div>

            <div style="margin-bottom:24px">
              <div style="font-weight:600;margin-bottom:10px;color:#1a1a2e">
                Annotated Sequence
                <span style="font-weight:400;font-size:0.8em;color:#888;margin-left:8px">
                  (hover residues for details)
                </span>
              </div>
              {('<div style="font-size:0.8em;color:#666;margin-bottom:8px">Row 1: sequence with liability highlights &nbsp;|&nbsp; Row 2: secondary structure &nbsp;|&nbsp; Row 3: RSA heatmap</div>') if has_struct else ''}
              <div class="seq-block">{annotated}</div>
            </div>

            <div>
              <div style="font-weight:600;margin-bottom:10px;color:#1a1a2e">
                Liabilities Detail
              </div>
              <div style="overflow-x:auto">
                <table>
                  <thead>
                    <tr>
                      <th>Pos.</th>
                      <th>Residue(s)</th>
                      <th>Category</th>
                      <th>Risk</th>
                      <th>Type</th>
                      <th>Note</th>
                      {struct_th}
                    </tr>
                  </thead>
                  <tbody>
                    {table_rows}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        </div>
        """)

    legend_rows = liability_legend_html()

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
  <style>{CSS}</style>
</head>
<body>
<div class="page-wrap">

  <div class="report-header">
    <h1>🔬 {title}</h1>
    <div class="subtitle">
      Post-translational modification risk &amp; sequence liability assessment
      {'— with Higher-Order Structure Analysis' if any(e.get('ss') is not None for e in sequences) else ''}
    </div>
    <div class="meta">
      <span>📅 Generated: {now}</span>
      <span>🔢 Chains analysed: {len(sequences)}</span>
      <span>⚗️  {generated_by}</span>
      <span>🔒 All computation local — sequences not transmitted externally</span>
    </div>
  </div>

  {''.join(seq_sections)}

  <div class="section">
    <div class="section-header">📖 Liability Legend &amp; References</div>
    <div class="section-body">
      <div style="overflow-x:auto">
        <table class="legend-table">
          <thead>
            <tr><th>Code</th><th>Liability</th><th>Risk</th><th>Description</th></tr>
          </thead>
          <tbody>{legend_rows}</tbody>
        </table>
      </div>
      <div style="margin-top:20px;font-size:0.82em;color:#888;line-height:1.6">
        <strong>Liability references:</strong>
        Geiger &amp; Clarke (1987) <em>J Biol Chem</em>;
        Robinson &amp; Robinson (2001) <em>PNAS</em>;
        Lam et al. (1997) <em>J Pharm Sci</em>;
        Clarke (1987) <em>J Biol Chem</em>;
        Piszkiewicz et al. (1970) <em>BBRC</em>;
        Dick et al. (2008) <em>Biotechnol Bioeng</em>;
        Harris (1995) <em>J Chromatogr</em>;
        Apweiler et al. (1999) <em>Biochim Biophys Acta</em>.
        <br><br>
        <strong>Structure prediction references:</strong>
        Chou &amp; Fasman (1974) <em>Biochemistry</em> 13(2):211-222 (secondary structure);
        Rost &amp; Sander (1994) <em>Proteins</em> 20(3):216-226 (RSA statistics);
        Tien et al. (2013) <em>PLOS ONE</em> 8(11):e80635 (maximum ASA values);
        Shrake &amp; Rupley (1973) <em>J Mol Biol</em> 79(2):351-371 (SASA algorithm).
      </div>
    </div>
  </div>

  <div class="footer">
    Generated by Protein Liability Analyzer v2 · Letarte Scientific Consulting ·
    For research and development use only
  </div>

</div>
</body>
</html>"""


# ══════════════════════════════════════════════════════════════════════════════
# 12. CONSOLE OUTPUT
# ══════════════════════════════════════════════════════════════════════════════

ANSI = {
    "reset":  "\033[0m",
    "bold":   "\033[1m",
    "red":    "\033[91m",
    "orange": "\033[93m",
    "green":  "\033[92m",
    "blue":   "\033[94m",
    "purple": "\033[95m",
    "cyan":   "\033[96m",
    "grey":   "\033[90m",
    "teal":   "\033[36m",
    "gold":   "\033[33m",
}

RISK_ANSI = {
    "high":   ANSI["red"],
    "medium": ANSI["orange"],
    "low":    ANSI["green"],
    "info":   ANSI["blue"],
}

SS_ANSI = {
    "H": ANSI["cyan"],
    "E": ANSI["gold"],
    "T": ANSI["green"],
    "C": ANSI["grey"],
}


def print_console_report(name: str, seq: str, findings: list, summary: dict,
                          ss: Optional[list] = None, rsa: Optional[list] = None,
                          use_color: bool = True) -> None:
    C = ANSI  if use_color else {k: "" for k in ANSI}
    R = RISK_ANSI if use_color else {k: "" for k in RISK_ANSI}
    S = SS_ANSI   if use_color else {k: "" for k in SS_ANSI}

    sep = "─" * 72
    print(f"\n{C['bold']}{sep}{C['reset']}")
    print(f"{C['bold']}  {name}{C['reset']}   ({summary['seq_len']} residues)")
    print(f"{sep}")

    # Liability summary
    order = ["Deamidation", "Oxidation", "Isomerization", "Glycosylation", "Truncation"]
    for cat in order:
        items = summary["by_category"].get(cat, [])
        if not items:
            continue
        counts = defaultdict(int)
        for f in items:
            counts[f["risk"]] += 1
        detail = ", ".join(
            f'{R[r]}{n} {r}{C["reset"]}'
            for r, n in sorted(counts.items(), key=lambda x: RISK_ORDER[x[0]])
        )
        print(f"  {C['bold']}{cat:<18}{C['reset']}  {len(items):>3} findings  [{detail}]")

    high_ct = len(summary["by_risk"].get("high", []))
    if high_ct:
        print(f"\n  {C['red']}{C['bold']}⚠  {high_ct} high-risk liabilities found!{C['reset']}")
    else:
        print(f"\n  {C['green']}✓  No high-risk liabilities detected.{C['reset']}")

    # HOS summary
    if ss is not None and rsa is not None:
        h_count = ss.count("H")
        e_count = ss.count("E")
        t_count = ss.count("T")
        c_count = ss.count("C")
        n       = len(ss)
        avg_rsa = sum(rsa) / len(rsa) if rsa else 0
        exposed = sum(1 for r in rsa if r > 0.25)
        buried  = sum(1 for r in rsa if r <= 0.10)

        print(f"\n  {C['bold']}Higher-Order Structure:{C['reset']}")
        print(
            f"  {S['H']}Helix {h_count:>4} ({h_count/n*100:4.1f}%){C['reset']}  "
            f"{S['E']}Strand {e_count:>3} ({e_count/n*100:4.1f}%){C['reset']}  "
            f"{S['T']}Turn {t_count:>4} ({t_count/n*100:4.1f}%){C['reset']}  "
            f"{S['C']}Coil {c_count:>4} ({c_count/n*100:4.1f}%){C['reset']}"
        )
        print(
            f"  Avg RSA: {avg_rsa:.3f}  |  "
            f"{C['red']}Exposed (RSA>0.25): {exposed}{C['reset']}  |  "
            f"{C['blue']}Buried (RSA≤0.10): {buried}{C['reset']}"
        )

    # Detailed findings
    print(f"\n  {C['bold']}Detailed Findings:{C['reset']}")
    sorted_f = sorted(findings, key=lambda f: (RISK_ORDER.get(f["risk"], 9), f["pos0"]))
    for f in sorted_f:
        risk_str = f"{R[f['risk']]}{f['risk'].upper():6}{C['reset']}"
        note     = f"  ← {f['note']}" if f.get("note") else ""

        struct_str = ""
        if ss is not None and rsa is not None:
            ss_char  = f.get("ss", "?")
            rsa_val  = f.get("rsa")
            exp_cls  = f.get("exposure_class", "")
            ss_color = S.get(ss_char, "")
            rsa_str  = f"{rsa_val:.2f}" if rsa_val is not None else "N/A"
            struct_str = (
                f"  [{ss_color}{ss_char}{C['reset']}]"
                f"  RSA={rsa_str}"
                f"  [{exp_cls}]"
            )

        print(
            f"  Pos {f['pos1']:>5}  {C['bold']}{f['residues']:<5}{C['reset']}"
            f"  {risk_str}  {f['category']:<16}  {f['label']}{note}{struct_str}"
        )

    print(f"{sep}\n")


# ══════════════════════════════════════════════════════════════════════════════
# 13. MAIN PROCESSING
# ══════════════════════════════════════════════════════════════════════════════

def run_structure_analysis_sequence(seq: str, verbose: bool = True):
    """
    Predict secondary structure and RSA from sequence alone.
    Returns (ss_list, rsa_list, source_string).
    """
    if verbose:
        print("  ▸ Predicting secondary structure (Chou-Fasman) …", end=" ", flush=True)
    ss = predict_secondary_structure(seq)
    if verbose:
        h_pct = ss.count("H") / len(ss) * 100
        e_pct = ss.count("E") / len(ss) * 100
        print(f"done  (H:{h_pct:.0f}%  E:{e_pct:.0f}%)")

    if verbose:
        print("  ▸ Estimating relative solvent accessibility …", end=" ", flush=True)
    rsa = predict_rsa_from_sequence(seq, ss)
    if verbose:
        avg = sum(rsa) / len(rsa)
        print(f"done  (avg RSA: {avg:.3f})")

    return ss, rsa, "Chou-Fasman (sequence-based prediction)"


def run_structure_analysis_pdb(pdb_text: str, seq: str,
                                chain_id: Optional[str] = None,
                                verbose: bool = True):
    """
    Analyse structure from a PDB file.
    Returns (ss_list, rsa_list, source_string).
    """
    # Parse atoms
    if verbose:
        print(f"  ▸ Parsing PDB atoms …", end=" ", flush=True)
    atoms = parse_pdb_atoms(pdb_text, chain_id=chain_id)
    if not atoms:
        print("  [Warning] No ATOM records found in PDB. Falling back to sequence prediction.")
        return run_structure_analysis_sequence(seq, verbose)
    if verbose:
        print(f"done  ({len(atoms)} heavy atoms)")

    # Map PDB residues to input sequence
    pdb_residues = extract_pdb_residues(atoms)
    seq_to_pdb, match_pct = map_pdb_to_sequence(pdb_residues, seq, chain_id)

    if verbose:
        print(f"  ▸ Sequence alignment: {match_pct:.1f}% identity to PDB chain")
    if match_pct < 50:
        print(
            f"  [Warning] Low sequence match ({match_pct:.1f}%).  "
            "Check --chain argument.  Falling back to sequence-based prediction."
        )
        return run_structure_analysis_sequence(seq, verbose)

    # Secondary structure from HELIX/SHEET records
    if verbose:
        print("  ▸ Assigning secondary structure from HELIX/SHEET records …", end=" ", flush=True)
    pdb_ss_map = parse_pdb_helix_sheet(pdb_text)
    # Build per-residue SS arrays aligned to input sequence
    n = len(seq)
    ss = ["C"] * n
    for seq_pos, (chain, res_num) in seq_to_pdb.items():
        if seq_pos < n:
            ss[seq_pos] = pdb_ss_map.get((chain, res_num), "C")
    h_pct = ss.count("H") / n * 100
    e_pct = ss.count("E") / n * 100
    if verbose:
        print(f"done  (H:{h_pct:.0f}%  E:{e_pct:.0f}%)")

    # Shrake-Rupley SASA
    if verbose:
        print("  ▸ Computing SASA (Shrake-Rupley) …")
    res_sasa = compute_sasa_per_residue(atoms, progress=verbose)
    if verbose:
        print()

    # Build RSA array aligned to input sequence
    rsa = []
    for i, aa in enumerate(seq):
        pdb_info = seq_to_pdb.get(i)
        if pdb_info:
            _, res_num = pdb_info
            sasa = res_sasa.get(res_num, None)
            max_asa = MAX_ASA.get(aa, 200.0)
            rsa_i = min(1.0, sasa / max_asa) if sasa is not None else 0.40
        else:
            rsa_i = COIL_AVG_RSA.get(aa, 0.40)   # fallback for unmapped residues
        rsa.append(round(max(0.0, min(1.0, rsa_i)), 3))

    pdb_name = pdb_residues[0][0] if pdb_residues else "?"
    source = (
        f"PDB structure — HELIX/SHEET records (secondary structure) + "
        f"Shrake-Rupley SASA, probe=1.4 Å (solvent exposure)"
    )
    return ss, rsa, source


def process_sequences(pairs: list, pdb_text: Optional[str] = None,
                       chain_id: Optional[str] = None,
                       run_hos: bool = False,
                       verbose: bool = True) -> list:
    """Process (name, raw_seq) pairs and return list of result dicts."""
    results = []
    for name, raw_seq in pairs:
        seq = clean_sequence(raw_seq)
        _, warns = validate_sequence(seq)
        for w in warns:
            print(f"  [Warning] {name}: {w}", file=sys.stderr)

        if verbose:
            print(f"\n  Analysing: {name}  ({len(seq)} residues)")

        findings = find_liabilities(seq)
        summary  = summarize(findings, len(seq))

        entry = {"name": name, "seq": seq, "findings": findings, "summary": summary}

        # Run structure analysis if requested
        if run_hos or pdb_text is not None:
            if pdb_text is not None:
                ss, rsa, source = run_structure_analysis_pdb(
                    pdb_text, seq, chain_id=chain_id, verbose=verbose
                )
            else:
                ss, rsa, source = run_structure_analysis_sequence(seq, verbose=verbose)

            add_structural_context(findings, ss, rsa, source=source)
            entry["ss"]        = ss
            entry["rsa"]       = rsa
            entry["hos_stats"] = hos_stats(ss, rsa)
            entry["hos_source"] = source

        results.append(entry)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# 14. MAIN ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Protein Sequence Liability Analyzer v2 — with HOS & Solvent Exposure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # ── Input ──────────────────────────────────────────────────────────────────
    group = parser.add_mutually_exclusive_group()
    group.add_argument("fasta_file", nargs="?",
                       help="FASTA file with one or more protein sequences")
    group.add_argument("--sequence", "-s", metavar="SEQ",
                       help="Raw amino-acid sequence (single-letter code)")
    group.add_argument("--interactive", "-i", action="store_true",
                       help="Prompt for sequence(s) interactively")

    # ── Structure ──────────────────────────────────────────────────────────────
    parser.add_argument(
        "--hos", action="store_true",
        help=(
            "Run sequence-based Higher Order Structure analysis "
            "(Chou-Fasman secondary structure + statistical RSA prediction). "
            "Enabled automatically when --pdb is provided."
        ),
    )
    parser.add_argument(
        "--pdb", metavar="FILE",
        help=(
            "PDB structure file for accurate secondary-structure assignment "
            "(HELIX/SHEET records) and Shrake-Rupley SASA calculation."
        ),
    )
    parser.add_argument(
        "--chain", metavar="ID", default=None,
        help="PDB chain ID to use (default: all chains, best match selected).",
    )

    # ── Output ──────────────────────────────────────────────────────────────────
    parser.add_argument("--output", "-o", metavar="FILE", default=None,
                        help="Output HTML report path")
    parser.add_argument("--no-color", action="store_true",
                        help="Disable ANSI colours in console output")
    parser.add_argument("--no-html",  action="store_true",
                        help="Skip HTML report generation")

    args = parser.parse_args()

    # ── Load input sequences ───────────────────────────────────────────────────
    pairs: list = []
    out_stem = "protein"

    if args.fasta_file:
        path = Path(args.fasta_file)
        if not path.exists():
            print(f"Error: file not found: {path}", file=sys.stderr)
            sys.exit(1)
        text = path.read_text(encoding="utf-8", errors="replace")
        pairs = parse_fasta(text) if text.lstrip().startswith(">") else [(path.stem, text)]
        out_stem = path.stem

    elif args.sequence:
        pairs = [("Input Sequence", args.sequence)]

    elif args.interactive:
        print("═" * 60)
        print("  Protein Liability Analyzer v2 — Interactive Mode")
        print("═" * 60)
        print("Paste one or more sequences (FASTA format or raw).")
        print("Enter a blank line followed by END to finish.\n")
        lines = []
        while True:
            try:
                line = input()
            except EOFError:
                break
            if line.strip().upper() == "END":
                break
            lines.append(line)
        text  = "\n".join(lines)
        pairs = parse_fasta(text) if ">" in text else [("Input Sequence", text)]

    else:
        if sys.stdin.isatty():
            print("Usage: protein_liability_analyzer_v2.py [fasta_file] [--sequence SEQ]")
            print("       protein_liability_analyzer_v2.py --interactive")
            print("       protein_liability_analyzer_v2.py --help")
            sys.exit(0)
        text  = sys.stdin.read()
        pairs = parse_fasta(text) if ">" in text else [("stdin", text)]

    if not pairs:
        print("No sequences found.", file=sys.stderr)
        sys.exit(1)

    # ── Load PDB if provided ───────────────────────────────────────────────────
    pdb_text = None
    if args.pdb:
        pdb_path = Path(args.pdb)
        if not pdb_path.exists():
            print(f"Error: PDB file not found: {pdb_path}", file=sys.stderr)
            sys.exit(1)
        pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
        print(f"📂  PDB loaded: {pdb_path}  ({len(pdb_text.splitlines())} lines)")

    # ── Decide whether to run HOS ──────────────────────────────────────────────
    run_hos = args.hos or (pdb_text is not None)

    # Interactive: ask about HOS if not already enabled and no PDB
    if args.interactive and not run_hos:
        print(
            "\nNo PDB structure provided.\n"
            "Would you like to predict higher-order structure (Chou-Fasman\n"
            "secondary structure + sequence-based RSA estimation)?\n"
            "All computation is local — no internet connection required."
        )
        choice = input("Run HOS analysis? [Y/n]: ").strip().lower()
        if choice in ("", "y", "yes"):
            run_hos = True

    # ── Process ────────────────────────────────────────────────────────────────
    use_color = not args.no_color and sys.stdout.isatty()
    print()

    results = process_sequences(
        pairs,
        pdb_text=pdb_text,
        chain_id=args.chain,
        run_hos=run_hos,
        verbose=True,
    )

    # ── Console output ─────────────────────────────────────────────────────────
    for entry in results:
        print_console_report(
            entry["name"], entry["seq"], entry["findings"], entry["summary"],
            ss=entry.get("ss"), rsa=entry.get("rsa"),
            use_color=use_color,
        )

    # ── HTML report ────────────────────────────────────────────────────────────
    if not args.no_html:
        if args.output:
            html_path = Path(args.output)
        else:
            suffix    = "_hos" if run_hos else ""
            html_path = Path(f"{out_stem}_liabilities{suffix}.html")

        title = "Protein Sequence Liability Analysis"
        if run_hos:
            if pdb_text:
                title += " + PDB Structure"
            else:
                title += " + HOS Prediction"

        html = build_html_report(results, title=title)
        html_path.write_text(html, encoding="utf-8")
        print(f"📄  HTML report saved → {html_path}")

    if run_hos and pdb_text is None:
        print(
            "\n  ℹ️   Note: HOS results are sequence-based PREDICTIONS.\n"
            "      Accuracy can be improved by providing a crystal structure\n"
            "      or homology model via --pdb.\n"
        )


if __name__ == "__main__":
    main()
