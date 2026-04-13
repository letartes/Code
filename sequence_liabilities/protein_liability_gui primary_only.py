#!/usr/bin/env python3
"""
Protein Sequence Liability Analyzer — GUI
==========================================
Graphical interface for identifying PTM risks and sequence liabilities
in biopharmaceutical protein sequences.

Detected liabilities:
  Deamidation · Oxidation · Isomerization
  Glycosylation · Pyroglutamate · Truncation/Clipping

CDR annotation (VH / VL) is performed automatically when antibody
variable-domain sequences are detected.

Usage:
  python protein_liability_gui.py
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import webbrowser
import tempfile
import os
import re
import sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime


# ═══════════════════════════════════════════════════════════════════════════════
# ❶  ANALYSIS ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

LIABILITIES = {
    # ── Deamidation ──────────────────────────────────────────────────────────
    "deamid_high": {
        "label":         "Deamidation – High Risk",
        "display_label": "Deamidation (NG / NS)",
        "short":         "Deam-H",
        "category":      "Deamidation",
        "risk":          "high",
        "color_text":    "#c0392b",
        "color_bg":      "#fde8e8",
        "color_border":  "#e74c3c",
        "description": (
            "NG and NS motifs. Asn–Gly is the most susceptible sequence to "
            "deamidation, forming succinimide intermediates rapidly at "
            "physiological pH. Asn–Ser is second most prone."
        ),
        "span_length": 2,
    },
    "deamid_medium": {
        "label":         "Deamidation – Medium Risk",
        "display_label": "Deamidation (NA / NT / NQ / NR / QG / QS)",
        "short":         "Deam-M",
        "category":      "Deamidation",
        "risk":          "medium",
        "color_text":    "#d35400",
        "color_bg":      "#fef0e6",
        "color_border":  "#e67e22",
        "description": (
            "NA, NT, NQ, NR, NH motifs and Gln-based equivalents (QG, QS). "
            "Glutamine can deamidate, though more slowly than asparagine."
        ),
        "span_length": 2,
    },

    # ── Oxidation ─────────────────────────────────────────────────────────────
    "oxid_high": {
        "label":         "Oxidation – High Risk",
        "display_label": "Oxidation (Met, Trp)",
        "short":         "Ox-H",
        "category":      "Oxidation",
        "risk":          "high",
        "color_text":    "#6c3483",
        "color_bg":      "#f4ecf7",
        "color_border":  "#9b59b6",
        "description": (
            "Met (M) and Trp (W). Methionine sulfoxide formation is a primary "
            "oxidative liability; tryptophan oxidation is common under "
            "photo-oxidative stress."
        ),
        "span_length": 1,
    },
    "oxid_medium": {
        "label":         "Oxidation – Medium Risk",
        "display_label": "Oxidation (His, Cys)",
        "short":         "Ox-M",
        "category":      "Oxidation",
        "risk":          "medium",
        "color_text":    "#1a5276",
        "color_bg":      "#eaf4fb",
        "color_border":  "#2980b9",
        "description": (
            "His (H) and Cys (C). Histidine oxidised to 2-oxo-histidine under "
            "metal-catalysed conditions; free cysteines prone to disulfide "
            "scrambling and sulfonic acid formation."
        ),
        "span_length": 1,
    },
    "oxid_low": {
        "label":         "Oxidation – Low Risk",
        "display_label": "Oxidation (Phe, Tyr)",
        "short":         "Ox-L",
        "category":      "Oxidation",
        "risk":          "low",
        "color_text":    "#0e6655",
        "color_bg":      "#e8f8f5",
        "color_border":  "#1abc9c",
        "description": (
            "Phe (F) and Tyr (Y). Can oxidise under extreme photo- or "
            "metal-catalysed conditions; generally lower concern."
        ),
        "span_length": 1,
    },

    # ── Isomerization ─────────────────────────────────────────────────────────
    "isom_high": {
        "label":         "Isomerization – High Risk",
        "display_label": "Isomerization (DG / DS)",
        "short":         "Isom-H",
        "category":      "Isomerization",
        "risk":          "high",
        "color_text":    "#1e8449",
        "color_bg":      "#eafaf1",
        "color_border":  "#27ae60",
        "description": (
            "DG and DS motifs. Asp–Gly undergoes succinimide-mediated "
            "isomerization to isoAsp most readily; Asp–Ser also highly prone."
        ),
        "span_length": 2,
    },
    "isom_medium": {
        "label":         "Isomerization – Medium Risk",
        "display_label": "Isomerization (DT / DA / DN / DD / DH)",
        "short":         "Isom-M",
        "category":      "Isomerization",
        "risk":          "medium",
        "color_text":    "#196f3d",
        "color_bg":      "#d5f5e3",
        "color_border":  "#239b56",
        "description": (
            "DT, DA, DN, DD, DH motifs — slower isomerization relative to DG/DS."
        ),
        "span_length": 2,
    },

    # ── Glycosylation ─────────────────────────────────────────────────────────
    "n_glycan": {
        "label":         "N-Glycosylation Sequon (NXS/T, X≠P)",
        "display_label": "N-Glycosylation Sequon (NXS/T)",
        "short":         "N-Glycan",
        "category":      "Glycosylation",
        "risk":          "info",
        "color_text":    "#154360",
        "color_bg":      "#d6eaf8",
        "color_border":  "#2e86c1",
        "description": (
            "Canonical N-linked glycosylation sequon N–X–S/T (X ≠ Pro). "
            "May represent intended or adventitious glycosylation."
        ),
        "span_length": 3,
    },
    "o_glycan": {
        "label":         "O-Glycosylation Site (S/T-P motif)",
        "display_label": "O-Glycosylation (S/T–P motif)",
        "short":         "O-Glycan",
        "category":      "Glycosylation",
        "risk":          "info",
        "color_text":    "#4a235a",
        "color_bg":      "#f5eef8",
        "color_border":  "#8e44ad",
        "description": (
            "Ser–Pro and Thr–Pro motifs are characteristic O-glycosylation "
            "sites (mucin-type)."
        ),
        "span_length": 2,
    },

    # ── Pyroglutamate ─────────────────────────────────────────────────────────
    "pyroglu": {
        "label":         "Pyroglutamate Formation",
        "display_label": "Pyroglutamate (N-term Q/E cyclisation)",
        "short":         "pyroGlu",
        "category":      "Pyroglutamate",
        "risk":          "medium",
        "color_text":    "#7e5109",
        "color_bg":      "#fef5e7",
        "color_border":  "#f0a500",
        "description": (
            "N-terminal Gln (Q) spontaneously cyclises to pyroglutamate "
            "(–17 Da, loss of NH₃). N-terminal Glu (E) can also cyclise "
            "(–18 Da). Both cause charge heterogeneity and are CQAs under "
            "ICH Q6B. Distinct from truncation – the protein chain is intact."
        ),
        "span_length": 1,
    },

    # ── Truncation / Clipping ─────────────────────────────────────────────────
    "trunc_high": {
        "label":         "Truncation – High Risk",
        "display_label": "Asp–Pro Bond Cleavage",
        "short":         "Trunc-H",
        "category":      "Truncation",
        "risk":          "high",
        "color_text":    "#78281f",
        "color_bg":      "#fdedec",
        "color_border":  "#c0392b",
        "description": (
            "Asp–Pro (DP) bonds are acid-labile and cleave readily under "
            "mildly acidic purification or formulation conditions."
        ),
        "span_length": 2,
    },
    "trunc_medium": {
        "label":         "Truncation – Medium Risk",
        "display_label": "C-term Lys Clipping / N-term Met Removal",
        "short":         "Trunc-M",
        "category":      "Truncation",
        "risk":          "medium",
        "color_text":    "#784212",
        "color_bg":      "#fef9e7",
        "color_border":  "#f39c12",
        "description": (
            "C-terminal Lys clipping (common in CHO-expressed IgGs) and "
            "N-terminal Met removal by Met aminopeptidase."
        ),
        "span_length": 1,
    },
}

RISK_ORDER = {"high": 0, "medium": 1, "low": 2, "info": 3}
VALID_AA   = set("ACDEFGHIKLMNPQRSTVWY")


# ── Sequence parsing ──────────────────────────────────────────────────────────

def parse_fasta(text: str) -> list:
    sequences, current_name, current_seq = [], "Sequence", []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_seq:
                sequences.append((current_name, "".join(current_seq).upper()))
            current_name = line[1:].strip() or "Sequence"
            current_seq  = []
        else:
            current_seq.append(re.sub(r"[^A-Za-z]", "", line))
    if current_seq:
        sequences.append((current_name, "".join(current_seq).upper()))
    return sequences


def clean_sequence(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper()


# ── CDR detection ─────────────────────────────────────────────────────────────

def find_cdr_regions(seq: str) -> list:
    """
    Detect CDR regions in antibody variable-domain sequences using conserved
    anchor residues (approximate Chothia-style boundaries).

    Returns a sorted list of (start_0, end_0_inclusive, name) tuples, or []
    when the sequence does not contain a recognisable antibody V-domain.

    NOTE: This is a structural approximation. For precise Kabat/Chothia/IMGT
    numbering use ANARCI or similar tools.
    """
    s = seq.upper()
    cdrs = []

    # ── VH ───────────────────────────────────────────────────────────────────
    # FR4 anchor: WGxGT/S  (highly conserved in all VH domains)
    m_fr4 = re.search(r'WG[QKHR]G[TS]', s)
    if m_fr4:
        fr4 = m_fr4.start()

        # CDR-H3: from (last C within 30 aa of FR4) + 1  →  FR4 - 1
        win3   = s[max(0, fr4 - 30): fr4]
        rc3    = win3.rfind('C')
        h3_found = None
        if rc3 >= 0:
            h3_s = max(0, fr4 - 30) + rc3 + 1
            h3_e = fr4 - 1
            if 2 <= h3_e - h3_s + 1 <= 28:
                cdrs.append((h3_s, h3_e, 'CDR-H3'))
                h3_found = h3_s

        # FR2 anchor: W[VIL..][RQH]  (conserved VH tryptophan)
        fr2_m = None
        for m in re.finditer(r'W[VILMAF][RQHK]', s[:fr4 - 35]):
            fr2_m = m
        if fr2_m:
            fr2 = fr2_m.start()

            # CDR-H1: (last C within 22 aa of FR2) + 3  →  FR2 - 1
            win1  = s[max(0, fr2 - 22): fr2]
            rc1   = win1.rfind('C')
            if rc1 >= 0:
                c_abs = max(0, fr2 - 22) + rc1
                h1_s  = c_abs + 3          # skip 2 FR1 residues after conserved C
                h1_e  = fr2 - 1
                if 4 <= h1_e - h1_s + 1 <= 15:
                    cdrs.append((h1_s, h1_e, 'CDR-H1'))

            # CDR-H2: FR2 body ≈ 14 aa after WxR, then ~17 aa CDR
            h2_s = fr2 + 14
            h2_e = h2_s + 16
            if h3_found is not None:
                h2_e = min(h2_e, h3_found - 29)
            if 5 <= h2_e - h2_s + 1 <= 22:
                cdrs.append((h2_s, h2_e, 'CDR-H2'))

        return sorted(cdrs, key=lambda x: x[0])

    # ── VL kappa / lambda ────────────────────────────────────────────────────
    # FR4 anchor: FGxGT/S
    m_fr4_vl = re.search(r'FG[QNRST]G[TS]', s)
    if m_fr4_vl:
        fr4 = m_fr4_vl.start()

        # CDR-L3: (last C within 18 aa of FR4) + 1  →  FR4 - 1
        win3  = s[max(0, fr4 - 18): fr4]
        rc3   = win3.rfind('C')
        l3_found = None
        if rc3 >= 0:
            l3_s = max(0, fr4 - 18) + rc3 + 1
            l3_e = fr4 - 1
            if 3 <= l3_e - l3_s + 1 <= 14:
                cdrs.append((l3_s, l3_e, 'CDR-L3'))
                l3_found = l3_s

        # FR2 anchor: W[YFH][LQP]  (VL conserved tryptophan)
        fr2_m = None
        for m in re.finditer(r'W[YFH][LQP]', s[:fr4 - 35]):
            fr2_m = m
        if fr2_m:
            fr2 = fr2_m.start()

            # CDR-L1: (last C within 28 aa of FR2) + 1  →  FR2 - 1
            win1 = s[max(0, fr2 - 28): fr2]
            rc1  = win1.rfind('C')
            if rc1 >= 0:
                c_abs = max(0, fr2 - 28) + rc1
                l1_s  = c_abs + 1
                l1_e  = fr2 - 1
                if 6 <= l1_e - l1_s + 1 <= 18:
                    cdrs.append((l1_s, l1_e, 'CDR-L1'))

            # CDR-L2: FR2 body ≈ 13 aa, then ~7 aa CDR
            l2_s = fr2 + 13
            l2_e = l2_s + 6
            if l3_found is not None:
                l2_e = min(l2_e, l3_found - 29)
            if 5 <= l2_e - l2_s + 1 <= 12:
                cdrs.append((l2_s, l2_e, 'CDR-L2'))

        return sorted(cdrs, key=lambda x: x[0])

    return []


def cdr_name_at(pos0: int, cdrs: list) -> str | None:
    """Return the CDR name if pos0 falls within any CDR, else None."""
    for (start, end, name) in cdrs:
        if start <= pos0 <= end:
            return name
    return None


# ── Liability detection ───────────────────────────────────────────────────────

def find_liabilities(seq: str, cdrs: list | None = None) -> list:
    """
    Scan seq for all liability motifs. If cdrs is provided (from
    find_cdr_regions), each finding is annotated with 'cdr' (name or None).
    """
    if cdrs is None:
        cdrs = []
    findings = []
    n = len(seq)

    def add(key, pos0, span, note=""):
        f = {
            "type":     key,
            "pos0":     pos0,
            "span":     span,
            "pos1":     pos0 + 1,
            "residues": seq[pos0: pos0 + span],
            "note":     note,
            "cdr":      None,
            **{k: v for k, v in LIABILITIES[key].items()
               if k not in ("description", "span_length")},
        }
        # Annotate CDR membership (use start position of the motif)
        f["cdr"] = cdr_name_at(pos0, cdrs)
        findings.append(f)

    # Deamidation
    HIGH_DEAMID = {"NG", "NS"}
    MED_DEAMID  = {"NA", "NT", "NQ", "NR", "NH", "QG", "QS"}
    for i in range(n - 1):
        pair = seq[i:i+2]
        if pair in HIGH_DEAMID:
            add("deamid_high", i, 2)
        elif pair in MED_DEAMID:
            add("deamid_medium", i, 2)

    # Oxidation
    for i, aa in enumerate(seq):
        if aa in "MW":
            add("oxid_high", i, 1)
        elif aa in "HC":
            add("oxid_medium", i, 1)
        elif aa in "FY":
            add("oxid_low", i, 1)

    # Isomerization
    HIGH_ISOM = {"DG", "DS"}
    MED_ISOM  = {"DT", "DA", "DN", "DD", "DH"}
    for i in range(n - 1):
        pair = seq[i:i+2]
        if pair in HIGH_ISOM:
            add("isom_high", i, 2)
        elif pair in MED_ISOM:
            add("isom_medium", i, 2)

    # N-Glycosylation sequon NXS/T (X ≠ P)
    for i in range(n - 2):
        if seq[i] == "N" and seq[i+1] != "P" and seq[i+2] in "ST":
            add("n_glycan", i, 3)

    # O-Glycosylation (S/T–P)
    for i in range(n - 1):
        if seq[i] in "ST" and seq[i+1] == "P":
            add("o_glycan", i, 2)

    # Pyroglutamate — N-terminal Q or E (separate from truncation)
    if seq and seq[0] in "QE":
        note = (
            "N-term Gln → pyroglutamate (–17 Da, loss of NH₃)"
            if seq[0] == "Q"
            else "N-term Glu → pyroglutamate (–18 Da, loss of H₂O)"
        )
        add("pyroglu", 0, 1, note)

    # Truncation: Asp–Pro acid-labile bond
    for i in range(n - 1):
        if seq[i:i+2] == "DP":
            add("trunc_high", i, 2, "Acid-labile Asp–Pro bond")

    # Truncation: C-terminal Lys clipping (CHO-expressed IgGs)
    if n > 0 and seq[-1] == "K":
        add("trunc_medium", n - 1, 1, "C-terminal Lys clipping (CHO expression)")

    # Truncation: N-terminal Met removal
    if seq and seq[0] == "M":
        add("trunc_medium", 0, 1, "N-terminal Met — Met aminopeptidase removal")

    return findings


def summarize(findings: list, seq_len: int) -> dict:
    by_category = defaultdict(list)
    by_risk     = defaultdict(list)
    for f in findings:
        by_category[f["category"]].append(f)
        by_risk[f["risk"]].append(f)
    return {
        "total":       len(findings),
        "by_category": dict(by_category),
        "by_risk":     dict(by_risk),
        "seq_len":     seq_len,
    }


def process_sequences(pairs: list) -> list:
    """Return list of (name, seq, findings, summary, cdrs)."""
    results = []
    for name, raw in pairs:
        seq      = clean_sequence(raw)
        cdrs     = find_cdr_regions(seq)
        findings = find_liabilities(seq, cdrs)
        summary  = summarize(findings, len(seq))
        results.append((name, seq, findings, summary, cdrs))
    return results


# ── HTML report ───────────────────────────────────────────────────────────────

CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  background: #f0f2f5; color: #1a1a2e; font-size: 15px;
}
.page-wrap { max-width: 1200px; margin: 0 auto; padding: 24px 16px 60px; }
.report-header {
  background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
  color: #fff; border-radius: 12px; padding: 32px 36px; margin-bottom: 28px;
}
.report-header h1 { font-size: 1.9em; font-weight: 700; margin-bottom: 6px; }
.report-header .subtitle { color: #a8b2c8; font-size: 0.95em; }
.report-header .meta { margin-top: 16px; font-size: 0.85em; color: #c0c8dc;
  display:flex; gap:24px; flex-wrap:wrap; }
.summary-grid {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(155px, 1fr));
  gap: 14px; margin-bottom: 28px;
}
.summary-card {
  background: #fff; border-radius: 10px; padding: 20px 16px;
  box-shadow: 0 1px 6px rgba(0,0,0,.08); text-align: center;
}
.card-icon  { font-size: 1.8em; margin-bottom: 6px; }
.card-title { font-weight: 600; font-size: 0.9em; color: #555; margin-bottom: 4px; }
.card-count { font-size: 2em; font-weight: 700; color: #1a1a2e; margin-bottom: 8px; }
.card-bars  { display:flex; flex-wrap:wrap; gap:4px; justify-content:center; }
.risk-badge {
  display: inline-block; padding: 2px 7px; border-radius: 20px;
  font-size: 0.72em; font-weight: 700; letter-spacing: .4px; white-space: nowrap;
}
.risk-high   { background:#fde8e8; color:#c0392b; border:1px solid #e74c3c; }
.risk-medium { background:#fef0e6; color:#d35400; border:1px solid #e67e22; }
.risk-low    { background:#e8f8f5; color:#0e6655; border:1px solid #1abc9c; }
.risk-info   { background:#d6eaf8; color:#154360; border:1px solid #2e86c1; }
.section {
  background: #fff; border-radius: 10px;
  box-shadow: 0 1px 6px rgba(0,0,0,.08); margin-bottom: 24px; overflow: hidden;
}
.section-header {
  padding: 16px 24px; border-bottom: 1px solid #eee;
  font-weight: 700; font-size: 1.05em; color: #1a1a2e; background: #fafbfc;
}
.section-body { padding: 20px 24px; }
.seq-block {
  font-family: "Courier New", monospace; font-size: 0.92em; line-height: 2.2;
  letter-spacing: 0.05em; background: #f8f9fa; border-radius: 8px;
  padding: 16px; overflow-x: auto; white-space: pre-wrap; word-break: break-all;
}
.aa { display: inline-block; }
.liability { border-radius: 2px; cursor: help; }
.cdr-region { text-decoration: underline; text-decoration-thickness: 2px;
              text-underline-offset: 3px; }
.pos-label { color: #9ba5b0; user-select: none; font-size: 0.88em; }
.cdr-label { color: #1a5276; font-size: 0.78em; font-weight: 600;
             margin-left: 6px; }
table { width: 100%; border-collapse: collapse; font-size: 0.9em; }
thead th {
  background: #f4f5f7; padding: 10px 12px; text-align: left;
  font-weight: 600; color: #444; border-bottom: 2px solid #e0e3e8;
}
tbody tr { border-bottom: 1px solid #f0f0f0; }
tbody tr:hover { background: #fafbff; }
tbody td { padding: 8px 12px; vertical-align: middle; }
.pos-cell { font-family: monospace; color: #7f8c8d; font-weight: 600; }
.cdr-badge {
  display:inline-block; padding:1px 6px; border-radius:10px;
  background:#e3f2fd; color:#1565c0; font-size:0.75em;
  font-weight:700; margin-left:4px; white-space:nowrap;
}
.legend-table td { padding: 6px 10px; vertical-align: top; }
.legend-table tr { border-bottom: 1px solid #f0f0f0; }
.footer { text-align: center; margin-top: 40px; font-size: 0.82em; color: #aaa; }
"""


def _html_annotated_seq(seq: str, findings: list, cdrs: list, wrap: int = 60) -> str:
    n = len(seq)

    # Per-position liability map
    pos_f = [[] for _ in range(n)]
    for f in findings:
        for j in range(f["pos0"], min(f["pos0"] + f["span"], n)):
            pos_f[j].append(f)

    def best(flist):
        return min(flist, key=lambda f: RISK_ORDER.get(f["risk"], 9))

    # CDR position set for underline
    cdr_pos = {}
    for (cs, ce, cn) in cdrs:
        for p in range(cs, ce + 1):
            cdr_pos[p] = cn

    html_chars = []
    for i, aa in enumerate(seq):
        flist   = pos_f[i]
        in_cdr  = i in cdr_pos
        classes = "aa"
        style   = ""
        title   = ""

        if flist:
            bf      = best(flist)
            types   = " | ".join(sorted({f["short"] for f in flist}))
            notes   = "; ".join(filter(None, [f.get("note", "") for f in flist]))
            cdr_str = f" [{cdr_pos[i]}]" if in_cdr else ""
            title   = f"Pos {i+1}: {types}{cdr_str}" + (f" — {notes}" if notes else "")
            classes = "aa liability" + (" cdr-region" if in_cdr else "")
            style   = (
                f"background:{bf['color_bg']};color:{bf['color_text']};"
                f"border-bottom:2px solid {bf['color_border']};"
            )
        elif in_cdr:
            classes = "aa cdr-region"
            title   = f"Pos {i+1}: {cdr_pos[i]}"

        if title:
            html_chars.append(
                f'<span class="{classes}" style="{style}" title="{title}">{aa}</span>'
            )
        else:
            html_chars.append(f'<span class="aa">{aa}</span>')

    # Wrap with position labels; annotate CDR spans
    lines = []
    for start in range(0, n, wrap):
        chunk  = html_chars[start: start + wrap]
        end    = min(start + wrap, n)
        # Find CDR spans that overlap this line
        cdr_tags = []
        for (cs, ce, cn) in cdrs:
            os_ = max(cs, start)
            oe  = min(ce, end - 1)
            if os_ <= oe:
                cdr_tags.append(f'<span class="cdr-label">{cn}:{os_+1}–{oe+1}</span>')
        cdr_row = "".join(cdr_tags)
        lines.append(
            f'<span class="pos-label">{start+1:>5} </span>'
            + "".join(chunk)
            + f' <span class="pos-label">{end}</span>'
            + (f'  {cdr_row}' if cdr_row else "")
        )
    return "\n".join(lines)


def build_html_report(results: list) -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    CAT_ICONS = {
        "Deamidation": "🧪", "Oxidation": "⚡", "Isomerization": "🔄",
        "Glycosylation": "🍬", "Pyroglutamate": "🔵", "Truncation": "✂️",
    }
    CAT_ORDER = ["Deamidation", "Oxidation", "Isomerization",
                 "Glycosylation", "Pyroglutamate", "Truncation"]

    sections_html = []
    for name, seq, findings, summary, cdrs in results:
        # Summary cards
        cards = []
        for cat in CAT_ORDER:
            items  = summary["by_category"].get(cat, [])
            if not items and cat not in ("Deamidation", "Oxidation",
                                         "Isomerization", "Glycosylation",
                                         "Pyroglutamate", "Truncation"):
                continue
            high   = sum(1 for f in items if f["risk"] == "high")
            medium = sum(1 for f in items if f["risk"] == "medium")
            low    = sum(1 for f in items if f["risk"] == "low")
            info   = sum(1 for f in items if f["risk"] == "info")
            bars   = ""
            if high:   bars += f'<span class="risk-badge risk-high">{high} High</span> '
            if medium: bars += f'<span class="risk-badge risk-medium">{medium} Med</span> '
            if low:    bars += f'<span class="risk-badge risk-low">{low} Low</span> '
            if info:   bars += f'<span class="risk-badge risk-info">{info} Info</span>'
            cards.append(
                f'<div class="summary-card">'
                f'<div class="card-icon">{CAT_ICONS.get(cat,"•")}</div>'
                f'<div class="card-title">{cat}</div>'
                f'<div class="card-count">{len(items)}</div>'
                f'<div class="card-bars">{bars if bars else "—"}</div>'
                f'</div>'
            )

        # CDR summary line
        cdr_info = ""
        if cdrs:
            cdr_info = (
                "<p style='margin-bottom:14px;font-size:.88em;color:#1565c0;"
                "background:#e3f2fd;padding:8px 14px;border-radius:6px'>"
                "📐 <strong>CDR regions detected:</strong> "
                + ", ".join(f"<strong>{cn}</strong> ({cs+1}–{ce+1})"
                            for cs, ce, cn in cdrs)
                + "  <em style='color:#7f8c8d'>(Chothia-like approximation)</em></p>"
            )

        # Findings table
        sorted_f = sorted(findings, key=lambda f: (RISK_ORDER.get(f["risk"], 9), f["pos0"]))
        rows = []
        for f in sorted_f:
            rb = f'<span class="risk-badge risk-{f["risk"]}">{f["risk"].upper()}</span>'
            rs = (
                f'<span style="background:{f["color_bg"]};color:{f["color_text"]};'
                f'border:1px solid {f["color_border"]};padding:2px 5px;'
                f'border-radius:3px;font-family:monospace;font-weight:700">'
                f'{f["residues"][0]}</span>'
            )
            cdr_badge = (
                f'<span class="cdr-badge">{f["cdr"]}</span>'
                if f.get("cdr") else ""
            )
            note_text = (f.get("note") or "—") + (f' {cdr_badge}' if cdr_badge else "")
            rows.append(
                f"<tr><td class='pos-cell'>{f['pos1']}</td><td>{rs}</td>"
                f"<td>{f['category']}</td><td>{rb}</td>"
                f"<td style='font-size:.88em'>{f['display_label']}</td>"
                f"<td style='font-size:.82em;color:#555'>{note_text}</td></tr>"
            )

        high_ct    = len(summary["by_risk"].get("high", []))
        alert_bg   = "#fde8e8" if high_ct else "#d6eaf8"
        alert_col  = "#c0392b" if high_ct else "#154360"
        alert_text = (
            f"⚠️  {high_ct} high-risk liabilit{'y' if high_ct==1 else 'ies'} detected."
            if high_ct else "✅  No high-risk liabilities detected."
        )

        sections_html.append(f"""
        <div class="section">
          <div class="section-header">
            {name}
            <span style="font-weight:400;font-size:.88em;color:#888;margin-left:10px">
              {summary['seq_len']} residues · {summary['total']} findings
              {f'· {len(cdrs)} CDR(s) detected' if cdrs else ''}
            </span>
          </div>
          <div class="section-body">
            <div style="background:{alert_bg};color:{alert_col};padding:10px 16px;
                        border-radius:6px;margin-bottom:18px;font-weight:600;font-size:.9em">
              {alert_text}
            </div>
            {cdr_info}
            <div class="summary-grid">{"".join(cards)}</div>
            <div style="font-weight:600;margin-bottom:10px">Annotated Sequence
              <span style="font-weight:400;font-size:.8em;color:#888;margin-left:8px">
                hover for details · underlined = CDR region
              </span>
            </div>
            <div class="seq-block">{_html_annotated_seq(seq, findings, cdrs)}</div>
            <div style="font-weight:600;margin:20px 0 10px">Liabilities Detail</div>
            <div style="overflow-x:auto">
              <table>
                <thead><tr>
                  <th>Pos.</th><th>Residue(s)</th><th>Category</th>
                  <th>Risk</th><th>Liability Type</th><th>Note / CDR</th>
                </tr></thead>
                <tbody>{"".join(rows)}</tbody>
              </table>
            </div>
          </div>
        </div>
        """)

    # Legend
    legend_rows = []
    for key, lib in LIABILITIES.items():
        rb = (f'<span class="risk-badge risk-{lib["risk"]}">'
              f'{lib["risk"].upper()}</span>')
        legend_rows.append(
            f'<tr><td><span style="background:{lib["color_bg"]};color:{lib["color_text"]};'
            f'border:1px solid {lib["color_border"]};padding:2px 6px;border-radius:3px;'
            f'font-family:monospace;font-weight:600">{lib["short"]}</span></td>'
            f'<td>{lib["label"]}</td><td>{rb}</td>'
            f'<td style="color:#555;font-size:.85em">{lib["description"]}</td></tr>'
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>Protein Sequence Liability Analysis</title>
  <style>{CSS}</style>
</head>
<body>
<div class="page-wrap">
  <div class="report-header">
    <h1>🔬 Protein Sequence Liability Analysis</h1>
    <div class="subtitle">PTM risk · deamidation · oxidation · isomerization · glycosylation · pyroglutamate · truncation</div>
    <div class="meta">
      <span>📅 {now}</span>
      <span>🔢 {len(results)} chain(s)</span>
      <span>⚗️ Letarte Scientific Consulting</span>
    </div>
  </div>
  {"".join(sections_html)}
  <div class="section">
    <div class="section-header">📖 Liability Legend</div>
    <div class="section-body">
      <table class="legend-table">
        <thead><tr><th>Code</th><th>Liability</th><th>Risk</th><th>Description</th></tr></thead>
        <tbody>{"".join(legend_rows)}</tbody>
      </table>
      <p style="margin-top:16px;font-size:.82em;color:#888">
        <strong>CDR note:</strong> CDR boundaries are approximate (Chothia-like anchor residues).
        Tryptic cleavage sites (K/R) are excluded from truncation reporting.
      </p>
    </div>
  </div>
  <div class="footer">Protein Liability Analyzer · Letarte Scientific Consulting · For R&amp;D use only</div>
</div>
</body>
</html>"""


# ═══════════════════════════════════════════════════════════════════════════════
# ❷  COLOUR THEME & CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

THEME = {
    "bg_dark":    "#1a1a2e",
    "bg_mid":     "#16213e",
    "bg_accent":  "#0f3460",
    "bg_light":   "#f0f2f5",
    "bg_white":   "#ffffff",
    "bg_panel":   "#f8f9fa",
    "fg_white":   "#ffffff",
    "fg_dark":    "#1a1a2e",
    "fg_muted":   "#6b7280",
    "fg_subtle":  "#9ba5b0",
    "border":     "#e5e7eb",
    # ANALYZE button — dark royal blue
    "btn_analyze":      "#0d2d6b",
    "btn_analyze_hov":  "#1a4080",
}

PLACEHOLDER = (
    ">MyProtein_HC\n"
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIH...\n\n"
    "Accepts:\n"
    "  • FASTA format  (single or multi-chain)\n"
    "  • Raw amino-acid sequence\n"
    "  • Multiple FASTA entries in one paste"
)

CAT_ORDER = ["Deamidation", "Oxidation", "Isomerization",
             "Glycosylation", "Pyroglutamate", "Truncation"]
CAT_ICONS = {
    "Deamidation":   "🧪",
    "Oxidation":     "⚡",
    "Isomerization": "🔄",
    "Glycosylation": "🍬",
    "Pyroglutamate": "🔵",
    "Truncation":    "✂️",
}

# Liability key → (bg, fg) for tkinter text tags
TAG_COLORS = {k: (v["color_bg"], v["color_text"]) for k, v in LIABILITIES.items()}


# ═══════════════════════════════════════════════════════════════════════════════
# ❸  MAIN APPLICATION WINDOW
# ═══════════════════════════════════════════════════════════════════════════════

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Protein Sequence Liability Analyzer")
        self.geometry("1300x840")
        self.minsize(920, 620)
        self.configure(bg=THEME["bg_dark"])

        self._results            = []
        self._html_path          = None
        self._placeholder_active = True

        self._configure_styles()
        self._build_header()
        self._build_body()
        self._build_footer()

    # ── Styles ────────────────────────────────────────────────────────────────

    def _configure_styles(self):
        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("TFrame",        background=THEME["bg_light"])
        style.configure("Panel.TFrame",  background=THEME["bg_white"])
        style.configure("TLabel",
            background=THEME["bg_light"],
            foreground=THEME["fg_dark"],
            font=("Helvetica", 11))
        style.configure("TButton",
            font=("Helvetica", 10), padding=(10, 5))
        style.configure("Treeview",
            background=THEME["bg_white"],
            fieldbackground=THEME["bg_white"],
            foreground=THEME["fg_dark"],
            font=("Helvetica", 10),
            rowheight=26)
        style.configure("Treeview.Heading",
            background=THEME["bg_panel"],
            foreground=THEME["fg_dark"],
            font=("Helvetica", 10, "bold"))
        style.map("Treeview", background=[("selected", "#dbeafe")])

    # ── Header ────────────────────────────────────────────────────────────────

    def _build_header(self):
        hdr = tk.Frame(self, bg=THEME["bg_accent"], pady=14)
        hdr.pack(fill="x")
        left = tk.Frame(hdr, bg=THEME["bg_accent"])
        left.pack(side="left", padx=22)
        tk.Label(left,
                 text="🔬  Protein Sequence Liability Analyzer",
                 font=("Helvetica", 17, "bold"),
                 fg=THEME["fg_white"], bg=THEME["bg_accent"]).pack(anchor="w")
        tk.Label(left,
                 text="deamidation · oxidation · isomerization · glycosylation · pyroglutamate · truncation · CDR annotation",
                 font=("Helvetica", 9), fg="#a8b2c8", bg=THEME["bg_accent"]).pack(anchor="w")
        tk.Label(hdr, text="Letarte Scientific Consulting",
                 font=("Helvetica", 10, "italic"),
                 fg="#7b8fb5", bg=THEME["bg_accent"]).pack(side="right", padx=22)

    # ── Body ─────────────────────────────────────────────────────────────────

    def _build_body(self):
        body = tk.Frame(self, bg=THEME["bg_light"])
        body.pack(fill="both", expand=True)

        # ─ Left panel ────────────────────────────────────────────────────────
        left_outer = tk.Frame(body, bg=THEME["bg_light"], width=450)
        left_outer.pack(side="left", fill="both", padx=(14, 6), pady=14)
        left_outer.pack_propagate(False)

        left_card = tk.Frame(left_outer, bg=THEME["bg_white"],
                             highlightbackground=THEME["border"],
                             highlightthickness=1)
        left_card.pack(fill="both", expand=True)

        ch = tk.Frame(left_card, bg=THEME["bg_panel"], pady=10)
        ch.pack(fill="x")
        tk.Label(ch, text="Sequence Input", font=("Helvetica", 11, "bold"),
                 bg=THEME["bg_panel"], fg=THEME["fg_dark"]).pack(side="left", padx=14)

        tb = tk.Frame(left_card, bg=THEME["bg_white"], pady=8)
        tb.pack(fill="x", padx=10)
        ttk.Button(tb, text="📁  Load FASTA", command=self._load_fasta).pack(side="left", padx=(0, 6))
        ttk.Button(tb, text="🗑  Clear",      command=self._clear_input).pack(side="left")

        txt_frame = tk.Frame(left_card, bg=THEME["bg_panel"],
                             highlightbackground=THEME["border"],
                             highlightthickness=1)
        txt_frame.pack(fill="both", expand=True, padx=10, pady=(0, 8))

        self.seq_text = tk.Text(
            txt_frame,
            font=("Courier New", 11),
            bg=THEME["bg_panel"],
            fg=THEME["fg_muted"],
            insertbackground=THEME["fg_dark"],
            relief="flat",
            wrap="word",
            padx=10, pady=8,
        )
        vsb = ttk.Scrollbar(txt_frame, orient="vertical", command=self.seq_text.yview)
        self.seq_text.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        self.seq_text.pack(fill="both", expand=True)

        self.seq_text.insert("1.0", PLACEHOLDER)
        self.seq_text.bind("<FocusIn>",  self._on_focus_in)
        self.seq_text.bind("<FocusOut>", self._on_focus_out)

        info = tk.Frame(left_card, bg="#f0f7ff", pady=6)
        info.pack(fill="x", padx=10, pady=(0, 10))
        tk.Label(info,
                 text="ℹ  FASTA or raw sequence · single or multi-chain · CDR auto-detected",
                 font=("Helvetica", 9), bg="#f0f7ff", fg="#2563eb").pack(padx=8)

        # ── ANALYZE button — dark royal blue, prominent ───────────────────
        self.analyze_btn = tk.Button(
            left_outer,
            text="  ANALYZE  ",
            font=("Helvetica", 14, "bold"),
            bg=THEME["btn_analyze"],
            fg="#000000",
            activebackground=THEME["btn_analyze_hov"],
            activeforeground="#ffffff",
            relief="flat",
            cursor="hand2",
            padx=10, pady=14,
            command=self._run_analysis,
        )
        self.analyze_btn.pack(fill="x", pady=(8, 0))
        self.analyze_btn.bind("<Enter>",
            lambda e: self.analyze_btn.config(bg=THEME["btn_analyze_hov"]))
        self.analyze_btn.bind("<Leave>",
            lambda e: self.analyze_btn.config(bg=THEME["btn_analyze"]))

        # ─ Right panel ───────────────────────────────────────────────────────
        right_outer = tk.Frame(body, bg=THEME["bg_light"])
        right_outer.pack(side="right", fill="both", expand=True, padx=(6, 14), pady=14)

        self.notebook = ttk.Notebook(right_outer)
        self.notebook.pack(fill="both", expand=True)

        self._build_summary_tab()
        self._build_sequence_tab()
        self._build_findings_tab()

    # ── Summary tab ───────────────────────────────────────────────────────────

    def _build_summary_tab(self):
        tab = tk.Frame(self.notebook, bg=THEME["bg_white"])
        self.notebook.add(tab, text="  Summary  ")

        self.summary_welcome = tk.Frame(tab, bg=THEME["bg_white"])
        self.summary_welcome.pack(fill="both", expand=True)
        tk.Label(self.summary_welcome,
                 text="🔬\n\nPaste a sequence and click  ANALYZE",
                 font=("Helvetica", 14), justify="center",
                 bg=THEME["bg_white"], fg=THEME["fg_muted"]).pack(expand=True)

        self.summary_content = tk.Frame(tab, bg=THEME["bg_white"])

        self.alert_frame = tk.Frame(self.summary_content, bg=THEME["bg_white"])
        self.alert_frame.pack(fill="x", padx=14, pady=(14, 6))
        self.alert_label = tk.Label(self.alert_frame, text="",
                                     font=("Helvetica", 11, "bold"),
                                     anchor="w", padx=12, pady=8)
        self.alert_label.pack(fill="x")

        # CDR info banner
        self.cdr_label = tk.Label(self.summary_content, text="",
                                   font=("Helvetica", 10),
                                   bg="#e3f2fd", fg="#1565c0",
                                   anchor="w", padx=12, pady=5,
                                   wraplength=800, justify="left")
        self.cdr_label.pack(fill="x", padx=14, pady=(0, 6))

        sel_row = tk.Frame(self.summary_content, bg=THEME["bg_white"])
        sel_row.pack(fill="x", padx=14, pady=(4, 6))
        tk.Label(sel_row, text="Chain:", font=("Helvetica", 10),
                 bg=THEME["bg_white"], fg=THEME["fg_muted"]).pack(side="left")
        self.chain_var  = tk.StringVar()
        self.chain_menu = ttk.Combobox(sel_row, textvariable=self.chain_var,
                                        state="readonly", width=40)
        self.chain_menu.pack(side="left", padx=8)
        self.chain_menu.bind("<<ComboboxSelected>>", self._on_chain_change)

        self.cards_frame = tk.Frame(self.summary_content, bg=THEME["bg_white"])
        self.cards_frame.pack(fill="x", padx=10, pady=6)
        self._card_widgets = {}
        for cat in CAT_ORDER:
            card = tk.Frame(self.cards_frame, bg=THEME["bg_white"],
                            highlightbackground=THEME["border"],
                            highlightthickness=1, padx=10, pady=8)
            card.pack(side="left", fill="both", expand=True, padx=3)
            icon_lbl   = tk.Label(card, text=CAT_ICONS[cat],
                                   font=("Helvetica", 16), bg=THEME["bg_white"])
            count_lbl  = tk.Label(card, text="—",
                                   font=("Helvetica", 20, "bold"),
                                   bg=THEME["bg_white"], fg=THEME["fg_dark"])
            cat_lbl    = tk.Label(card, text=cat,
                                   font=("Helvetica", 8, "bold"),
                                   bg=THEME["bg_white"], fg=THEME["fg_muted"])
            detail_lbl = tk.Label(card, text="",
                                   font=("Helvetica", 8),
                                   bg=THEME["bg_white"], fg=THEME["fg_muted"],
                                   wraplength=85, justify="center")
            icon_lbl.pack()
            count_lbl.pack()
            cat_lbl.pack()
            detail_lbl.pack()
            self._card_widgets[cat] = (count_lbl, detail_lbl)

        self.stats_label = tk.Label(self.summary_content, text="",
                                     font=("Helvetica", 10),
                                     bg=THEME["bg_white"], fg=THEME["fg_muted"],
                                     anchor="w")
        self.stats_label.pack(fill="x", padx=14, pady=(4, 8))

    # ── Annotated sequence tab ────────────────────────────────────────────────

    def _build_sequence_tab(self):
        tab = tk.Frame(self.notebook, bg=THEME["bg_white"])
        self.notebook.add(tab, text="  Annotated Sequence  ")

        # Legend strip
        leg = tk.Frame(tab, bg=THEME["bg_panel"], pady=6)
        leg.pack(fill="x")
        tk.Label(leg, text="Colour key:", font=("Helvetica", 9, "bold"),
                 bg=THEME["bg_panel"], fg=THEME["fg_muted"]).pack(side="left", padx=10)
        for key, lib in LIABILITIES.items():
            bg, fg = lib["color_bg"], lib["color_text"]
            tk.Label(leg, text=f" {lib['short']} ",
                     font=("Courier New", 9, "bold"),
                     bg=bg, fg=fg, relief="flat", padx=2).pack(side="left", padx=2)
        # CDR underline sample
        tk.Label(leg, text=" ── ",
                 font=("Helvetica", 9), bg=THEME["bg_panel"],
                 fg=THEME["fg_muted"]).pack(side="left", padx=4)
        cdr_sample = tk.Label(leg, text=" CDR ",
                               font=("Courier New", 9),
                               bg=THEME["bg_panel"], fg="#1565c0")
        cdr_sample.pack(side="left", padx=2)
        tk.Label(leg, text="(underlined)",
                 font=("Helvetica", 9, "italic"),
                 bg=THEME["bg_panel"], fg=THEME["fg_muted"]).pack(side="left")

        txt_wrap = tk.Frame(tab, bg=THEME["bg_white"])
        txt_wrap.pack(fill="both", expand=True, padx=10, pady=10)

        self.ann_text = tk.Text(
            txt_wrap,
            font=("Courier New", 12),
            bg=THEME["bg_panel"],
            fg=THEME["fg_dark"],
            relief="flat",
            wrap="char",
            padx=12, pady=10,
            state="disabled",
        )
        vsb2 = ttk.Scrollbar(txt_wrap, orient="vertical", command=self.ann_text.yview)
        self.ann_text.configure(yscrollcommand=vsb2.set)
        vsb2.pack(side="right", fill="y")
        self.ann_text.pack(fill="both", expand=True)

        # Liability colour tags
        for key, (bg, fg) in TAG_COLORS.items():
            self.ann_text.tag_configure(
                key,
                background=bg,
                foreground=fg,
                font=("Courier New", 12, "bold"),
            )
        # CDR underline tag — no colour override, just underline
        self.ann_text.tag_configure(
            "cdr_underline",
            underline=True,
            font=("Courier New", 12),
        )
        # Combine: CDR + liability = bold + underline
        for key in TAG_COLORS:
            self.ann_text.tag_configure(
                f"{key}_cdr",
                background=TAG_COLORS[key][0],
                foreground=TAG_COLORS[key][1],
                font=("Courier New", 12, "bold"),
                underline=True,
            )
        self.ann_text.tag_configure("pos_label",
            foreground=THEME["fg_subtle"],
            font=("Courier New", 10))
        self.ann_text.tag_configure("cdr_header",
            foreground="#1565c0",
            font=("Courier New", 10, "bold"))

    # ── All Findings tab ──────────────────────────────────────────────────────

    def _build_findings_tab(self):
        tab = tk.Frame(self.notebook, bg=THEME["bg_white"])
        self.notebook.add(tab, text="  All Findings  ")

        fbar = tk.Frame(tab, bg=THEME["bg_panel"], pady=8)
        fbar.pack(fill="x")
        tk.Label(fbar, text="Risk:", font=("Helvetica", 10),
                 bg=THEME["bg_panel"], fg=THEME["fg_muted"]).pack(side="left", padx=10)
        self.filter_var = tk.StringVar(value="All")
        for opt in ["All", "High", "Medium", "Low", "Info"]:
            tk.Radiobutton(fbar, text=opt, variable=self.filter_var, value=opt,
                           command=self._apply_filter,
                           bg=THEME["bg_panel"], fg=THEME["fg_dark"],
                           selectcolor=THEME["bg_panel"],
                           font=("Helvetica", 10)).pack(side="left", padx=4)

        tk.Label(fbar, text="Category:", font=("Helvetica", 10),
                 bg=THEME["bg_panel"], fg=THEME["fg_muted"]).pack(side="left", padx=(16, 4))
        self.cat_filter_var  = tk.StringVar(value="All")
        self.cat_filter_menu = ttk.Combobox(
            fbar, textvariable=self.cat_filter_var,
            values=["All"] + CAT_ORDER, state="readonly", width=16)
        self.cat_filter_menu.pack(side="left")
        self.cat_filter_menu.bind("<<ComboboxSelected>>",
                                   lambda e: self._apply_filter())

        # CDR-only toggle
        self.cdr_only_var = tk.BooleanVar(value=False)
        tk.Checkbutton(fbar, text="CDR only", variable=self.cdr_only_var,
                       command=self._apply_filter,
                       bg=THEME["bg_panel"], fg=THEME["fg_dark"],
                       selectcolor=THEME["bg_panel"],
                       font=("Helvetica", 10)).pack(side="left", padx=(16, 4))

        cols   = ("pos", "residues", "category", "risk", "type", "cdr", "note")
        hdrs   = ("Pos.", "Residue", "Category", "Risk", "Liability Type", "CDR", "Note")
        widths = (55, 60, 110, 70, 255, 70, 275)

        tree_wrap = tk.Frame(tab, bg=THEME["bg_white"])
        tree_wrap.pack(fill="both", expand=True, padx=10, pady=8)

        self.findings_tree = ttk.Treeview(tree_wrap, columns=cols, show="headings")
        for col, hdr, w in zip(cols, hdrs, widths):
            self.findings_tree.heading(col, text=hdr,
                command=lambda c=col: self._sort_tree(c))
            anchor = "center" if col in ("pos", "residues", "risk", "cdr") else "w"
            self.findings_tree.column(col, width=w, anchor=anchor,
                                       stretch=(col == "note"))

        vsb3 = ttk.Scrollbar(tree_wrap, orient="vertical",
                              command=self.findings_tree.yview)
        hsb3 = ttk.Scrollbar(tree_wrap, orient="horizontal",
                              command=self.findings_tree.xview)
        self.findings_tree.configure(yscrollcommand=vsb3.set,
                                     xscrollcommand=hsb3.set)
        vsb3.pack(side="right", fill="y")
        hsb3.pack(side="bottom", fill="x")
        self.findings_tree.pack(fill="both", expand=True)

        self.findings_tree.tag_configure("high",   background="#fde8e8")
        self.findings_tree.tag_configure("medium", background="#fef0e6")
        self.findings_tree.tag_configure("low",    background="#e8f8f5")
        self.findings_tree.tag_configure("info",   background="#d6eaf8")

        self._all_findings_rows = []
        self._sort_col = "pos"
        self._sort_rev = False

    # ── Footer ────────────────────────────────────────────────────────────────

    def _build_footer(self):
        foot = tk.Frame(self, bg=THEME["bg_mid"], pady=8)
        foot.pack(fill="x", side="bottom")

        self.status_lbl = tk.Label(foot, text="Ready.",
                                    font=("Helvetica", 10),
                                    bg=THEME["bg_mid"], fg="#7b8fb5")
        self.status_lbl.pack(side="left", padx=14)

        self.open_btn = ttk.Button(foot, text="🌐  Open HTML Report",
                                    command=self._open_html)
        self.open_btn.pack(side="right", padx=8)
        self.open_btn.state(["disabled"])

        self.save_btn = ttk.Button(foot, text="💾  Save HTML Report",
                                    command=self._save_html)
        self.save_btn.pack(side="right", padx=4)
        self.save_btn.state(["disabled"])

        ttk.Separator(foot, orient="vertical").pack(side="right", fill="y", padx=8)

        self.count_lbl = tk.Label(foot, text="",
                                   font=("Helvetica", 10, "bold"),
                                   bg=THEME["bg_mid"], fg="#7b8fb5")
        self.count_lbl.pack(side="right", padx=8)

    # ── Placeholder ───────────────────────────────────────────────────────────

    def _on_focus_in(self, _event):
        if self._placeholder_active:
            self.seq_text.delete("1.0", "end")
            self.seq_text.config(fg=THEME["fg_dark"])
            self._placeholder_active = False

    def _on_focus_out(self, _event):
        if not self.seq_text.get("1.0", "end").strip():
            self.seq_text.config(fg=THEME["fg_muted"])
            self.seq_text.insert("1.0", PLACEHOLDER)
            self._placeholder_active = True

    def _clear_input(self):
        self.seq_text.config(fg=THEME["fg_muted"])
        self.seq_text.delete("1.0", "end")
        self.seq_text.insert("1.0", PLACEHOLDER)
        self._placeholder_active = True

    # ── File loading ──────────────────────────────────────────────────────────

    def _load_fasta(self):
        path = filedialog.askopenfilename(
            title="Open FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.faa *.txt"),
                       ("All files", "*.*")])
        if not path:
            return
        try:
            text = Path(path).read_text(encoding="utf-8", errors="replace")
        except Exception as exc:
            messagebox.showerror("Error", f"Could not read file:\n{exc}")
            return
        self.seq_text.delete("1.0", "end")
        self.seq_text.config(fg=THEME["fg_dark"])
        self.seq_text.insert("1.0", text.strip())
        self._placeholder_active = False
        self._set_status(f"Loaded: {Path(path).name}")

    # ── Analysis ─────────────────────────────────────────────────────────────

    def _run_analysis(self):
        if self._placeholder_active:
            messagebox.showwarning("No Sequence",
                                   "Please paste or load a protein sequence first.")
            return
        raw = self.seq_text.get("1.0", "end").strip()
        if not raw:
            messagebox.showwarning("No Sequence", "The sequence input is empty.")
            return

        self.analyze_btn.config(state="disabled", text="  Analyzing…  ")
        self._set_status("Analyzing…")

        def worker():
            pairs   = parse_fasta(raw) if ">" in raw else [("Sequence", raw)]
            if not pairs:
                self.after(0, lambda: messagebox.showerror(
                    "Parse Error", "No valid sequences found."))
                self.after(0, self._reset_btn)
                return
            results = process_sequences(pairs)
            self.after(0, lambda: self._display_results(results))

        threading.Thread(target=worker, daemon=True).start()

    def _reset_btn(self):
        self.analyze_btn.config(state="normal", text="  ANALYZE  ")

    def _display_results(self, results: list):
        self._results = results

        html     = build_html_report(results)
        tmp      = tempfile.NamedTemporaryFile(
            delete=False, suffix=".html",
            prefix="protein_liabilities_",
            mode="w", encoding="utf-8")
        tmp.write(html)
        tmp.close()
        self._html_path = tmp.name

        self._populate_summary(results)
        self._populate_sequence_tab(results)
        self._populate_findings_tab(results)

        self.save_btn.state(["!disabled"])
        self.open_btn.state(["!disabled"])

        total = sum(s["total"] for _, _, _, s, _ in results)
        highs = sum(len(s["by_risk"].get("high", [])) for _, _, _, s, _ in results)
        self.count_lbl.config(
            text=f"{total} liabilities  ·  {highs} high-risk")
        self._set_status(
            f"Done — {len(results)} chain(s), {total} liabilities flagged.")
        self._reset_btn()
        self.notebook.select(0)

    # ── Populate Summary ──────────────────────────────────────────────────────

    def _populate_summary(self, results: list):
        self.summary_welcome.pack_forget()
        self.summary_content.pack(fill="both", expand=True)

        names = [name for name, *_ in results]
        self.chain_menu["values"] = names
        self.chain_var.set(names[0])
        self._show_chain_summary(0)

    def _on_chain_change(self, _event=None):
        names = [name for name, *_ in self._results]
        idx   = names.index(self.chain_var.get()) if self.chain_var.get() in names else 0
        self._show_chain_summary(idx)

    def _show_chain_summary(self, idx: int):
        name, seq, findings, summary, cdrs = self._results[idx]

        high_ct = len(summary["by_risk"].get("high", []))
        self.alert_label.config(
            text=(f"⚠️   {high_ct} high-risk liabilit{'y' if high_ct==1 else 'ies'} detected"
                  if high_ct else "✅   No high-risk liabilities detected"),
            bg="#fde8e8" if high_ct else "#d6eaf8",
            fg="#c0392b" if high_ct else "#154360")

        if cdrs:
            cdr_str = "   ".join(f"{cn}: {cs+1}–{ce+1}" for cs, ce, cn in cdrs)
            self.cdr_label.config(
                text=f"📐  CDRs detected (Chothia approx.):  {cdr_str}",
                bg="#e3f2fd", fg="#1565c0")
            self.cdr_label.pack(fill="x", padx=14, pady=(0, 6))
        else:
            self.cdr_label.pack_forget()

        for cat in CAT_ORDER:
            items  = summary["by_category"].get(cat, [])
            counts = defaultdict(int)
            for f in items:
                counts[f["risk"]] += 1
            count_lbl, detail_lbl = self._card_widgets[cat]
            count_lbl.config(text=str(len(items)) if items else "0")
            parts = [f"{counts[r]} {r}" for r in ("high","medium","low","info")
                     if counts[r]]
            detail_lbl.config(text="\n".join(parts) if parts else "—")

        self.stats_label.config(
            text=(
                f"Length: {summary['seq_len']} aa  ·  Total: {summary['total']}  ·  "
                f"High: {len(summary['by_risk'].get('high',[]))}  "
                f"Med: {len(summary['by_risk'].get('medium',[]))}  "
                f"Low: {len(summary['by_risk'].get('low',[]))}  "
                f"Info: {len(summary['by_risk'].get('info',[]))}"
            )
        )

    # ── Populate annotated sequence ───────────────────────────────────────────

    def _populate_sequence_tab(self, results: list):
        WRAP = 60
        self.ann_text.config(state="normal")
        self.ann_text.delete("1.0", "end")

        for chain_idx, (name, seq, findings, summary, cdrs) in enumerate(results):
            if chain_idx > 0:
                self.ann_text.insert("end", "\n\n")

            # CDR header
            if cdrs:
                cdr_str = "  ".join(f"{cn}:{cs+1}-{ce+1}" for cs, ce, cn in cdrs)
                self.ann_text.insert(
                    "end",
                    f"▶  {name}   [{cdr_str}]\n",
                    "cdr_header")
            else:
                self.ann_text.insert("end", f"▶  {name}\n", "pos_label")

            n = len(seq)

            # Per-position liability map
            pos_f = [[] for _ in range(n)]
            for f in findings:
                for j in range(f["pos0"], min(f["pos0"] + f["span"], n)):
                    pos_f[j].append(f)

            def best_key(flist):
                if not flist:
                    return None
                return min(flist, key=lambda f: RISK_ORDER.get(f["risk"], 9))["type"]

            # CDR position set
            cdr_pos = set()
            for (cs, ce, _cn) in cdrs:
                cdr_pos.update(range(cs, ce + 1))

            for start in range(0, n, WRAP):
                chunk_end = min(start + WRAP, n)
                self.ann_text.insert("end", f"{start+1:>5}  ", "pos_label")
                for i in range(start, chunk_end):
                    aa    = seq[i]
                    flist = pos_f[i]
                    tag   = best_key(flist)
                    in_cdr = i in cdr_pos
                    if tag:
                        final_tag = f"{tag}_cdr" if in_cdr else tag
                        self.ann_text.insert("end", aa, final_tag)
                    elif in_cdr:
                        self.ann_text.insert("end", aa, "cdr_underline")
                    else:
                        self.ann_text.insert("end", aa)
                self.ann_text.insert("end", f"  {chunk_end}\n", "pos_label")

        self.ann_text.config(state="disabled")

    # ── Populate findings table ───────────────────────────────────────────────

    def _populate_findings_tab(self, results: list):
        self.findings_tree.delete(*self.findings_tree.get_children())
        self._all_findings_rows = []

        for name, seq, findings, summary, cdrs in results:
            for f in findings:
                cdr_val = f.get("cdr") or ""
                row = (
                    f["pos1"],
                    f["residues"][0],
                    f["category"],
                    f["risk"].upper(),
                    f["display_label"],        # ← no risk level in this column
                    cdr_val,                   # ← CDR column
                    f.get("note", "") or "",
                )
                self._all_findings_rows.append((row, f["risk"], name))

        self._sort_col = "pos"
        self._sort_rev = False
        self._apply_filter()

    def _apply_filter(self):
        risk_filter = self.filter_var.get().lower()
        cat_filter  = self.cat_filter_var.get()
        cdr_only    = self.cdr_only_var.get()

        self.findings_tree.delete(*self.findings_tree.get_children())
        for row, risk, chain in self._all_findings_rows:
            if risk_filter != "all" and risk != risk_filter:
                continue
            if cat_filter != "All" and row[2] != cat_filter:
                continue
            if cdr_only and not row[5]:          # col 5 = CDR
                continue
            self.findings_tree.insert("", "end", values=row, tags=(risk,))

    def _sort_tree(self, col: str):
        col_idx = {"pos":0,"residues":1,"category":2,"risk":3,
                   "type":4,"cdr":5,"note":6}[col]
        reverse = (self._sort_col == col) and not self._sort_rev
        self._sort_col = col
        self._sort_rev = reverse

        def key(item):
            v = item[0][col_idx]
            try:    return (0, int(v))
            except: return (1, str(v).lower())

        self._all_findings_rows.sort(key=key, reverse=reverse)
        self._apply_filter()

    # ── HTML report ───────────────────────────────────────────────────────────

    def _open_html(self):
        if self._html_path and os.path.exists(self._html_path):
            webbrowser.open(f"file://{self._html_path}")
        else:
            messagebox.showinfo("No Report", "Run an analysis first.")

    def _save_html(self):
        if not self._results:
            messagebox.showinfo("No Results", "Run an analysis first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".html",
            filetypes=[("HTML file", "*.html"), ("All files", "*.*")],
            initialfile="protein_liabilities.html",
            title="Save HTML Report")
        if not path:
            return
        html = build_html_report(self._results)
        Path(path).write_text(html, encoding="utf-8")
        self._set_status(f"Report saved → {path}")
        messagebox.showinfo("Saved", f"HTML report saved to:\n{path}")

    # ── Utility ───────────────────────────────────────────────────────────────

    def _set_status(self, msg: str):
        self.status_lbl.config(text=msg)


# ═══════════════════════════════════════════════════════════════════════════════
# ❹  ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    app = App()
    app.mainloop()
