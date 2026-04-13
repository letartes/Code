#!/usr/bin/env python3
"""
Protein Sequence Liability Analyzer
====================================
Identifies post-translational modification risks and sequence liabilities
in biopharmaceutical protein sequences.

Detected liabilities:
  - Deamidation   : NG (high), NS (high), NA/NT/NQ/QG/QS (medium)
  - Oxidation     : M/W (high), H/C (medium), F/Y (low)
  - Isomerization : DG/DS (high), DT/DA/DN/DD (medium)
  - Glycosylation : N-X-S/T sequon (N-linked), S-P / T-P (O-linked)
  - Truncation    : DP bond, N-term E/Q pyroglutamate, K/R tryptic sites,
                    C-term K clipping, N-term Met

Usage:
  python protein_liability_analyzer.py sequence.fasta
  python protein_liability_analyzer.py --sequence EVQLVESGGGLVQ...
  python protein_liability_analyzer.py --interactive
  python protein_liability_analyzer.py sequence.fasta --output report.html
"""

import re
import sys
import os
import argparse
from collections import defaultdict
from datetime import datetime
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Liability catalogue
# ─────────────────────────────────────────────────────────────────────────────

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
            "N-terminal Met is subject to Met aminopeptidase removal and "
            "oxidation."
        ),
        "references": "Dick et al. (2008); Harris (1995)",
        "span_length": 1,
    },
}

# Risk level ordering for sorting
RISK_ORDER = {"high": 0, "medium": 1, "low": 2, "info": 3}


# ─────────────────────────────────────────────────────────────────────────────
# Sequence parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_fasta(text: str) -> list[tuple[str, str]]:
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
            # Strip whitespace and digits (common in numbered FASTA)
            current_seq.append(re.sub(r"[^A-Za-z]", "", line))
    if current_seq:
        sequences.append((current_name, "".join(current_seq).upper()))
    return sequences


def clean_sequence(seq: str) -> str:
    """Remove whitespace and digits; uppercase."""
    return re.sub(r"[^A-Za-z]", "", seq).upper()


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def validate_sequence(seq: str) -> tuple[bool, list[str]]:
    """Return (is_valid, list_of_warnings)."""
    warnings = []
    invalid = sorted(set(seq) - VALID_AA)
    if invalid:
        warnings.append(f"Non-standard characters found and ignored: {', '.join(invalid)}")
    if len(seq) < 5:
        warnings.append("Sequence is very short (< 5 residues).")
    return len(invalid) == 0, warnings


# ─────────────────────────────────────────────────────────────────────────────
# Detection logic
# ─────────────────────────────────────────────────────────────────────────────

def find_liabilities(seq: str, n_term_window: int = 10) -> list[dict]:
    """
    Scan sequence and return a list of finding dicts:
      { type, position (1-based), span, residues, label, ... }
    """
    findings = []
    n = len(seq)

    def add(liability_key: str, pos0: int, span: int, note: str = ""):
        """pos0 is 0-based start."""
        findings.append({
            "type": liability_key,
            "pos0": pos0,
            "span": span,
            "pos1": pos0 + 1,          # 1-based
            "residues": seq[pos0: pos0 + span],
            "note": note,
            **{k: v for k, v in LIABILITIES[liability_key].items()
               if k not in ("description", "references", "span_length")},
        })

    # ── Deamidation ───────────────────────────────────────────────────────────
    HIGH_DEAMID  = {"NG", "NS"}
    MED_DEAMID   = {"NA", "NT", "NQ", "NR", "NH", "QG", "QS"}
    for i in range(n - 1):
        pair = seq[i:i+2]
        if pair in HIGH_DEAMID:
            add("deamid_high", i, 2)
        elif pair in MED_DEAMID:
            add("deamid_medium", i, 2)

    # ── Oxidation ─────────────────────────────────────────────────────────────
    for i, aa in enumerate(seq):
        if aa in "MW":
            note = ""
            if aa == "M" and i == 0:
                note = "N-terminal Met (also truncation risk)"
            add("oxid_high", i, 1, note)
        elif aa in "HC":
            add("oxid_medium", i, 1)
        elif aa in "FY":
            add("oxid_low", i, 1)

    # ── Isomerization ─────────────────────────────────────────────────────────
    HIGH_ISOM  = {"DG", "DS"}
    MED_ISOM   = {"DT", "DA", "DN", "DD", "DH"}
    for i in range(n - 1):
        pair = seq[i:i+2]
        if pair in HIGH_ISOM:
            add("isom_high", i, 2)
        elif pair in MED_ISOM:
            add("isom_medium", i, 2)

    # ── N-Glycosylation sequon NXS/T (X ≠ P) ─────────────────────────────────
    for i in range(n - 2):
        if seq[i] == "N" and seq[i+1] != "P" and seq[i+2] in "ST":
            add("n_glycan", i, 3)

    # ── O-Glycosylation (S/T–P) ───────────────────────────────────────────────
    for i in range(n - 1):
        if seq[i] in "ST" and seq[i+1] == "P":
            add("o_glycan", i, 2)

    # ── Truncation: Asp–Pro (acid-labile) ─────────────────────────────────────
    for i in range(n - 1):
        if seq[i:i+2] == "DP":
            add("trunc_high", i, 2, "Acid-labile Asp–Pro bond")

    # ── Truncation: N-term pyroglutamate (Q or E at position 1) ──────────────
    if seq and seq[0] in "QE":
        note = (
            "N-terminal Gln → pyroglutamate (cyclisation, –17 Da)"
            if seq[0] == "Q"
            else "N-terminal Glu → pyroglutamate (cyclisation, –18 Da)"
        )
        add("trunc_high", 0, 1, note)

    # ── Truncation: tryptic sites K / R (not K-P or R-P) ────────────────────
    for i, aa in enumerate(seq):
        if aa in "KR":
            # C-terminal K clipping is a well-known mAb liability
            if i == n - 1 and aa == "K":
                add("trunc_medium", i, 1, "C-terminal Lys clipping (CHO expression)")
            else:
                next_aa = seq[i+1] if i + 1 < n else ""
                if next_aa != "P":
                    add("trunc_medium", i, 1, f"Tryptic site ({aa})")
                # K-P / R-P are trypsin-resistant, still flag lightly
                else:
                    add("trunc_medium", i, 1, f"Trypsin-resistant {aa}–Pro (low risk)")

    # ── N-terminal Met ────────────────────────────────────────────────────────
    if seq and seq[0] == "M":
        add("trunc_medium", 0, 1,
            "N-terminal Met: subject to Met aminopeptidase removal")

    return findings


# ─────────────────────────────────────────────────────────────────────────────
# Statistics helpers
# ─────────────────────────────────────────────────────────────────────────────

def summarize(findings: list[dict], seq_len: int) -> dict:
    by_category = defaultdict(list)
    by_risk     = defaultdict(list)
    by_type     = defaultdict(list)
    for f in findings:
        by_category[f["category"]].append(f)
        by_risk[f["risk"]].append(f)
        by_type[f["type"]].append(f)
    return {
        "total": len(findings),
        "by_category": dict(by_category),
        "by_risk": dict(by_risk),
        "by_type": dict(by_type),
        "seq_len": seq_len,
    }


# ─────────────────────────────────────────────────────────────────────────────
# HTML report generation
# ─────────────────────────────────────────────────────────────────────────────

def build_annotated_sequence(seq: str, findings: list[dict], wrap: int = 60) -> str:
    """Build HTML for colour-annotated sequence, wrapped at `wrap` residues."""
    n = len(seq)

    # Build a per-position list of findings (a residue can have multiple)
    pos_findings: list[list[dict]] = [[] for _ in range(n)]
    for f in findings:
        for j in range(f["pos0"], min(f["pos0"] + f["span"], n)):
            pos_findings[j].append(f)

    # For each position, pick the highest-priority finding for the background
    # colour (ties broken by first encounter).
    def best_finding(flist: list[dict]) -> dict | None:
        if not flist:
            return None
        return min(flist, key=lambda f: (RISK_ORDER.get(f["risk"], 9),))

    # Build HTML character by character, creating <span> elements
    html_chars = []
    for i, aa in enumerate(seq):
        flist = pos_findings[i]
        if not flist:
            html_chars.append(f'<span class="aa">{aa}</span>')
        else:
            bf = best_finding(flist)
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

    # Wrap into lines with position numbers
    lines_html = []
    for start in range(0, n, wrap):
        chunk = html_chars[start: start + wrap]
        pos_label = f'<span class="pos-label">{start+1:>5} </span>'
        end_label  = f' <span class="pos-label">{min(start+wrap, n)}</span>'
        lines_html.append(pos_label + "".join(chunk) + end_label)

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


def findings_table_html(findings: list[dict]) -> str:
    if not findings:
        return "<p><em>No liabilities found.</em></p>"

    sorted_findings = sorted(
        findings,
        key=lambda f: (RISK_ORDER.get(f["risk"], 9), f["pos0"])
    )

    rows = []
    for f in sorted_findings:
        risk_badge = (
            f'<span class="risk-badge risk-{f["risk"]}">{f["risk"].upper()}</span>'
        )
        res_span = (
            f'<span style="background:{f["color_bg"]};color:{f["color_text"]};'
            f'border:1px solid {f["color_border"]};padding:2px 5px;'
            f'border-radius:3px;font-family:monospace;font-weight:700">'
            f'{f["residues"]}</span>'
        )
        note = f["note"] if f["note"] else "—"
        rows.append(
            f"<tr>"
            f"<td class='pos-cell'>{f['pos1']}</td>"
            f"<td>{res_span}</td>"
            f"<td>{f['category']}</td>"
            f"<td>{risk_badge}</td>"
            f"<td style='font-size:0.88em'>{f['label']}</td>"
            f"<td style='font-size:0.82em;color:#555'>{note}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def summary_cards_html(summary: dict) -> str:
    order = ["Deamidation", "Oxidation", "Isomerization", "Glycosylation", "Truncation"]
    cards = []
    for cat in order:
        items = summary["by_category"].get(cat, [])
        high   = sum(1 for f in items if f["risk"] == "high")
        medium = sum(1 for f in items if f["risk"] == "medium")
        low    = sum(1 for f in items if f["risk"] == "low")
        info   = sum(1 for f in items if f["risk"] == "info")
        total  = len(items)
        bars = ""
        if high:   bars += f'<span class="risk-badge risk-high">  {high} High</span> '
        if medium: bars += f'<span class="risk-badge risk-medium">{medium} Med</span> '
        if low:    bars += f'<span class="risk-badge risk-low">   {low} Low</span> '
        if info:   bars += f'<span class="risk-badge risk-info">  {info} Info</span>'
        icon = {
            "Deamidation": "🧪",
            "Oxidation": "⚡",
            "Isomerization": "🔄",
            "Glycosylation": "🍬",
            "Truncation": "✂️",
        }.get(cat, "•")
        cards.append(
            f'<div class="summary-card">'
            f'<div class="card-icon">{icon}</div>'
            f'<div class="card-title">{cat}</div>'
            f'<div class="card-count">{total}</div>'
            f'<div class="card-bars">{bars}</div>'
            f'</div>'
        )
    return "\n".join(cards)


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

/* Cards */
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
.card-icon { font-size: 1.8em; margin-bottom: 6px; }
.card-title { font-weight: 600; font-size: 0.9em; color: #555; margin-bottom: 4px; }
.card-count { font-size: 2em; font-weight: 700; color: #1a1a2e; margin-bottom: 8px; }
.card-bars { display:flex; flex-wrap:wrap; gap:4px; justify-content:center; }

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
  font-size: 0.92em;
  line-height: 2.0;
  letter-spacing: 0.05em;
  background: #f8f9fa;
  border-radius: 8px;
  padding: 16px;
  overflow-x: auto;
  white-space: pre-wrap;
  word-break: break-all;
}
.aa { display: inline-block; }
.liability { border-radius: 2px; cursor: help; }
.pos-label { color: #9ba5b0; user-select: none; font-size: 0.88em; }

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

/* Footer */
.footer {
  text-align: center;
  margin-top: 40px;
  font-size: 0.82em;
  color: #aaa;
}
.no-print { }
@media print {
  .no-print { display: none; }
  body { background: #fff; }
  .section { box-shadow: none; border: 1px solid #ddd; }
}
"""


def build_html_report(
    sequences: list[tuple[str, list[dict], dict]],
    title: str = "Protein Sequence Liability Analysis",
    generated_by: str = "Protein Liability Analyzer",
) -> str:
    """
    sequences: list of (name, findings, summary) tuples
    """
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    seq_sections = []
    for idx, (name, seq, findings, summary) in enumerate(sequences):
        annotated = build_annotated_sequence(seq, findings)
        table_rows = findings_table_html(findings)
        cards = summary_cards_html(summary)
        total = summary["total"]
        high_ct = len(summary["by_risk"].get("high", []))

        alert_color = "#c0392b" if high_ct > 0 else "#2e86c1"
        alert_bg    = "#fde8e8" if high_ct > 0 else "#d6eaf8"
        alert_text  = (
            f"⚠️  {high_ct} high-risk liabilit{'y' if high_ct == 1 else 'ies'} detected."
            if high_ct > 0
            else "✅  No high-risk liabilities detected."
        )

        seq_sections.append(f"""
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

    html = f"""<!DOCTYPE html>
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
    </div>
    <div class="meta">
      <span>📅 Generated: {now}</span>
      <span>🔢 Chains analysed: {len(sequences)}</span>
      <span>⚗️  {generated_by}</span>
    </div>
  </div>

  {''.join(seq_sections)}

  <div class="section">
    <div class="section-header">📖 Liability Legend &amp; References</div>
    <div class="section-body">
      <div style="overflow-x:auto">
        <table class="legend-table">
          <thead>
            <tr>
              <th>Code</th><th>Liability</th><th>Risk</th><th>Description</th>
            </tr>
          </thead>
          <tbody>
            {legend_rows}
          </tbody>
        </table>
      </div>
      <div style="margin-top:20px;font-size:0.82em;color:#888;line-height:1.6">
        <strong>References:</strong> Geiger &amp; Clarke (1987) <em>J Biol Chem</em>;
        Robinson &amp; Robinson (2001) <em>PNAS</em>; Lam et al. (1997) <em>J Pharm Sci</em>;
        Clarke (1987) <em>J Biol Chem</em>; Piszkiewicz et al. (1970) <em>BBRC</em>;
        Dick et al. (2008) <em>Biotechnol Bioeng</em>; Harris (1995) <em>J Chromatogr</em>;
        Apweiler et al. (1999) <em>Biochim Biophys Acta</em>.
      </div>
    </div>
  </div>

  <div class="footer">
    Generated by Protein Liability Analyzer · Letarte Scientific Consulting ·
    For research and development use only
  </div>

</div>
</body>
</html>"""
    return html


# ─────────────────────────────────────────────────────────────────────────────
# Console (text) output
# ─────────────────────────────────────────────────────────────────────────────

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
}

RISK_ANSI = {
    "high":   ANSI["red"],
    "medium": ANSI["orange"],
    "low":    ANSI["green"],
    "info":   ANSI["blue"],
}


def print_console_report(
    name: str,
    seq: str,
    findings: list[dict],
    summary: dict,
    use_color: bool = True,
) -> None:
    C = ANSI if use_color else {k: "" for k in ANSI}
    R = RISK_ANSI if use_color else {k: "" for k in RISK_ANSI}

    sep = "─" * 72
    print(f"\n{C['bold']}{sep}{C['reset']}")
    print(f"{C['bold']}  {name}{C['reset']}   ({summary['seq_len']} residues)")
    print(f"{sep}")

    # Category summary
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

    print(f"\n{C['bold']}  Detailed Findings:{C['reset']}")
    sorted_f = sorted(findings, key=lambda f: (RISK_ORDER.get(f["risk"], 9), f["pos0"]))
    for f in sorted_f:
        risk_str = f"{R[f['risk']]}{f['risk'].upper():6}{C['reset']}"
        note = f"  ← {f['note']}" if f["note"] else ""
        print(
            f"  Pos {f['pos1']:>5}  {C['bold']}{f['residues']:<5}{C['reset']}"
            f"  {risk_str}  {f['category']:<16}  {f['label']}{note}"
        )

    print(f"{sep}\n")


# ─────────────────────────────────────────────────────────────────────────────
# Main entry-point
# ─────────────────────────────────────────────────────────────────────────────

def process_sequences(pairs: list[tuple[str, str]]) -> list[tuple]:
    """Run analysis on (name, raw_seq) pairs; return (name, seq, findings, summary)."""
    results = []
    for name, raw_seq in pairs:
        seq = clean_sequence(raw_seq)
        _, warns = validate_sequence(seq)
        for w in warns:
            print(f"  [Warning] {name}: {w}", file=sys.stderr)
        findings = find_liabilities(seq)
        summary  = summarize(findings, len(seq))
        results.append((name, seq, findings, summary))
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Protein Sequence Liability Analyzer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "fasta_file",
        nargs="?",
        help="FASTA file containing one or more protein sequences",
    )
    group.add_argument(
        "--sequence", "-s",
        metavar="SEQ",
        help="Raw amino-acid sequence string (single-letter code)",
    )
    group.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Prompt for sequence(s) interactively",
    )
    parser.add_argument(
        "--output", "-o",
        metavar="FILE",
        default=None,
        help="Output HTML report path (default: <input_basename>_liabilities.html)",
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable ANSI colours in console output",
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Skip HTML report generation",
    )
    args = parser.parse_args()

    pairs: list[tuple[str, str]] = []
    out_stem = "protein"

    if args.fasta_file:
        path = Path(args.fasta_file)
        if not path.exists():
            print(f"Error: file not found: {path}", file=sys.stderr)
            sys.exit(1)
        text = path.read_text(encoding="utf-8", errors="replace")
        pairs = parse_fasta(text) if text.lstrip().startswith(">") else [
            (path.stem, text)
        ]
        out_stem = path.stem

    elif args.sequence:
        pairs = [("Input Sequence", args.sequence)]

    elif args.interactive:
        print("Protein Sequence Liability Analyzer")
        print("─" * 40)
        print("Paste one or more sequences (FASTA format or raw sequence).")
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
        text = "\n".join(lines)
        pairs = parse_fasta(text) if ">" in text else [("Input Sequence", text)]

    else:
        # Read from stdin
        if sys.stdin.isatty():
            print("Usage: protein_liability_analyzer.py [fasta_file] [--sequence SEQ]")
            print("       protein_liability_analyzer.py --interactive")
            print("       protein_liability_analyzer.py --help")
            sys.exit(0)
        text = sys.stdin.read()
        pairs = parse_fasta(text) if ">" in text else [("stdin", text)]

    if not pairs:
        print("No sequences found.", file=sys.stderr)
        sys.exit(1)

    use_color = not args.no_color and sys.stdout.isatty()
    results = process_sequences(pairs)

    # Console output
    for name, seq, findings, summary in results:
        print_console_report(name, seq, findings, summary, use_color=use_color)

    # HTML output
    if not args.no_html:
        if args.output:
            html_path = Path(args.output)
        else:
            html_path = Path(f"{out_stem}_liabilities.html")

        html = build_html_report(results)
        html_path.write_text(html, encoding="utf-8")
        print(f"📄  HTML report saved → {html_path}")


if __name__ == "__main__":
    main()
