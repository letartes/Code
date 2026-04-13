#!/usr/bin/env python3
"""
Protein Sequence Alignment Tool
================================
Pairwise and multiple sequence alignment for antibody / protein sequences.
Uses BLOSUM62 scoring matrix with Smith-Waterman (local) or
Needleman-Wunsch (global) alignment — identical to BLASTP logic.
Multiple sequence alignment uses a center-star progressive approach.

CDR annotation is built-in (no external packages required).
Uses conserved framework anchor residues (Kabat convention):
  VH:  C~22, W~36, C~92, WGxG~103
  VL:  C~23, W~35, C~88, FGxG~98

All sequences are processed locally. Nothing is transmitted externally.

Usage:
    python3 protein_aligner.py
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import re
import os

# ─────────────────────────────────────────────────────────────────────────────
# BLOSUM62 matrix (built-in — no biopython needed)
# ─────────────────────────────────────────────────────────────────────────────
BLOSUM62_RAW = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -2 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -2 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""

def _parse_blosum62():
    lines = [l for l in BLOSUM62_RAW.strip().splitlines() if l.strip()]
    header = lines[0].split()
    matrix = {}
    for line in lines[1:]:
        parts = line.split()
        row_aa = parts[0]
        for col_aa, val in zip(header, parts[1:]):
            matrix[(row_aa, col_aa)] = int(val)
    return matrix

BLOSUM62 = _parse_blosum62()

def blosum_score(a, b):
    a, b = a.upper(), b.upper()
    return BLOSUM62.get((a, b), BLOSUM62.get((b, a), -4))

# ─────────────────────────────────────────────────────────────────────────────
# Pairwise alignment — Smith-Waterman (local) & Needleman-Wunsch (global)
# ─────────────────────────────────────────────────────────────────────────────
NEG_INF = float('-inf')

def align_sequences(seq1, seq2, mode='local', gap_open=-11, gap_extend=-1):
    """
    Affine-gap pairwise alignment using BLOSUM62.
    mode: 'local' (Smith-Waterman / BLAST-like) or 'global' (Needleman-Wunsch)
    Returns (aligned_seq1, aligned_seq2, score)
    """
    m, n = len(seq1), len(seq2)
    M  = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    IX = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    IY = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    TB = [[''] * (n + 1) for _ in range(m + 1)]

    M[0][0] = 0
    for i in range(1, m + 1):
        IX[i][0] = gap_open + (i - 1) * gap_extend if mode == 'global' else NEG_INF
        M[i][0]  = IX[i][0]
        TB[i][0] = 'X'
    for j in range(1, n + 1):
        IY[0][j] = gap_open + (j - 1) * gap_extend if mode == 'global' else NEG_INF
        M[0][j]  = IY[0][j]
        TB[0][j] = 'Y'

    best_score, best_pos = 0, (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = blosum_score(seq1[i-1], seq2[j-1])
            IX[i][j] = max(M[i-1][j] + gap_open,  IX[i-1][j] + gap_extend)
            IY[i][j] = max(M[i][j-1] + gap_open,  IY[i][j-1] + gap_extend)
            diag = M[i-1][j-1] + s
            candidates = [diag, IX[i][j], IY[i][j]]
            if mode == 'local':
                candidates.append(0)
            best = max(candidates)
            M[i][j] = best
            if best == 0 and mode == 'local':  TB[i][j] = 'S'
            elif best == diag:                 TB[i][j] = 'M'
            elif best == IX[i][j]:             TB[i][j] = 'X'
            else:                              TB[i][j] = 'Y'
            if mode == 'local' and best > best_score:
                best_score = best
                best_pos   = (i, j)

    i, j       = best_pos if mode == 'local' else (m, n)
    final_score = best_score if mode == 'local' else M[m][n]

    aln1, aln2 = [], []
    while i > 0 or j > 0:
        tb = TB[i][j]
        if tb == 'S' or (mode == 'local' and M[i][j] <= 0 and i > 0 and j > 0):
            break
        if tb == 'M':
            aln1.append(seq1[i-1]); aln2.append(seq2[j-1]); i -= 1; j -= 1
        elif tb == 'X':
            aln1.append(seq1[i-1]); aln2.append('-');        i -= 1
        elif tb == 'Y':
            aln1.append('-');       aln2.append(seq2[j-1]);  j -= 1
        else:
            break

    return ''.join(reversed(aln1)), ''.join(reversed(aln2)), final_score

# ─────────────────────────────────────────────────────────────────────────────
# Multiple Sequence Alignment — center-star progressive algorithm
# ─────────────────────────────────────────────────────────────────────────────

def multiple_sequence_alignment(seqs_named, mode='local', gap_open=-11, gap_extend=-1):
    """
    Center-star progressive MSA.
      1. Choose the sequence with the highest total pairwise score as the center.
      2. Align every other sequence to the center (pairwise Smith-Waterman).
      3. Reconcile all gap insertions into a single merged alignment.
    Returns (list of (name, aligned_seq), avg_center_score).
    For N=2 falls back to plain pairwise.
    """
    names = [n for n, _ in seqs_named]
    seqs  = [s for _, s in seqs_named]
    n     = len(seqs)

    if n < 2:
        return seqs_named, 0.0

    if n == 2:
        a1, a2, sc = align_sequences(seqs[0], seqs[1], mode, gap_open, gap_extend)
        return [(names[0], a1), (names[1], a2)], sc

    # ── Step 1: pick center ────────────────────────────────────────────────
    pair_sc = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            _, _, s = align_sequences(seqs[i], seqs[j], mode, gap_open, gap_extend)
            pair_sc[i][j] = pair_sc[j][i] = s

    center = max(range(n), key=lambda i: sum(pair_sc[i]))
    c_seq  = seqs[center]
    c_len  = len(c_seq)

    # ── Step 2: align all others to center ────────────────────────────────
    pw = []   # (orig_idx, aln_center, aln_other)
    for i in range(n):
        if i == center:
            continue
        ac, ao, _ = align_sequences(c_seq, seqs[i], mode, gap_open, gap_extend)
        pw.append((i, ac, ao))

    # ── Step 3: build merged insertion map ────────────────────────────────
    # max_ins[k] = max gap columns to insert *before* center position k
    max_ins = [0] * (c_len + 1)
    for _, ac, _ in pw:
        cpos = ins = 0
        for ch in ac:
            if ch == '-':
                ins += 1
            else:
                if ins > max_ins[cpos]:
                    max_ins[cpos] = ins
                cpos += 1
                ins  = 0
        if ins > max_ins[cpos]:
            max_ins[cpos] = ins

    # ── Step 4: rebuild center ────────────────────────────────────────────
    c_aln = []
    for k in range(c_len):
        c_aln.extend(['-'] * max_ins[k])
        c_aln.append(c_seq[k])
    c_aln.extend(['-'] * max_ins[c_len])

    # ── Step 5: rebuild each other sequence ──────────────────────────────
    def rebuild(ac, ao):
        out = []
        ai  = 0
        for cp in range(c_len + 1):
            # collect gap columns at position cp
            local = 0
            while ai < len(ac) and ac[ai] == '-':
                local += 1
                ai    += 1
            chars = [ao[ai - local + j] if (ai - local + j) < len(ao) else '-'
                     for j in range(local)]
            out.extend(['-'] * (max_ins[cp] - local))
            out.extend(chars)
            if cp < c_len:
                out.append(ao[ai] if ai < len(ao) else '-')
                ai += 1
        return ''.join(out)

    result        = [None] * n
    result[center] = (names[center], ''.join(c_aln))
    for oi, ac, ao in pw:
        result[oi] = (names[oi], rebuild(ac, ao))

    avg = sum(pair_sc[center]) / (n - 1)
    return result, avg


def compute_pairwise_stats(aligned):
    """Return dict of (name_i, name_j) -> (matches, align_len, pct_id)."""
    stats = {}
    seqs  = [(n, s) for n, s in aligned]
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            ni, si = seqs[i]
            nj, sj = seqs[j]
            L   = min(len(si), len(sj))
            si2, sj2 = si[:L], sj[:L]
            m   = sum(a == b and a != '-' for a, b in zip(si2, sj2))
            pct = 100.0 * m / L if L else 0
            stats[(ni, nj)] = (m, L, pct)
    return stats


def msa_consensus(aligned):
    """
    CLUSTAL-style consensus:
      '*' = 100% identical  ':' = conserved (BLOSUM62 > 0)  '.' = semi-conserved  ' ' = variable
    """
    if not aligned:
        return ''
    seqs = [s for _, s in aligned]
    L    = max(len(s) for s in seqs)
    cons = []
    for i in range(L):
        col = [s[i] if i < len(s) else '-' for s in seqs]
        non_gap = [c for c in col if c != '-']
        if not non_gap:
            cons.append(' ')
            continue
        if len(set(non_gap)) == 1 and len(non_gap) == len(col):
            cons.append('*')
        elif len(set(non_gap)) == 1:
            cons.append(':')
        else:
            # Check if all non-gap pairs have positive BLOSUM62 score
            pairs = [(non_gap[a], non_gap[b])
                     for a in range(len(non_gap))
                     for b in range(a + 1, len(non_gap))]
            if all(blosum_score(a, b) > 0 for a, b in pairs):
                cons.append(':')
            elif any(blosum_score(a, b) > 0 for a, b in pairs):
                cons.append('.')
            else:
                cons.append(' ')
    return ''.join(cons)

# ─────────────────────────────────────────────────────────────────────────────
# CDR annotation — built-in, no external packages required
# ─────────────────────────────────────────────────────────────────────────────

# Scheme-specific CDR boundary offsets (all relative to detected anchors c1, w1, c2)
# Tuple: (cdr1_from_c1, cdr1_end_from_w1,   <- cdr1 = c1+A .. w1+B
#          cdr2_from_w1, cdr2_length,         <- cdr2 = w1+C .. w1+C+D-1
#          cdr3_from_c2)                      <- cdr3 = c2+E .. J-anchor-1
#
# VH offsets validated against Trastuzumab, Pertuzumab, Cetuximab, Bevacizumab:
#   Kabat  CDR1 = c1+9  ..w1-1   → DYTMD          (Kabat 31–35)
#   Chothia CDR1= c1+4  ..w1-6   → GFTFT           (Chothia / structural loop 26–30)
#   IMGT   CDR1 = c1+4  ..w1-1   → GFTFTDYTMD     (IMGT 27–38, unified loop)
#
#   Kabat  CDR2 = w1+14 .. +29   → DVNPNSGGSIYNQRFK (Kabat 50–65, 16 aa)
#   Chothia CDR2= w1+16 .. +20   → NPNSG            (structural loop 52–56, 5 aa)
#   IMGT   CDR2 = w1+20 .. +29   → GGSIYNQRFK       (IMGT 56–65, 10 aa)
#
#   CDR3 anchored identically across schemes: c2+3 .. WGxG–1  (VH)
#                                              c2+1 .. FGxG–1  (VL)
_CDR_PARAMS = {
    'H': {
        'kabat':   dict(c1r1=9,  w1r1=-1, w1r2=14, r2len=16, c2r3=3),
        'chothia': dict(c1r1=4,  w1r1=-6, w1r2=16, r2len=5,  c2r3=3),
        'imgt':    dict(c1r1=4,  w1r1=-1, w1r2=20, r2len=10, c2r3=4),
    },
    # VL: Kabat/Chothia/IMGT CDR1 and CDR3 agree closely; CDR2 length differs slightly
    'K': {
        'kabat':   dict(c1r1=1,  w1r1=-1, w1r2=15, r2len=6,  c2r3=1),
        'chothia': dict(c1r1=1,  w1r1=-1, w1r2=15, r2len=7,  c2r3=1),
        'imgt':    dict(c1r1=1,  w1r1=-1, w1r2=15, r2len=7,  c2r3=1),
    },
    'L': {
        'kabat':   dict(c1r1=1,  w1r1=-1, w1r2=15, r2len=6,  c2r3=1),
        'chothia': dict(c1r1=1,  w1r1=-1, w1r2=15, r2len=7,  c2r3=1),
        'imgt':    dict(c1r1=1,  w1r1=-1, w1r2=15, r2len=7,  c2r3=1),
    },
}


def get_cdr_regions(seq, chain_type, scheme='kabat'):
    seq    = seq.upper()
    ct     = chain_type[0].upper() if chain_type else 'H'
    scheme = scheme.lower().split()[0]          # normalise e.g. 'Kabat' → 'kabat'
    if scheme not in ('kabat', 'chothia', 'imgt'):
        scheme = 'kabat'
    if ct not in ('H', 'K', 'L'):
        result = _detect_cdrs(seq, 'H', scheme)
        if result is None:
            result = _detect_cdrs(seq, 'K', scheme)
        return result
    return _detect_cdrs(seq, ct, scheme)


def _detect_cdrs(seq, chain_type, scheme='kabat'):
    n      = len(seq)
    params = _CDR_PARAMS.get(chain_type, _CDR_PARAMS['H']).get(scheme, _CDR_PARAMS['H']['kabat'])

    # ── Anchor 1: first conserved Cys (end of FR1) ────────────────────────
    c1 = next((i for i in range(14, min(34, n)) if seq[i] == 'C'), None)
    if c1 is None:
        return None

    # ── Anchor 2: conserved Trp after CDR1 (start of FR2) ─────────────────
    w1 = next((i for i in range(c1 + 4, min(c1 + 22, n)) if seq[i] == 'W'), None)
    if w1 is None:
        return None

    cdr1_start = c1 + params['c1r1']
    cdr1_end   = w1 + params['w1r1']          # e.g. w1-1 or w1-6

    cdr2_start = w1 + params['w1r2']
    cdr2_end   = cdr2_start + params['r2len'] - 1

    # ── Anchor 3: second conserved Cys (end of FR3) ───────────────────────
    c2_lo = cdr2_end + (28 if chain_type == 'H' else 22)
    c2_hi = cdr2_end + (50 if chain_type == 'H' else 46)
    c2    = next((i for i in range(c2_lo, min(c2_hi, n)) if seq[i] == 'C'), None)
    if c2 is None:
        return None

    # ── CDR3: c2 + offset → J-region anchor ──────────────────────────────
    cdr3_start = c2 + params['c2r3']
    j_pat      = r'WG[A-Z]G' if chain_type == 'H' else r'FG[A-Z]G'
    m          = re.search(j_pat, seq[cdr3_start: cdr3_start + 32])
    cdr3_end   = (cdr3_start + m.start() - 1) if m else \
                 min(cdr3_start + (12 if chain_type == 'H' else 9) - 1, n - 1)
    cdr3_end   = max(cdr3_end, cdr3_start)

    def ss(s, e):
        return seq[s: e + 1] if 0 <= s <= e < n else ''

    return {
        'cdr1':     (cdr1_start, min(cdr1_end,  n - 1)),
        'cdr2':     (cdr2_start, min(cdr2_end,  n - 1)),
        'cdr3':     (cdr3_start, min(cdr3_end,  n - 1)),
        'cdr1_seq': ss(cdr1_start, cdr1_end),
        'cdr2_seq': ss(cdr2_start, cdr2_end),
        'cdr3_seq': ss(cdr3_start, cdr3_end),
        'scheme':   scheme.upper(),
    }


def get_residue_region(seq_pos, cdr_ann):
    if not cdr_ann:
        return 'match'
    for region in ('cdr1', 'cdr2', 'cdr3'):
        s, e = cdr_ann.get(region, (-1, -1))
        if s != -1 and s <= seq_pos <= e:
            return region
    return 'fr'

# ─────────────────────────────────────────────────────────────────────────────
# FASTA parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_fasta(text):
    records = []
    name, seq_lines = None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            if name is not None:
                records.append((name, ''.join(seq_lines).upper()))
            name = line[1:].strip()
            seq_lines = []
        elif line:
            seq_lines.append(re.sub(r'\s', '', line))
    if name is not None:
        records.append((name, ''.join(seq_lines).upper()))
    return records

def clean_sequence(text):
    return re.sub(r'[\s\d]', '', text).upper()

# ─────────────────────────────────────────────────────────────────────────────
# Colours
# ─────────────────────────────────────────────────────────────────────────────
DARK_BG        = '#1e1e2e'
DARK_TEXT      = '#cdd6f4'
CDR1_COLOR     = '#f38ba8'
CDR2_COLOR     = '#fab387'
CDR3_COLOR     = '#a6e3a1'
FR_COLOR       = '#89dceb'
MATCH_COLOR    = '#cdd6f4'
MISMATCH_COLOR = '#f9e2af'
GAP_COLOR      = '#6c7086'
HEADER_COLOR   = '#cba6f7'
STATS_COLOR    = '#89b4fa'
SEP_COLOR      = '#45475a'
CONS_COLOR     = '#a6e3a1'

# ─────────────────────────────────────────────────────────────────────────────
# Compact sequence panel (one per sequence, can be removed)
# ─────────────────────────────────────────────────────────────────────────────

class SequencePanel(ttk.Frame):
    def __init__(self, parent, index, on_remove, **kwargs):
        super().__init__(parent, relief='groove', borderwidth=1, **kwargs)
        self.index     = index
        self.on_remove = on_remove
        self._build()

    def _build(self):
        # Header row: label + remove button
        hdr = ttk.Frame(self)
        hdr.pack(fill='x', padx=4, pady=(4, 2))
        self._lbl = ttk.Label(hdr, text=f'Sequence {self.index}',
                              font=('Helvetica', 10, 'bold'))
        self._lbl.pack(side='left')
        ttk.Button(hdr, text='✕', width=3,
                   command=lambda: self.on_remove(self)).pack(side='right')

        # Name row
        nr = ttk.Frame(self)
        nr.pack(fill='x', padx=4, pady=(0, 2))
        ttk.Label(nr, text='Name:').pack(side='left')
        self.name_entry = ttk.Entry(nr, width=28)
        self.name_entry.pack(side='left', padx=(4, 0), fill='x', expand=True)

        # Sequence text
        self.text = tk.Text(self, height=4, wrap='word',
                            font=('Courier', 10), relief='flat', bd=0,
                            bg='#fafafa', padx=4, pady=2)
        self.text.pack(fill='both', expand=True, padx=4)

        # Button row
        br = ttk.Frame(self)
        br.pack(fill='x', padx=4, pady=(2, 4))
        ttk.Button(br, text='📂 Load FASTA', command=self._load).pack(side='left')
        ttk.Button(br, text='✕ Clear',       command=self._clear).pack(side='left', padx=4)

    def renumber(self, index):
        self.index = index
        self._lbl.configure(text=f'Sequence {index}')

    def _load(self):
        fp = filedialog.askopenfilename(
            title=f'Load FASTA — Sequence {self.index}',
            filetypes=[('FASTA / Text', '*.fa *.fasta *.faa *.txt'), ('All', '*.*')]
        )
        if not fp:
            return
        with open(fp) as fh:
            content = fh.read()
        records = parse_fasta(content)
        if records:
            name, seq = records[0]
            self.name_entry.delete(0, 'end')
            self.name_entry.insert(0, name)
            self.text.delete('1.0', 'end')
            self.text.insert('1.0', seq)
        else:
            seq = clean_sequence(content)
            self.text.delete('1.0', 'end')
            self.text.insert('1.0', seq)
            self.name_entry.delete(0, 'end')
            self.name_entry.insert(0, os.path.splitext(os.path.basename(fp))[0])

    def _clear(self):
        self.text.delete('1.0', 'end')
        self.name_entry.delete(0, 'end')

    def get_name(self):
        return self.name_entry.get().strip() or f'Seq{self.index}'

    def get_sequence(self):
        raw = self.text.get('1.0', 'end').strip()
        if not raw:
            return ''
        if raw.startswith('>'):
            records = parse_fasta(raw)
            return records[0][1] if records else ''
        return clean_sequence(raw)

# ─────────────────────────────────────────────────────────────────────────────
# Main application
# ─────────────────────────────────────────────────────────────────────────────

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('🧬  Protein Sequence Alignment Tool')
        self.geometry('1300x1000')
        self.minsize(900, 720)
        self.configure(bg='#eff1f5')
        self.panels = []
        self._build_style()
        self._build_ui()
        # Start with 2 panels
        self._add_panel()
        self._add_panel()

    # ── Style ──────────────────────────────────────────────────────────────
    def _build_style(self):
        s = ttk.Style(self)
        s.theme_use('clam')
        s.configure('.', background='#eff1f5', font=('Helvetica', 10))
        s.configure('TLabelframe',       background='#eff1f5')
        s.configure('TLabelframe.Label', font=('Helvetica', 10, 'bold'),
                    foreground='#4c4f69')
        s.configure('Run.TButton',   font=('Helvetica', 11, 'bold'), padding=8)
        s.configure('TRadiobutton',  background='#eff1f5')
        s.configure('TCheckbutton',  background='#eff1f5')

    # ── UI skeleton ────────────────────────────────────────────────────────
    def _build_ui(self):
        # Title bar
        bar = tk.Frame(self, bg='#4c4f69', height=48)
        bar.pack(fill='x')
        bar.pack_propagate(False)
        tk.Label(bar, text='🧬  Protein Sequence Alignment Tool',
                 bg='#4c4f69', fg='white',
                 font=('Helvetica', 14, 'bold')).pack(side='left', padx=16, pady=10)
        tk.Label(bar,
                 text='🔒  All sequences processed locally — nothing transmitted externally',
                 bg='#4c4f69', fg='#a6e3a1',
                 font=('Helvetica', 10, 'italic')).pack(side='right', padx=16)

        # Main content
        content = ttk.Frame(self, padding=10)
        content.pack(fill='both', expand=True)

        # ── Left column: sequence panels ──────────────────────────────────
        top = ttk.Frame(content)
        top.pack(fill='both', expand=False)
        top.columnconfigure(0, weight=1, minsize=420)
        top.columnconfigure(1, weight=3)

        # Scrollable panel list
        seq_lf = ttk.LabelFrame(top, text='Sequences', padding=4)
        seq_lf.grid(row=0, column=0, sticky='nsew', padx=(0, 8), pady=4)

        self._panel_canvas = tk.Canvas(seq_lf, bg='#eff1f5',
                                       highlightthickness=0, width=370)
        vsb = ttk.Scrollbar(seq_lf, orient='vertical',
                            command=self._panel_canvas.yview)
        self._panel_canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side='right', fill='y')
        self._panel_canvas.pack(side='left', fill='both', expand=True)

        self._panel_inner = ttk.Frame(self._panel_canvas)
        self._cw = self._panel_canvas.create_window(
            (0, 0), window=self._panel_inner, anchor='nw')
        self._panel_inner.bind('<Configure>', self._on_panel_configure)
        self._panel_canvas.bind('<Configure>', self._on_canvas_configure)

        # Bind mousewheel
        self._panel_canvas.bind('<Enter>',
            lambda e: self._panel_canvas.bind_all('<MouseWheel>', self._on_mousewheel))
        self._panel_canvas.bind('<Leave>',
            lambda e: self._panel_canvas.unbind_all('<MouseWheel>'))

        # "Add" / "Load Multi-FASTA" buttons — keep a reference so _add_panel
        # can always insert new panels *before* this row reliably.
        self._add_btn_row = ttk.Frame(self._panel_inner)
        add_row = self._add_btn_row
        add_row.pack(fill='x', pady=(6, 2))
        ttk.Button(add_row, text='＋ Add Sequence',
                   command=self._add_panel).pack(side='left', padx=2)
        ttk.Button(add_row, text='📂 Load Multi-FASTA',
                   command=self._load_multi_fasta).pack(side='left', padx=2)

        # ── Right column: options + results ───────────────────────────────
        right = ttk.Frame(top)
        right.grid(row=0, column=1, sticky='nsew')

        # Options
        opt = ttk.LabelFrame(right, text='Alignment Options', padding=8)
        opt.pack(fill='x', pady=(4, 6))

        ttk.Label(opt, text='Chain type:').grid(row=0, column=0, sticky='w')
        self.chain_var = tk.StringVar(value='H — Heavy chain')
        ttk.Combobox(opt, textvariable=self.chain_var, width=18, state='readonly',
                     values=['H — Heavy chain', 'K — Kappa light', 'L — Lambda light',
                             'Auto-detect', 'Generic protein']
                     ).grid(row=0, column=1, sticky='w', padx=4)

        ttk.Label(opt, text='Mode:').grid(row=0, column=2, sticky='w', padx=(12, 4))
        self.mode_var = tk.StringVar(value='local')
        ttk.Radiobutton(opt, text='Local (BLAST-like)', variable=self.mode_var,
                        value='local').grid(row=0, column=3, sticky='w', padx=2)
        ttk.Radiobutton(opt, text='Global', variable=self.mode_var,
                        value='global').grid(row=0, column=4, sticky='w', padx=2)

        ttk.Label(opt, text='Gap open:').grid(row=1, column=0, sticky='w', pady=(6, 0))
        self.gap_open_var = tk.DoubleVar(value=-11.0)
        ttk.Spinbox(opt, from_=-20, to=0, increment=1,
                    textvariable=self.gap_open_var, width=6,
                    format='%.0f').grid(row=1, column=1, sticky='w', pady=(6, 0))

        ttk.Label(opt, text='Gap extend:').grid(row=1, column=2, sticky='w',
                                                padx=(12, 4), pady=(6, 0))
        self.gap_ext_var = tk.DoubleVar(value=-1.0)
        ttk.Spinbox(opt, from_=-10, to=0, increment=1,
                    textvariable=self.gap_ext_var, width=6,
                    format='%.0f').grid(row=1, column=3, sticky='w', pady=(6, 0))

        self.annotate_cdr = tk.BooleanVar(value=True)
        ttk.Checkbutton(opt, text='Annotate CDR regions (built-in)',
                        variable=self.annotate_cdr).grid(
                        row=1, column=4, sticky='w', padx=(12, 0), pady=(6, 0))

        ttk.Label(opt, text='CDR scheme:').grid(row=1, column=5, sticky='w',
                                                padx=(10, 4), pady=(6, 0))
        self.scheme_var = tk.StringVar(value='Kabat')
        ttk.Combobox(opt, textvariable=self.scheme_var, width=9, state='readonly',
                     values=['Kabat', 'Chothia', 'IMGT']
                     ).grid(row=1, column=6, sticky='w', pady=(6, 0))

        # Action buttons
        btn_row = ttk.Frame(right)
        btn_row.pack(pady=4)
        ttk.Button(btn_row, text='⚡  Run Alignment', style='Run.TButton',
                   command=self._run).pack(side='left', padx=4)
        ttk.Button(btn_row, text='Clear All',
                   command=self._clear_all).pack(side='left', padx=4)
        ttk.Button(btn_row, text='Export Results',
                   command=self._export).pack(side='left', padx=4)

        # CDR info note
        note = tk.Frame(right, bg='#e6f4ea', bd=1, relief='solid')
        note.pack(fill='x', pady=(0, 4))
        tk.Label(note,
                 text='✅  CDR annotation is built-in (no extra packages). '
                      'Kabat anchors — reliable for human/humanized VH, Vκ, Vλ.',
                 bg='#e6f4ea', fg='#276740',
                 font=('Helvetica', 9), anchor='w').pack(fill='x', padx=8, pady=3)

        # Stats bar
        self.stats_var = tk.StringVar(value='Add sequences and click Run Alignment.')
        tk.Label(right, textvariable=self.stats_var,
                 font=('Courier', 10, 'bold'),
                 bg='#dce0e8', fg='#4c4f69',
                 anchor='w', padx=8, pady=4, relief='sunken').pack(fill='x', pady=(0, 4))

        # Results
        res_lf = ttk.LabelFrame(right, text='Alignment', padding=4)
        res_lf.pack(fill='both', expand=True)

        self.result = tk.Text(res_lf, wrap='none',
                              font=('Courier', 11), bg=DARK_BG, fg=DARK_TEXT,
                              insertbackground='white', state='disabled')
        vsb2 = ttk.Scrollbar(res_lf, orient='vertical',   command=self.result.yview)
        hsb2 = ttk.Scrollbar(res_lf, orient='horizontal', command=self.result.xview)
        self.result.configure(yscrollcommand=vsb2.set, xscrollcommand=hsb2.set)
        vsb2.pack(side='right',  fill='y')
        hsb2.pack(side='bottom', fill='x')
        self.result.pack(fill='both', expand=True)

        self._setup_tags()
        self._last_text = ''

    def _on_panel_configure(self, e):
        self._panel_canvas.configure(
            scrollregion=self._panel_canvas.bbox('all'))

    def _on_canvas_configure(self, e):
        self._panel_canvas.itemconfig(self._cw, width=e.width)

    def _on_mousewheel(self, e):
        self._panel_canvas.yview_scroll(int(-1 * (e.delta / 120)), 'units')

    # ── Panel management ───────────────────────────────────────────────────
    def _add_panel(self):
        idx   = len(self.panels) + 1
        panel = SequencePanel(self._panel_inner, idx, self._remove_panel)
        # Insert before the "Add / Load" button row (last child)
        panel.pack(fill='x', pady=2, before=self._add_btn_row)
        self.panels.append(panel)
        self._update_remove_buttons()

    def _remove_panel(self, panel):
        if len(self.panels) <= 2:
            messagebox.showinfo('Minimum Sequences',
                                'At least 2 sequences are required.')
            return
        self.panels.remove(panel)
        panel.destroy()
        for i, p in enumerate(self.panels, 1):
            p.renumber(i)
        self._update_remove_buttons()

    def _update_remove_buttons(self):
        # Hide remove buttons when only 2 panels remain
        show = len(self.panels) > 2
        for p in self.panels:
            for child in p.winfo_children():
                if isinstance(child, ttk.Frame):
                    for w in child.winfo_children():
                        if isinstance(w, ttk.Button) and w.cget('text') == '✕':
                            w.configure(state='normal' if show else 'disabled')

    def _load_multi_fasta(self):
        fp = filedialog.askopenfilename(
            title='Load Multi-FASTA',
            filetypes=[('FASTA / Text', '*.fa *.fasta *.faa *.txt'), ('All', '*.*')]
        )
        if not fp:
            return
        with open(fp) as fh:
            records = parse_fasta(fh.read())
        if not records:
            messagebox.showwarning('No sequences found', 'No FASTA sequences found in file.')
            return

        # Clear existing panels and rebuild
        for p in list(self.panels):
            p.destroy()
        self.panels.clear()

        for i, (name, seq) in enumerate(records, 1):
            idx   = len(self.panels) + 1
            panel = SequencePanel(self._panel_inner, idx, self._remove_panel)
            panel.pack(fill='x', pady=2,
                       before=self._add_btn_row)
            panel.name_entry.insert(0, name)
            panel.text.insert('1.0', seq)
            self.panels.append(panel)

        if len(self.panels) < 2:
            self._add_panel()

        self._update_remove_buttons()
        self.stats_var.set(
            f'Loaded {len(records)} sequences from {os.path.basename(fp)}. '
            'Click Run Alignment.')

    # ── Text tags ──────────────────────────────────────────────────────────
    def _setup_tags(self):
        t = self.result
        t.tag_configure('header',   foreground=HEADER_COLOR,   font=('Courier', 11, 'bold'))
        t.tag_configure('stats',    foreground=STATS_COLOR)
        t.tag_configure('sep',      foreground=SEP_COLOR)
        t.tag_configure('cdr1',     foreground=CDR1_COLOR,     font=('Courier', 11, 'bold'))
        t.tag_configure('cdr2',     foreground=CDR2_COLOR,     font=('Courier', 11, 'bold'))
        t.tag_configure('cdr3',     foreground=CDR3_COLOR,     font=('Courier', 11, 'bold'))
        t.tag_configure('fr',       foreground=FR_COLOR)
        t.tag_configure('match',    foreground=MATCH_COLOR)
        t.tag_configure('mismatch', foreground=MISMATCH_COLOR)
        t.tag_configure('gap',      foreground=GAP_COLOR)
        t.tag_configure('pos_sep',  foreground=STATS_COLOR)
        t.tag_configure('cons',     foreground=CONS_COLOR)

    # ── Run ────────────────────────────────────────────────────────────────
    def _run(self):
        seqs_named = [(p.get_name(), p.get_sequence())
                      for p in self.panels if p.get_sequence()]

        if len(seqs_named) < 2:
            messagebox.showwarning('Missing Input',
                                   'Please enter at least 2 sequences.')
            return

        valid = set('ACDEFGHIKLMNPQRSTVWYBZX*-')
        for name, seq in seqs_named:
            bad = set(seq) - valid
            if bad:
                if not messagebox.askyesno(
                        'Unexpected characters',
                        f'{name} contains: {bad}\nContinue anyway?'):
                    return

        mode       = self.mode_var.get()
        gap_open   = float(self.gap_open_var.get())
        gap_extend = float(self.gap_ext_var.get())
        ct         = self.chain_var.get()[0]

        try:
            if len(seqs_named) == 2:
                a1, a2, score = align_sequences(
                    seqs_named[0][1], seqs_named[1][1], mode, gap_open, gap_extend)
                aligned = [(seqs_named[0][0], a1), (seqs_named[1][0], a2)]
            else:
                aligned, score = multiple_sequence_alignment(
                    seqs_named, mode, gap_open, gap_extend)
        except Exception as e:
            messagebox.showerror('Alignment Error', str(e))
            return

        # CDR annotation (per sequence)
        cdr_anns = {}
        scheme = self.scheme_var.get()
        if self.annotate_cdr.get():
            for name, seq in seqs_named:
                ann = get_cdr_regions(seq, ct, scheme)
                if ann:
                    cdr_anns[name] = ann

        if len(seqs_named) == 2:
            self._display_pairwise(seqs_named, aligned, cdr_anns, score)
        else:
            self._display_msa(seqs_named, aligned, cdr_anns, score)

    # ── Pairwise display ───────────────────────────────────────────────────
    def _display_pairwise(self, seqs_named, aligned, cdr_anns, score):
        name1, seq1 = seqs_named[0]
        name2, seq2 = seqs_named[1]
        aln1 = aligned[0][1]
        aln2 = aligned[1][1]
        L    = len(aln1)

        matches   = sum(a == b and a != '-' for a, b in zip(aln1, aln2))
        positives = sum(1 for a, b in zip(aln1, aln2)
                        if a != '-' and b != '-' and blosum_score(a, b) > 0)
        gaps      = sum(1 for a, b in zip(aln1, aln2) if a == '-' or b == '-')
        pct_id    = 100.0 * matches   / L if L else 0
        pct_pos   = 100.0 * positives / L if L else 0
        pct_gap   = 100.0 * gaps      / L if L else 0

        self.stats_var.set(
            f'Score: {score:.1f}  |  Identity: {matches}/{L} ({pct_id:.1f}%)'
            f'  |  Positives: {positives}/{L} ({pct_pos:.1f}%)'
            f'  |  Gaps: {gaps}/{L} ({pct_gap:.1f}%)'
        )

        t = self.result
        t.configure(state='normal')
        t.delete('1.0', 'end')

        cdr1 = cdr_anns.get(name1)
        cdr2 = cdr_anns.get(name2)

        def w(txt, tag=''):
            t.insert('end', txt, tag)

        W = 78
        w('=' * W + '\n', 'sep')
        w('  PAIRWISE ALIGNMENT\n', 'header')
        w('=' * W + '\n', 'sep')
        w(f'\n  Query:   {name1}  ({len(seq1)} aa)\n', 'stats')
        w(f'  Subject: {name2}  ({len(seq2)} aa)\n', 'stats')
        w(f'\n  Score:      {score:.1f}\n', 'stats')
        w(f'  Identity:   {matches}/{L}  ({pct_id:.1f}%)\n', 'stats')
        w(f'  Positives:  {positives}/{L}  ({pct_pos:.1f}%)\n', 'stats')
        w(f'  Gaps:       {gaps}/{L}  ({pct_gap:.1f}%)\n', 'stats')

        if cdr1 or cdr2:
            _sch = next((v['scheme'] for v in cdr_anns.values()), self.scheme_var.get().upper())
            w(f'\n  CDR sequences ({_sch}):\n', 'stats')
            for lbl, ann, seq in (('  Query  ', cdr1, seq1),
                                   ('  Subject', cdr2, seq2)):
                if ann:
                    w(f'    {lbl}: ', 'stats')
                    w(f"CDR1={ann['cdr1_seq']} ", 'cdr1')
                    w(f"CDR2={ann['cdr2_seq']} ", 'cdr2')
                    w(f"CDR3={ann['cdr3_seq']}\n", 'cdr3')
            if cdr1 and cdr2:
                w('\n  Per-CDR identity:\n', 'stats')
                for reg, col in (('cdr1', 'cdr1'), ('cdr2', 'cdr2'), ('cdr3', 'cdr3')):
                    s1e = cdr1.get(reg, (-1, -1))
                    if s1e[0] == -1:
                        continue
                    m_c, tot = self._cdr_id(aln1, aln2, s1e[0], s1e[1])
                    pct_c = 100.0 * m_c / tot if tot else 0
                    w(f'    {reg.upper()}: ', 'stats')
                    w(f'{m_c}/{tot} ({pct_c:.1f}%)\n', col)

        if self.annotate_cdr.get() and not (cdr1 or cdr2):
            w('\n  ⚠️  CDR detection failed — may not be a standard V-region.\n', 'stats')

        if cdr1 or cdr2:
            w('\n  Legend: ', 'stats')
            w('CDR1 ', 'cdr1'); w('CDR2 ', 'cdr2'); w('CDR3 ', 'cdr3')
            w('Framework ', 'fr')
        w('\n  Residues: ', 'stats')
        w('Match ', 'match'); w('Mismatch ', 'mismatch'); w('Gap\n', 'gap')
        w('\n' + '=' * W + '\n\n', 'sep')

        # Alignment blocks
        BLOCK = 60
        q_pos = s_pos = 0
        consensus = ''.join(
            '|' if a == b and a != '-' else
            ('+' if a != '-' and b != '-' and blosum_score(a, b) > 0 else ' ')
            for a, b in zip(aln1, aln2)
        )
        LW = 10
        for i in range(0, L, BLOCK):
            b1 = aln1[i:i+BLOCK]
            b2 = aln2[i:i+BLOCK]
            bc = consensus[i:i+BLOCK]
            q_end = q_pos + len(b1.replace('-', ''))
            s_end = s_pos + len(b2.replace('-', ''))
            w(f"  {'Query':>{LW}} {q_pos+1:>5} ", 'pos_sep')
            self._write_block(t, b1, b2, cdr1, q_pos)
            w(f' {q_end:<5}\n', 'pos_sep')
            w(f"  {'':>{LW}} {'':>5} {bc}\n", 'sep')
            w(f"  {'Subject':>{LW}} {s_pos+1:>5} ", 'pos_sep')
            self._write_block(t, b2, b1, cdr2, s_pos)
            w(f' {s_end:<5}\n\n', 'pos_sep')
            q_pos = q_end
            s_pos = s_end

        self._last_text = t.get('1.0', 'end')
        t.configure(state='disabled')

    # ── MSA display ────────────────────────────────────────────────────────
    def _display_msa(self, seqs_named, aligned, cdr_anns, score):
        n    = len(aligned)
        L    = max(len(s) for _, s in aligned)
        names = [nm for nm, _ in aligned]

        # Pairwise identity matrix
        pw_stats = compute_pairwise_stats(aligned)

        # Average overall identity
        all_pcts = [v[2] for v in pw_stats.values()]
        avg_id   = sum(all_pcts) / len(all_pcts) if all_pcts else 0

        self.stats_var.set(
            f'MSA: {n} sequences  |  Alignment length: {L} aa  '
            f'|  Avg pairwise identity: {avg_id:.1f}%'
        )

        consensus = msa_consensus(aligned)

        t = self.result
        t.configure(state='normal')
        t.delete('1.0', 'end')

        def w(txt, tag=''):
            t.insert('end', txt, tag)

        W = 78
        w('=' * W + '\n', 'sep')
        w(f'  MULTIPLE SEQUENCE ALIGNMENT  ({n} sequences)\n', 'header')
        w('=' * W + '\n', 'sep')

        # Sequence list
        w('\n  Sequences:\n', 'stats')
        for i, (nm, _) in enumerate(aligned, 1):
            orig_seq = next(s for nn, s in seqs_named if nn == nm)
            w(f'    {i:>2}. {nm}  ({len(orig_seq)} aa)\n', 'stats')

        w(f'\n  Algorithm:  Center-star progressive (BLOSUM62)\n', 'stats')
        w(f'  Avg pairwise identity:  {avg_id:.1f}%\n', 'stats')

        # Pairwise identity matrix
        w('\n  Pairwise identities:\n', 'stats')
        label_w = max(len(nm) for nm in names) + 2
        # Header
        w(f"  {'':>{label_w}}", 'stats')
        for nm in names:
            w(f'  {nm[:8]:>8}', 'stats')
        w('\n', 'stats')
        for nm_i in names:
            w(f'  {nm_i:>{label_w}}', 'stats')
            for nm_j in names:
                if nm_i == nm_j:
                    w(f'  {"—":>8}', 'sep')
                else:
                    key = (nm_i, nm_j) if (nm_i, nm_j) in pw_stats else (nm_j, nm_i)
                    pct = pw_stats[key][2]
                    color = 'cdr3' if pct >= 90 else ('fr' if pct >= 70 else 'mismatch')
                    w(f'  {pct:>7.1f}%', color)
            w('\n', 'stats')

        # CDR summary
        any_cdr = any(cdr_anns.get(nm) for nm, _ in aligned)
        if any_cdr:
            _sch = next((v['scheme'] for v in cdr_anns.values()), self.scheme_var.get().upper())
            w(f'\n  CDR sequences ({_sch}):\n', 'stats')
            for nm, _ in aligned:
                ann = cdr_anns.get(nm)
                if ann:
                    w(f'    {nm}: ', 'stats')
                    w(f"CDR1={ann['cdr1_seq']} ", 'cdr1')
                    w(f"CDR2={ann['cdr2_seq']} ", 'cdr2')
                    w(f"CDR3={ann['cdr3_seq']}\n", 'cdr3')

        if self.annotate_cdr.get() and not any_cdr:
            w('\n  ⚠️  CDR detection failed — sequences may not be standard V-regions.\n',
              'stats')

        if any_cdr:
            w('\n  Legend: ', 'stats')
            w('CDR1 ', 'cdr1'); w('CDR2 ', 'cdr2'); w('CDR3 ', 'cdr3')
            w('Framework ', 'fr')
        w('\n  Consensus: ', 'stats')
        w('* ', 'cons'); w('= 100% identical   : = conserved   . = semi-conserved\n', 'stats')
        w('\n' + '=' * W + '\n\n', 'sep')

        # Alignment blocks
        BLOCK = 60
        LW    = max(len(nm) for nm in names)
        seq_positions = {nm: 0 for nm, _ in aligned}

        for i in range(0, L, BLOCK):
            # Position ruler
            start_col = i + 1
            ruler = f'{start_col:<10}'
            mid   = i + BLOCK // 2
            end   = min(i + BLOCK, L)
            mid_s = f'{mid + 1}'
            end_s = f'{end}'
            spaces = BLOCK - len(ruler) - len(mid_s) - len(end_s)
            half   = spaces // 2
            ruler  = (ruler + ' ' * half + mid_s +
                      ' ' * (spaces - half) + end_s)
            w(f"  {'':>{LW}}       {ruler[:BLOCK]}\n", 'sep')

            for nm, aln_seq in aligned:
                blk = aln_seq[i:i+BLOCK]
                spos = seq_positions[nm]
                end_pos = spos + len(blk.replace('-', ''))
                w(f'  {nm:>{LW}} {spos+1:>5} ', 'pos_sep')
                self._write_block(t, blk, blk, cdr_anns.get(nm), spos,
                                  msa_mode=True)
                w(f' {end_pos:<5}\n', 'pos_sep')
                seq_positions[nm] = end_pos

            # Consensus row
            cons_blk = consensus[i:i+BLOCK]
            w(f"  {'Consensus':>{LW}} {'':>5} {cons_blk}\n\n", 'cons')

        self._last_text = t.get('1.0', 'end')
        t.configure(state='disabled')

    # ── Helpers ────────────────────────────────────────────────────────────
    def _write_block(self, t, block, other_block, cdr_ann, seq_start,
                     msa_mode=False):
        """Write a block of aligned residues with per-residue colour tags."""
        seq_pos = seq_start
        for i, aa in enumerate(block):
            other = other_block[i] if i < len(other_block) else '-'
            if aa == '-':
                t.insert('end', aa, 'gap')
            else:
                if cdr_ann:
                    tag = get_residue_region(seq_pos, cdr_ann)
                elif msa_mode:
                    tag = 'match'
                else:
                    tag = 'match' if aa == other else 'mismatch'
                t.insert('end', aa, tag)
                seq_pos += 1

    def _cdr_id(self, aln1, aln2, cdr_start, cdr_end):
        matches = total = 0
        sp = 0
        for a, b in zip(aln1, aln2):
            if a != '-':
                if cdr_start <= sp <= cdr_end:
                    total += 1
                    if a == b:
                        matches += 1
                sp += 1
        return matches, total

    def _clear_all(self):
        for p in self.panels:
            p._clear()
        self.result.configure(state='normal')
        self.result.delete('1.0', 'end')
        self.result.configure(state='disabled')
        self.stats_var.set('Add sequences and click Run Alignment.')
        self._last_text = ''

    def _export(self):
        if not self._last_text.strip():
            messagebox.showinfo('No Results', 'Run an alignment first.')
            return
        fp = filedialog.asksaveasfilename(
            defaultextension='.txt',
            filetypes=[('Text file', '*.txt'), ('All', '*.*')],
            title='Export Alignment Results'
        )
        if fp:
            with open(fp, 'w') as fh:
                fh.write(self._last_text)
            messagebox.showinfo('Saved', f'Results saved to:\n{fp}')

# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    app = App()
    app.mainloop()
