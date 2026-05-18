#!/usr/bin/env python3
"""
Protein Liability Analyzer v2 - Desktop GUI  (Auto-HOS Edition)
================================================================
HOS is ALWAYS performed:
  - No PDB supplied  -> Chou-Fasman + sequence-based RSA estimation
  - PDB supplied     -> HELIX/SHEET records + Shrake-Rupley SASA

Three tabs:
  1.  Analysis Setup
  2.  Results & Log   (streaming console output)
  3.  Report          (interactive findings table, auto-populated after each run)

Run:  python protein_liability_gui_autohos.py
Requires protein_liability_analyzer_v2.py in the same folder.
All computation is local - sequences never transmitted.
"""

import re
import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
import sys
import os
import csv
import queue
import webbrowser
import platform
import subprocess
from pathlib import Path
from datetime import datetime

# ── Locate the analysis module ────────────────────────────────────────────────
# PyInstaller sets sys.frozen = True inside a bundled executable.
# Use sys.executable (the .exe / .app binary) so _HERE points to the folder
# the user actually installed into, not the temp extraction directory.
if getattr(sys, "frozen", False):
    _HERE = Path(sys.executable).resolve().parent
else:
    _HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

try:
    import protein_liability_analyzer_v2 as pla
    CORE_OK  = True
    CORE_ERR = ""
except ImportError as _e:
    CORE_OK  = False
    CORE_ERR = str(_e)

# ── Palette ───────────────────────────────────────────────────────────────────
DARK    = "#1a1a2e"
ACCENT  = "#0f3460"
ACCENT2 = "#2980b9"
LIGHT   = "#f0f2f5"
WHITE   = "#ffffff"
BORDER  = "#dde1e7"
GREY    = "#6c757d"
LOG_BG  = "#0d1117"
LOG_FG  = "#c9d1d9"

# Risk palette
R_HIGH_BG   = "#ffeaea";  R_HIGH_FG   = "#8b0000"
R_MED_BG    = "#fff3e0";  R_MED_FG    = "#7d4c00"
R_LOW_BG    = "#e8f4fd";  R_LOW_FG    = "#154360"
R_INFO_BG   = "#f5f5f5";  R_INFO_FG   = "#444444"

F_BODY  = ("Segoe UI", 10)
F_SMALL = ("Segoe UI",  9)
F_BOLD  = ("Segoe UI", 10, "bold")
F_TITLE = ("Segoe UI", 16, "bold")
F_MONO  = ("Courier New", 9)


# ── Thread-safe stdout redirector ─────────────────────────────────────────────
class _QWriter:
    def __init__(self, q: queue.Queue, tag: str = "normal"):
        self._q, self._tag = q, tag
    def write(self, s: str):
        if s: self._q.put((self._tag, s))
    def flush(self): pass


# ── Risk rank helper ──────────────────────────────────────────────────────────
def _risk_rank(level: str) -> int:
    return {"high": 0, "medium": 1, "low": 2, "info": 3}.get(level, 4)


# ── Pure-tkinter 3D molecular viewer ─────────────────────────────────────────
class Mol3DCanvas(tk.Canvas):
    """
    Software-rendered 3D Cα-trace viewer — no external dependencies.
    Draws backbone as connected lines with coloured spheres (ovals) for
    liability residues. Supports drag-to-rotate and scroll-to-zoom.
    """

    _CHAIN_COLORS = [
        "#4a90d9", "#e67e22", "#2ecc71", "#9b59b6",
        "#e74c3c", "#1abc9c", "#f39c12", "#3498db",
    ]
    _RISK_COLORS = {
        "high":   "#e74c3c",
        "medium": "#f39c12",
        "low":    "#3498db",
        "info":   "#95a5a6",
    }
    _CDR_COLORS_3D = {
        "CDR-H1": "#6ab0f5", "CDR-H2": "#f5c842", "CDR-H3": "#f07070",
        "CDR-L1": "#6fcf97", "CDR-L2": "#f5c842", "CDR-L3": "#f07070",
    }
    _AA_COLORS = {
        "A": "#e8a87c", "V": "#e8a87c", "I": "#e8a87c", "L": "#e8a87c",
        "M": "#e8a87c", "F": "#d4956a", "W": "#d4956a", "P": "#e8a87c",
        "S": "#88d8a3", "T": "#88d8a3", "C": "#f5e642", "Y": "#88d8a3",
        "N": "#88d8a3", "Q": "#88d8a3",
        "D": "#e07070", "E": "#e07070", "K": "#70a0e0", "R": "#70a0e0",
        "H": "#9b80d0", "G": "#cccccc",
    }
    _SS_COLORS = {"H": "#e74c3c", "E": "#3498db"}

    def __init__(self, parent, **kwargs):
        kwargs.setdefault("bg", DARK)
        kwargs.setdefault("highlightthickness", 0)
        super().__init__(parent, **kwargs)
        self._atoms: list      = []
        self._bonds: list      = []
        self._rot              = [[1.0, 0.0, 0.0],
                                  [0.0, 1.0, 0.0],
                                  [0.0, 0.0, 1.0]]
        self._scale            = 4.0
        self._drag_start       = None
        self._color_mode       = "risk"
        self._show_labels      = True
        self._show_cdr         = False
        self._risk_filter: set = {"high", "medium", "low", "info"}

        self.bind("<ButtonPress-1>", self._on_press)
        self.bind("<B1-Motion>",     self._on_drag)
        self.bind("<MouseWheel>",    self._on_scroll_mac)
        self.bind("<Button-4>",      lambda _: self._zoom(1.1))
        self.bind("<Button-5>",      lambda _: self._zoom(0.9))
        self.bind("<Configure>",     lambda _: self._redraw())
        self._draw_placeholder()

    # ── Public API ────────────────────────────────────────────────────────────
    def load(self, atoms: list, bonds: list):
        self._atoms = atoms
        self._bonds = bonds
        if atoms:
            max_r = max(math.sqrt(sum(v * v for v in a["xyz"])) for a in atoms)
            w = max(self.winfo_width(), 200)
            h = max(self.winfo_height(), 200)
            self._scale = min(w, h) * 0.38 / max(max_r, 1.0)
        self._rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        self._redraw()

    def clear(self):
        self._atoms = []
        self._bonds = []
        self.delete("all")
        self._draw_placeholder()

    def set_color_mode(self, mode: str):
        self._color_mode = mode
        self._redraw()

    def set_show_labels(self, val: bool):
        self._show_labels = val
        self._redraw()

    def set_show_cdr(self, val: bool):
        self._show_cdr = val
        self._redraw()

    def set_risk_filter(self, risks: set):
        self._risk_filter = risks
        self._redraw()

    def reset_view(self):
        self._rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        if self._atoms:
            max_r = max(math.sqrt(sum(v * v for v in a["xyz"])) for a in self._atoms)
            w = max(self.winfo_width(), 200)
            h = max(self.winfo_height(), 200)
            self._scale = min(w, h) * 0.38 / max(max_r, 1.0)
        self._redraw()

    # ── Projection ────────────────────────────────────────────────────────────
    def _proj(self, xyz):
        x, y, z = xyz
        r = self._rot
        px = r[0][0]*x + r[0][1]*y + r[0][2]*z
        py = r[1][0]*x + r[1][1]*y + r[1][2]*z
        pz = r[2][0]*x + r[2][1]*y + r[2][2]*z
        w = max(self.winfo_width(),  1)
        h = max(self.winfo_height(), 1)
        return w/2 + px * self._scale, h/2 - py * self._scale, pz

    # ── Color resolution ──────────────────────────────────────────────────────
    def _atom_color(self, atom: dict) -> str:
        mode = self._color_mode
        if mode == "chain":
            ci = ord(atom["chain"]) % len(self._CHAIN_COLORS) if atom["chain"] else 0
            return self._CHAIN_COLORS[ci]
        if mode == "cdr":
            return atom.get("color_cdr") or "#4a6fa5"
        if mode == "aa":
            return self._AA_COLORS.get(atom.get("aa", "G"), "#aaaaaa")
        if mode == "ss":
            return self._SS_COLORS.get(atom.get("ss_char", ""), "#666688")
        # default: "risk"
        if atom.get("is_liability"):
            return self._RISK_COLORS.get(atom.get("hl_risk", "info"), "#95a5a6")
        return "#4a6fa5"

    # ── Rendering ─────────────────────────────────────────────────────────────
    def _redraw(self):
        self.delete("all")
        if not self._atoms:
            self._draw_placeholder()
            return

        projected = [self._proj(a["xyz"]) for a in self._atoms]

        # Bonds — sorted back-to-front for painter's algorithm
        bond_list = []
        for i, j in self._bonds:
            x1, y1, z1 = projected[i]
            x2, y2, z2 = projected[j]
            col_i = self._atom_color(self._atoms[i])
            col_j = self._atom_color(self._atoms[j])
            col = col_i if z1 >= z2 else col_j
            if col == "#4a6fa5":
                col = "#2a4070"
            bond_list.append(((z1 + z2) / 2, x1, y1, x2, y2, col))
        bond_list.sort(key=lambda t: t[0])
        for _, x1, y1, x2, y2, col in bond_list:
            self.create_line(x1, y1, x2, y2, fill=col, width=2, tags="bond")

        # Atoms — back-to-front
        order = sorted(range(len(self._atoms)), key=lambda k: projected[k][2])
        for k in order:
            a = self._atoms[k]
            sx, sy, _ = projected[k]
            is_hl = a.get("is_liability", False)
            risk   = a.get("hl_risk", "")

            if is_hl and risk not in self._risk_filter:
                # Suppressed by filter — draw tiny dim dot
                self.create_oval(sx-2, sy-2, sx+2, sy+2,
                                 fill="#334466", outline="", tags="atom")
                continue

            color = self._atom_color(a)

            if is_hl:
                r = float(a.get("hl_size", 8))
                self.create_oval(sx-r, sy-r, sx+r, sy+r,
                                 fill=color, outline="#ffffff", width=1.5,
                                 tags="atom_hl")
                if self._show_labels and a.get("hl_label"):
                    self.create_text(sx, sy - r - 5,
                                     text=a["hl_label"], fill="#ffffff",
                                     font=("Courier New", 8), tags="label")
            else:
                # CDR backbone highlight
                if self._show_cdr and a.get("cdr_name"):
                    cdr_col = self._CDR_COLORS_3D.get(a["cdr_name"], "#888888")
                    self.create_oval(sx-4, sy-4, sx+4, sy+4,
                                     fill=cdr_col, outline="", tags="atom_cdr")
                else:
                    self.create_oval(sx-2.5, sy-2.5, sx+2.5, sy+2.5,
                                     fill=color, outline="", tags="atom")

        # CDR legend strip
        if self._show_cdr:
            self._draw_cdr_legend()

    def _draw_cdr_legend(self):
        present = {}
        for a in self._atoms:
            cdr = a.get("cdr_name")
            if cdr:
                present[cdr] = self._CDR_COLORS_3D.get(cdr, "#888888")
        if not present:
            return
        h = max(self.winfo_height(), 100)
        x, y = 10, h - 22
        self.create_text(x, y, text="CDR:", fill="#8899bb",
                         font=("Segoe UI", 8), anchor="w", tags="legend")
        x += 38
        for name in sorted(present):
            col = present[name]
            self.create_rectangle(x, y-6, x+10, y+6, fill=col, outline="", tags="legend")
            self.create_text(x + 14, y, text=name, fill="#ccccdd",
                             font=("Segoe UI", 8), anchor="w", tags="legend")
            x += 62

    def _draw_placeholder(self):
        w = max(self.winfo_width(),  400)
        h = max(self.winfo_height(), 300)
        self.create_text(w // 2, h // 2,
                         text="Run an analysis with a PDB file\nto enable 3D structure visualization",
                         fill="#445588", font=("Segoe UI", 13), justify="center")

    # ── Interaction ───────────────────────────────────────────────────────────
    def _on_press(self, event):
        self._drag_start = (event.x, event.y)

    def _on_drag(self, event):
        if not self._drag_start:
            return
        dx = event.x - self._drag_start[0]
        dy = event.y - self._drag_start[1]
        self._drag_start = (event.x, event.y)
        self._rotate(dy * 0.007, dx * 0.007)

    def _rotate(self, ax: float, ay: float):
        cx, sx = math.cos(ax), math.sin(ax)
        cy, sy = math.cos(ay), math.sin(ay)
        Rx = [[1, 0, 0], [0, cx, -sx], [0, sx, cx]]
        Ry = [[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]]

        def mm(A, B):
            return [[sum(A[i][k] * B[k][j] for k in range(3))
                     for j in range(3)] for i in range(3)]

        self._rot = mm(Rx, mm(Ry, self._rot))
        self._redraw()

    def _zoom(self, factor: float):
        self._scale *= factor
        self._redraw()

    def _on_scroll_mac(self, event):
        self._zoom(1.1 if event.delta > 0 else 0.9)


# ══════════════════════════════════════════════════════════════════════════════
class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Protein Liability Analyzer v2 - Auto-HOS - Letarte Scientific")
        self.root.configure(bg=LIGHT)
        self.root.minsize(860, 680)
        self._centre(1040, 860)

        # State
        self.v_mode     = tk.StringVar(value="paste")
        self.v_fasta    = tk.StringVar()
        self.v_pdb      = tk.StringVar()
        self.v_chain    = tk.StringVar(value="")
        self.v_outdir   = tk.StringVar(value=str(_HERE))
        self.v_outname  = tk.StringVar(value="protein_liabilities")
        self.v_nohtml   = tk.BooleanVar(value=False)
        self.v_seqlen   = tk.StringVar(value="")
        self.v_cdr_scheme = tk.StringVar(value="kabat")
        self.v_status   = tk.StringVar(value="Ready")

        # Report filter state
        self.v_rpt_risk   = tk.StringVar(value="All")
        self.v_rpt_cat    = tk.StringVar(value="All")
        self.v_rpt_search = tk.StringVar()

        self._html_path     = None
        self._running       = False
        self._last_pdb_used = False
        self._log_q: queue.Queue = queue.Queue()

        # Report data store
        self._rpt_all_rows: list = []
        self._rpt_sort_col  = None
        self._rpt_sort_rev  = False

        # 3D viewer state
        self._3d_results:   list         = []
        self._3d_pdb_text:  str | None   = None

        self._build_styles()
        self._build_header()
        self._build_notebook()
        self._build_statusbar()
        self.root.after(60, self._drain_queue)

        if not CORE_OK:
            self.root.after(300, lambda: messagebox.showerror(
                "Module not found",
                f"Cannot import protein_liability_analyzer_v2.py\n\n"
                f"Put it in the same folder as this script:\n{_HERE}\n\nError: {CORE_ERR}"
            ))

    # ── Window helper ─────────────────────────────────────────────────────────
    def _centre(self, w, h):
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        self.root.geometry(f"{w}x{h}+{max(0,(sw-w)//2)}+{max(0,(sh-h)//2-30)}")

    # ── Styles ────────────────────────────────────────────────────────────────
    def _build_styles(self):
        s = ttk.Style()
        s.theme_use("clam")
        s.configure("TFrame",       background=LIGHT)
        s.configure("TLabel",       background=LIGHT, foreground=DARK, font=F_BODY)
        s.configure("TEntry",       fieldbackground=WHITE, foreground=DARK,
                    insertcolor=DARK, font=F_BODY, padding=4)
        s.configure("TCheckbutton", background=WHITE, foreground=DARK, font=F_BODY)
        s.configure("TRadiobutton", background=WHITE, foreground=DARK, font=F_BODY)
        s.configure("TProgressbar", troughcolor=BORDER, background=ACCENT, thickness=6)
        s.configure("TNotebook",    background=LIGHT, borderwidth=0)
        s.configure("TNotebook.Tab", padding=[16, 8], font=F_BODY,
                    background=BORDER, foreground=GREY)
        s.map("TNotebook.Tab",
              background=[("selected", WHITE)],
              foreground=[("selected", DARK)],
              expand=[("selected", [1,1,1,0])])
        # Treeview
        s.configure("Treeview",
                    background=WHITE, fieldbackground=WHITE,
                    foreground=DARK, font=F_SMALL, rowheight=22)
        s.configure("Treeview.Heading",
                    background=ACCENT, foreground=WHITE,
                    font=("Segoe UI", 9, "bold"), relief="flat")
        s.map("Treeview",
              background=[("selected", ACCENT2)],
              foreground=[("selected", WHITE)])

    # ── Header bar ────────────────────────────────────────────────────────────
    def _build_header(self):
        hdr = tk.Frame(self.root, bg=DARK, height=66)
        hdr.pack(fill="x")
        hdr.pack_propagate(False)
        inner = tk.Frame(hdr, bg=DARK)
        inner.pack(fill="both", expand=True, padx=20, pady=10)

        # LSC logo drawn with Canvas primitives (no image dependency)
        lw, W, H, r, pad = 2, 74, 42, 19, 2
        logo = tk.Canvas(inner, width=W, height=H, bg=DARK, highlightthickness=0)
        logo.pack(side="left", padx=(0, 16))
        logo.create_arc(pad, pad, pad + 2*r, H - pad,
                        start=90, extent=180, style="arc", outline=WHITE, width=lw)
        logo.create_arc(W - pad - 2*r, pad, W - pad, H - pad,
                        start=-90, extent=180, style="arc", outline=WHITE, width=lw)
        logo.create_line(pad + r, pad, W - pad - r, pad, fill=WHITE, width=lw)
        logo.create_line(pad + r, H - pad, W - pad - r, H - pad, fill=WHITE, width=lw)
        logo.create_text(W // 2, H // 2, text="LSC", fill=WHITE,
                         font=("Futura", 17, "bold"))

        text_frame = tk.Frame(inner, bg=DARK)
        text_frame.pack(side="left")
        tk.Label(text_frame, text="Protein Liability Analyzer  v2  -  Auto-HOS",
                 bg=DARK, fg=WHITE, font=F_TITLE).pack(anchor="w")
        tk.Label(text_frame,
                 text="HOS always performed  ·  Letarte Scientific Consulting  ·  "
                      "All computation local - sequences never transmitted",
                 bg=DARK, fg="#7788aa", font=("Segoe UI", 8)).pack(anchor="w")

    # ── Notebook ──────────────────────────────────────────────────────────────
    def _build_notebook(self):
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill="both", expand=True, padx=10, pady=(6, 2))
        t1 = ttk.Frame(self.nb)
        t2 = ttk.Frame(self.nb)
        t3 = ttk.Frame(self.nb)
        t4 = ttk.Frame(self.nb)
        self.nb.add(t1, text="  Analysis Setup  ")
        self.nb.add(t2, text="  Results & Log  ")
        self.nb.add(t3, text="  Report  ")
        self.nb.add(t4, text="  3D Structure  ")
        self._build_setup(t1)
        self._build_results(t2)
        self._build_report(t3)
        self._build_3d_tab(t4)

    # ── Card helper ───────────────────────────────────────────────────────────
    def _card(self, parent, title, builder):
        wrap = tk.Frame(parent, bg=LIGHT)
        wrap.pack(fill="x", padx=14, pady=6)
        bar = tk.Frame(wrap, bg=ACCENT, height=28)
        bar.pack(fill="x")
        bar.pack_propagate(False)
        tk.Label(bar, text=title, bg=ACCENT, fg=WHITE,
                 font=("Segoe UI", 10, "bold"), padx=12).pack(side="left", pady=4)
        body = tk.Frame(wrap, bg=WHITE, padx=14, pady=12)
        body.pack(fill="x")
        builder(body)
        tk.Frame(wrap, bg=BORDER, height=1).pack(fill="x")

    def _flatbtn(self, parent, text, cmd, bg=BORDER, fg=DARK, state="normal", **kw):
        b = tk.Button(parent, text=text, command=cmd, bg=bg, fg=fg,
                      font=F_SMALL, relief="flat", bd=0, padx=10, pady=5,
                      state=state, cursor="hand2" if state == "normal" else "arrow",
                      activebackground=ACCENT2, activeforeground=WHITE, **kw)
        return b

    # ── Setup tab ─────────────────────────────────────────────────────────────
    def _build_setup(self, parent):
        cvs = tk.Canvas(parent, bg=LIGHT, highlightthickness=0)
        vsb = ttk.Scrollbar(parent, orient="vertical", command=cvs.yview)
        cvs.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        cvs.pack(side="left", fill="both", expand=True)
        inner = tk.Frame(cvs, bg=LIGHT)
        win   = cvs.create_window((0,0), window=inner, anchor="nw")
        cvs.bind("<Configure>", lambda e: cvs.itemconfig(win, width=e.width))
        inner.bind("<Configure>", lambda e: cvs.configure(scrollregion=cvs.bbox("all")))
        cvs.bind_all("<MouseWheel>", lambda e: cvs.yview_scroll(int(-1*(e.delta/120)), "units"))

        self._card(inner, "Input",
                   self._input_card)
        self._card(inner, "Higher-Order Structure Analysis (always performed)",
                   self._hos_card)
        self._card(inner, "Output",
                   self._output_card)

        run_wrap = tk.Frame(inner, bg=LIGHT)
        run_wrap.pack(fill="x", padx=14, pady=(10,18))

        self.progress = ttk.Progressbar(run_wrap, mode="indeterminate")

        self.run_btn = tk.Button(
            run_wrap, text="Run Analysis", command=self._start,
            bg=ACCENT, fg=ACCENT2, font=("Segoe UI", 13, "bold"),
            relief="flat", bd=0, padx=24, pady=14, cursor="hand2",
            activebackground=ACCENT2, activeforeground=WHITE,
        )
        self.run_btn.pack(fill="x")
        self.run_btn.bind("<Enter>", lambda _: self.run_btn.configure(bg=ACCENT2))
        self.run_btn.bind("<Leave>", lambda _: self.run_btn.configure(bg=ACCENT))

    def _input_card(self, p):
        row = tk.Frame(p, bg=WHITE)
        row.pack(fill="x", pady=(0,8))
        for lbl, val in [("Paste / type sequence", "paste"),
                         ("Load from FASTA file",  "file")]:
            tk.Radiobutton(row, text=lbl, variable=self.v_mode, value=val,
                           bg=WHITE, fg=DARK, font=F_BODY, activebackground=WHITE,
                           command=self._toggle_mode).pack(side="left", padx=(0,20))

        self.paste_frame = tk.Frame(p, bg=WHITE)
        hdr_row = tk.Frame(self.paste_frame, bg=WHITE)
        hdr_row.pack(fill="x", pady=(0,4))
        tk.Label(hdr_row, text="Sequence - FASTA or plain amino-acid string:",
                 bg=WHITE, fg=GREY, font=F_SMALL).pack(side="left")
        tk.Label(hdr_row, textvariable=self.v_seqlen,
                 bg=WHITE, fg=ACCENT2, font=F_SMALL).pack(side="right")
        self.seq_text = tk.Text(
            self.paste_frame, height=10, wrap="word",
            font=("Courier New", 10), fg=DARK, bg="#f8f9fa",
            relief="flat", bd=0, highlightthickness=1,
            highlightbackground=BORDER, highlightcolor=ACCENT,
            insertbackground=DARK,
        )
        self.seq_text.pack(fill="x")
        self.seq_text.bind("<KeyRelease>", self._update_len)

        self.file_frame = tk.Frame(p, bg=WHITE)
        tk.Label(self.file_frame, text="FASTA file:",
                 bg=WHITE, fg=DARK, font=F_BODY, width=12, anchor="w").pack(side="left")
        tk.Entry(self.file_frame, textvariable=self.v_fasta, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left", fill="x", expand=True, padx=(0,8))
        self._flatbtn(self.file_frame, "Browse...", self._pick_fasta).pack(side="left")
        self._toggle_mode()

    def _hos_card(self, p):
        status_frame = tk.Frame(p, bg="#eaf6ea")
        status_frame.pack(fill="x", pady=(0, 10))
        tk.Label(status_frame,
                 text="HOS analysis is always performed automatically:",
                 bg="#eaf6ea", fg="#1a5c1a", font=F_BOLD,
                 padx=10, pady=6).pack(anchor="w")
        modes_frame = tk.Frame(status_frame, bg="#eaf6ea")
        modes_frame.pack(fill="x", padx=20, pady=(0, 6))
        tk.Label(modes_frame,
                 text="No PDB file  ->  Chou-Fasman secondary structure prediction "
                      "+ sequence-based RSA estimation (pure Python, entirely local)",
                 bg="#eaf6ea", fg="#1a5c1a", font=F_SMALL,
                 wraplength=740, justify="left").pack(anchor="w", pady=2)
        tk.Label(modes_frame,
                 text="PDB file provided  ->  HELIX/SHEET records for secondary structure "
                      "+ Shrake-Rupley SASA (all-atom, 1.4 Å probe, no BioPython required)",
                 bg="#eaf6ea", fg="#1a5c1a", font=F_SMALL,
                 wraplength=740, justify="left").pack(anchor="w", pady=2)

        pdb_row = tk.Frame(p, bg=WHITE)
        pdb_row.pack(fill="x", pady=(0,6))
        tk.Label(pdb_row, text="PDB file\n(optional):", bg=WHITE, fg=DARK, font=F_BODY,
                 width=13, justify="left", anchor="w").pack(side="left")
        tk.Entry(pdb_row, textvariable=self.v_pdb, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left", fill="x", expand=True, padx=(0,8))
        self._flatbtn(pdb_row, "Browse...", self._pick_pdb).pack(side="left")

        ch_row = tk.Frame(p, bg=WHITE)
        ch_row.pack(fill="x", pady=(0,8))
        tk.Label(ch_row, text="Chain ID:", bg=WHITE, fg=DARK, font=F_BODY,
                 width=13, anchor="w").pack(side="left")
        tk.Entry(ch_row, textvariable=self.v_chain, width=6, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left")
        tk.Label(ch_row, text="  Leave blank for automatic selection",
                 bg=WHITE, fg=GREY, font=F_SMALL).pack(side="left")

        cdr_row = tk.Frame(p, bg=WHITE)
        cdr_row.pack(fill="x", pady=(6, 2))
        tk.Label(cdr_row, text="CDR scheme:", bg=WHITE, fg=DARK, font=F_BODY,
                 width=13, anchor="w").pack(side="left")
        for scheme, label in [("kabat", "Kabat"), ("chothia", "Chothia"), ("imgt", "IMGT")]:
            tk.Radiobutton(
                cdr_row, text=label, variable=self.v_cdr_scheme, value=scheme,
                bg=WHITE, fg=DARK, font=F_BODY, activebackground=WHITE,
            ).pack(side="left", padx=(0, 16))
        tk.Label(cdr_row, text="(antibody V-domain sequences only)",
                 bg=WHITE, fg=GREY, font=F_SMALL).pack(side="left")

        info = tk.Frame(p, bg="#e8f4fd")
        info.pack(fill="x", pady=(6, 0))
        tk.Label(info,
                 text="Supplying a PDB file upgrades the analysis from theoretical "
                      "(sequence-based) to experimental (structure-based) mode. "
                      "Either way, all computation runs locally - your sequence is never transmitted.",
                 bg="#e8f4fd", fg="#154360", font=F_SMALL,
                 wraplength=760, justify="left", padx=10, pady=8).pack(anchor="w")

    def _output_card(self, p):
        d_row = tk.Frame(p, bg=WHITE)
        d_row.pack(fill="x", pady=(0,6))
        tk.Label(d_row, text="Save folder:", bg=WHITE, fg=DARK, font=F_BODY,
                 width=13, anchor="w").pack(side="left")
        tk.Entry(d_row, textvariable=self.v_outdir, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left", fill="x", expand=True, padx=(0,8))
        self._flatbtn(d_row, "Browse...", self._pick_outdir).pack(side="left")

        n_row = tk.Frame(p, bg=WHITE)
        n_row.pack(fill="x", pady=(0,4))
        tk.Label(n_row, text="Report name:", bg=WHITE, fg=DARK, font=F_BODY,
                 width=13, anchor="w").pack(side="left")
        tk.Entry(n_row, textvariable=self.v_outname, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left", fill="x", expand=True)
        tk.Label(n_row, text=".html", bg=WHITE, fg=GREY, font=F_BODY).pack(side="left", padx=4)

        tk.Checkbutton(p, text="Console output only - skip HTML report generation",
                       variable=self.v_nohtml, bg=WHITE, fg=GREY, font=F_SMALL,
                       activebackground=WHITE).pack(anchor="w", pady=(6,0))

    # ── Results / Log tab ─────────────────────────────────────────────────────
    def _build_results(self, parent):
        toolbar = tk.Frame(parent, bg=LIGHT, pady=5)
        toolbar.pack(fill="x", padx=10)
        tk.Label(toolbar, text="Analysis Log", bg=LIGHT, fg=DARK, font=F_BOLD).pack(side="left")

        right = tk.Frame(toolbar, bg=LIGHT)
        right.pack(side="right")
        self.btn_html = self._flatbtn(right, "Open HTML Report", self._open_html,
                                      bg=ACCENT, fg=WHITE, state="disabled")
        self.btn_html.pack(side="left", padx=(0,6))
        self.btn_folder = self._flatbtn(right, "Open Folder", self._open_folder,
                                        state="disabled")
        self.btn_folder.pack(side="left", padx=(0,14))
        self._flatbtn(right, "Copy log", self._copy_log).pack(side="left", padx=(0,4))
        self._flatbtn(right, "Clear",    self._clear_log).pack(side="left")

        log_wrap = tk.Frame(parent, bg=LOG_BG)
        log_wrap.pack(fill="both", expand=True, padx=10, pady=(0,6))
        vsb = ttk.Scrollbar(log_wrap, orient="vertical")
        vsb.pack(side="right", fill="y")
        self.log = tk.Text(
            log_wrap, font=F_MONO, bg=LOG_BG, fg=LOG_FG,
            relief="flat", bd=0, state="disabled", wrap="word",
            padx=12, pady=10, cursor="arrow",
            yscrollcommand=vsb.set,
        )
        self.log.pack(fill="both", expand=True)
        vsb.configure(command=self.log.yview)

        self.log.tag_configure("normal",  foreground=LOG_FG)
        self.log.tag_configure("dim",     foreground="#484f58")
        self.log.tag_configure("cyan",    foreground="#79c0ff")
        self.log.tag_configure("green",   foreground="#56d364")
        self.log.tag_configure("yellow",  foreground="#e3b341")
        self.log.tag_configure("red",     foreground="#ff7b72",
                               font=("Courier New", 9, "bold"))
        self.log.tag_configure("bold",    foreground="#f0f6fc",
                               font=("Courier New", 9, "bold"))
        self.log.tag_configure("header",  foreground="#f0f6fc",
                               font=("Courier New", 9, "bold"),
                               background="#161b22")

        self._log(
            "\n  Welcome to Protein Liability Analyzer v2  ·  Auto-HOS Edition\n"
            "  ------------------------------------------------------------------------------------------------------------------------------------\n"
            "  HOS analysis runs automatically on every sequence:\n"
            "    - No PDB  ->  Chou-Fasman + sequence-based RSA (theoretical)\n"
            "    - PDB     ->  HELIX/SHEET + Shrake-Rupley SASA (experimental)\n\n"
            "  Set up your analysis in the Setup tab, then click Run Analysis.\n"
            "  Results will appear here and in the Report tab.\n\n",
            "dim"
        )

    # ── Report tab ────────────────────────────────────────────────────────────
    def _build_report(self, parent):
        # ── Summary bar (populated after analysis) ────────────────────────
        self.rpt_summary = tk.Frame(parent, bg=LIGHT, pady=6)
        self.rpt_summary.pack(fill="x", padx=10)
        self._rpt_placeholder = tk.Label(
            self.rpt_summary,
            text="No analysis results yet.  Run an analysis to populate this report.",
            bg=LIGHT, fg=GREY, font=F_BODY)
        self._rpt_placeholder.pack(pady=6)

        # ── Filter toolbar ────────────────────────────────────────────────
        fbar = tk.Frame(parent, bg=LIGHT, pady=4)
        fbar.pack(fill="x", padx=10)

        tk.Label(fbar, text="Risk:", bg=LIGHT, fg=DARK, font=F_SMALL).pack(side="left")
        risk_cb = ttk.Combobox(fbar, textvariable=self.v_rpt_risk,
                               values=["All", "High", "Medium", "Low", "Info"],
                               width=9, state="readonly", font=F_SMALL)
        risk_cb.pack(side="left", padx=(4, 12))
        risk_cb.bind("<<ComboboxSelected>>", lambda _: self._filter_report())

        tk.Label(fbar, text="Category:", bg=LIGHT, fg=DARK, font=F_SMALL).pack(side="left")
        self.rpt_cat_cb = ttk.Combobox(fbar, textvariable=self.v_rpt_cat,
                                       values=["All"], width=14,
                                       state="readonly", font=F_SMALL)
        self.rpt_cat_cb.pack(side="left", padx=(4, 12))
        self.rpt_cat_cb.bind("<<ComboboxSelected>>", lambda _: self._filter_report())

        tk.Label(fbar, text="Search:", bg=LIGHT, fg=DARK, font=F_SMALL).pack(side="left")
        self.v_rpt_search.trace_add("write", lambda *_: self._filter_report())
        tk.Entry(fbar, textvariable=self.v_rpt_search, width=20, font=F_SMALL,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left", padx=(4, 12))

        self.btn_export_csv = self._flatbtn(
            fbar, "Export CSV", self._export_csv, state="disabled")
        self.btn_export_csv.pack(side="right")

        self.rpt_count_lbl = tk.Label(fbar, text="", bg=LIGHT, fg=GREY, font=F_SMALL)
        self.rpt_count_lbl.pack(side="right", padx=10)

        # ── Treeview ──────────────────────────────────────────────────────
        # Column order mirrors the HTML "All Findings" table
        cols    = ("pos", "residue", "category", "risk", "secstr", "exposure", "rsa",
                   "type_label", "note")
        hdrs    = ("Pos.", "Res.", "Category", "Seq Risk", "Sec. Str.",
                   "Exposure", "RSA %", "Liability Type", "Note")
        widths  = (50, 55, 110, 80, 72, 88, 52, 210, 180)
        anchors = ("center", "center", "w", "center", "center",
                   "center", "center", "w", "w")

        tree_outer = tk.Frame(parent, bg=LIGHT)
        tree_outer.pack(fill="both", expand=True, padx=10, pady=(2,0))

        vsb = ttk.Scrollbar(tree_outer, orient="vertical")
        hsb = ttk.Scrollbar(tree_outer, orient="horizontal")
        vsb.pack(side="right", fill="y")
        hsb.pack(side="bottom", fill="x")

        self.rpt_tree = ttk.Treeview(
            tree_outer, columns=cols, show="headings",
            yscrollcommand=vsb.set, xscrollcommand=hsb.set,
            selectmode="browse")
        vsb.configure(command=self.rpt_tree.yview)
        hsb.configure(command=self.rpt_tree.xview)
        self.rpt_tree.pack(fill="both", expand=True)

        for col, hdr, w, anc in zip(cols, hdrs, widths, anchors):
            self.rpt_tree.heading(col, text=hdr,
                                  command=lambda c=col: self._sort_report(c))
            self.rpt_tree.column(col, width=w, minwidth=34, anchor=anc)

        # Row color tags - match HTML RISK_COLORS palette
        self.rpt_tree.tag_configure("high",
            background="#fde8e8", foreground="#c0392b")
        self.rpt_tree.tag_configure("medium",
            background="#fef0e6", foreground="#d35400")
        self.rpt_tree.tag_configure("low",
            background="#e8f8f5", foreground="#0e6655")
        self.rpt_tree.tag_configure("info",
            background="#d6eaf8", foreground="#154360")
        self.rpt_tree.bind("<<TreeviewSelect>>", self._on_rpt_select)

        # ── Structural columns note bar (mirrors HTML blue footer) ────────
        struct_note = tk.Frame(parent, bg="#e3f2fd")
        struct_note.pack(fill="x", padx=10, pady=(0,2))
        tk.Label(struct_note,
                 text="Sec. Str.:  H = α Helix  ·  E = β Sheet  ·  T = Turn  ·  C = Coil  "
                      "  |  Exposure from Chou-Fasman RSA (sequence) or Shrake-Rupley SASA (PDB)  "
                      "  |  Click column headers to sort  ·  Click a row for details",
                 bg="#e3f2fd", fg="#1565c0", font=("Segoe UI", 8),
                 padx=8, pady=4).pack(anchor="w")

        # ── Detail panel ──────────────────────────────────────────────────
        sep = tk.Frame(parent, bg=BORDER, height=1)
        sep.pack(fill="x", padx=10, pady=(4,0))

        detail_outer = tk.Frame(parent, bg=LIGHT)
        detail_outer.pack(fill="x", padx=10, pady=(0,6))

        tk.Label(detail_outer, text="Selected liability:",
                 bg=LIGHT, fg=GREY, font=F_SMALL).pack(anchor="w", pady=(4,2))

        self.rpt_detail = tk.Text(
            detail_outer, height=4,
            font=("Courier New", 10), bg=WHITE, fg=DARK,
            relief="flat", bd=0, highlightthickness=1,
            highlightbackground=BORDER, padx=10, pady=6,
            wrap="word", state="disabled")
        self.rpt_detail.pack(fill="x")

        self.rpt_detail.tag_configure("motif",   background="#ffe066", foreground="#000000")
        self.rpt_detail.tag_configure("dim",     foreground=GREY,   font=("Segoe UI", 9))
        self.rpt_detail.tag_configure("bold",    foreground=DARK,   font=("Segoe UI", 9, "bold"))
        self.rpt_detail.tag_configure("high_t",  foreground="#c0392b", font=("Segoe UI", 9, "bold"))
        self.rpt_detail.tag_configure("med_t",   foreground="#d35400", font=("Segoe UI", 9, "bold"))
        self.rpt_detail.tag_configure("low_t",   foreground="#0e6655", font=("Segoe UI", 9, "bold"))
        self.rpt_detail.tag_configure("info_t",  foreground="#154360", font=("Segoe UI", 9, "bold"))
        self.rpt_detail.tag_configure("desc",    foreground="#555555", font=("Segoe UI", 9, "italic"))

    # ── Populate report after analysis ────────────────────────────────────────
    # SS single-char -> human label (Chou-Fasman or PDB)
    _SS_LABEL = {"H": "Helix (α)", "E": "Sheet (β)", "T": "Turn", "C": "Coil",
                 "--": "--", "?": "--"}

    def _populate_report(self, results):
        # Clear treeview and summary
        for item in self.rpt_tree.get_children():
            self.rpt_tree.delete(item)
        for w in self.rpt_summary.winfo_children():
            w.destroy()

        self._rpt_all_rows = []
        categories_seen = set()

        cdr_scheme = self.v_cdr_scheme.get()

        for entry in results:
            name = entry["name"]
            seq  = entry["seq"]

            # Build per-position CDR map for this chain
            cdrs      = pla.annotate_cdrs(seq, scheme=cdr_scheme)
            cdr_at    = {}
            for s, e, cdr_name in cdrs:
                for i in range(s, e):
                    cdr_at[i] = cdr_name

            for f in entry["findings"]:
                pos0  = f.get("pos0", 0)
                pos1  = f.get("pos1", pos0 + 1)
                res   = f.get("residues", "")

                # SS and RSA written directly onto the finding by add_structural_context
                ss_char = f.get("ss", "--") or "--"
                ss_lbl  = self._SS_LABEL.get(ss_char, ss_char)
                rsa_val = f.get("rsa", None)
                rsa_str = f"{rsa_val * 100:.0f}" if rsa_val is not None else "--"
                exposure = f.get("exposure_class", "--")

                category  = f.get("category", "")
                risk_key  = f.get("risk", "info").lower()
                base_lbl  = re.sub(
                    r'\s*[–—-]\s*(High|Medium|Low|Info)\s*Risk\s*',
                    ' ', f.get("label", f.get("type", ""))
                ).strip()
                cdr_name  = cdr_at.get(pos0)
                type_lbl  = f"{base_lbl} · {cdr_name}" if cdr_name else base_lbl
                note      = f.get("note", "") or ""

                categories_seen.add(category)

                row = {
                    "chain":      name,
                    "pos":        pos1,
                    "_pos0":      pos0,
                    "residue":    res,
                    "category":   category,
                    "risk":       risk_key,
                    "secstr":     ss_lbl,
                    "exposure":   exposure,
                    "rsa":        rsa_str,
                    "type_label": type_lbl,
                    "note":       note,
                    # extras for detail panel
                    "_seq": seq,
                    "_f":   f,
                }
                self._rpt_all_rows.append(row)

        # Update category combobox with seen categories
        cat_list = ["All"] + sorted(categories_seen)
        self.rpt_cat_cb.configure(values=cat_list)
        self.v_rpt_cat.set("All")

        # ── Summary bar ──────────────────────────────────────────────────
        total  = len(self._rpt_all_rows)
        counts = {"high": 0, "medium": 0, "low": 0, "info": 0}
        for r in self._rpt_all_rows:
            counts[r["risk"]] = counts.get(r["risk"], 0) + 1

        tk.Label(self.rpt_summary,
                 text=f"  {total} liabilities  ·  {len(results)} chain(s)  ",
                 bg=LIGHT, fg=DARK, font=F_BOLD).pack(side="left", padx=(0, 10))

        badge_specs = [
            ("High",   counts["high"],   "#c0392b", "#fde8e8", "#e74c3c"),
            ("Medium", counts["medium"], "#d35400", "#fef0e6", "#e67e22"),
            ("Low",    counts["low"],    "#0e6655", "#e8f8f5", "#1abc9c"),
            ("Info",   counts["info"],   "#154360", "#d6eaf8", "#2e86c1"),
        ]
        for label, cnt, fg, bg, border in badge_specs:
            badge = tk.Frame(self.rpt_summary, bg=bg, padx=10, pady=3,
                             highlightthickness=1, highlightbackground=border)
            badge.pack(side="left", padx=4)
            tk.Label(badge, text=f"{label}:  {cnt}", bg=bg, fg=fg, font=F_BOLD).pack()

        src_lbl = ("Experimental (PDB)" if self._last_pdb_used
                   else "Theoretical (Chou-Fasman + RSA)")
        tk.Label(self.rpt_summary, text=f"  HOS source: {src_lbl}",
                 bg=LIGHT, fg=GREY, font=F_SMALL).pack(side="right", padx=6)

        # ── Render ────────────────────────────────────────────────────────
        self._render_rpt_rows(self._rpt_all_rows)
        self.btn_export_csv.configure(state="normal", cursor="hand2")

    def _render_rpt_rows(self, rows):
        for item in self.rpt_tree.get_children():
            self.rpt_tree.delete(item)
        for r in rows:
            vals = (r["pos"], r["residue"], r["category"],
                    r["risk"].upper(), r["secstr"], r["exposure"],
                    r["rsa"], r["type_label"], r["note"])
            self.rpt_tree.insert("", tk.END, values=vals, tags=(r["risk"],))
        n = len(rows)
        self.rpt_count_lbl.configure(
            text=f"{n} finding{'s' if n != 1 else ''}")

    def _get_filtered_rows(self):
        risk_f = self.v_rpt_risk.get().lower()
        cat_f  = self.v_rpt_cat.get()
        srch   = self.v_rpt_search.get().lower()
        return [
            r for r in self._rpt_all_rows
            if (risk_f == "all" or r["risk"] == risk_f)
            and (cat_f == "All" or r["category"] == cat_f)
            and (not srch
                 or srch in r["chain"].lower()
                 or srch in r["category"].lower()
                 or srch in r["type_label"].lower()
                 or srch in r["residue"].lower()
                 or srch in r["note"].lower())
        ]

    def _filter_report(self):
        self._render_rpt_rows(self._get_filtered_rows())

    def _sort_report(self, col):
        RISK_ORDER = {"high": 0, "medium": 1, "low": 2, "info": 3}
        if self._rpt_sort_col == col:
            self._rpt_sort_rev = not self._rpt_sort_rev
        else:
            self._rpt_sort_col = col
            self._rpt_sort_rev = False

        def key(r):
            v = r.get(col, "")
            if col == "pos":
                try: return int(v)
                except: return 0
            if col == "risk":
                return RISK_ORDER.get(str(v).lower(), 9)
            if col == "rsa":
                try: return float(v)
                except: return 9999
            return str(v).lower()

        self._render_rpt_rows(
            sorted(self._get_filtered_rows(), key=key, reverse=self._rpt_sort_rev))

    def _on_rpt_select(self, _=None):
        sel = self.rpt_tree.selection()
        if not sel:
            return
        vals = self.rpt_tree.item(sel[0], "values")
        # vals order: pos, residue, category, risk, secstr, exposure, rsa, type_label, note
        try:
            pos_val  = int(vals[0])
            type_lbl = vals[7]
            row = next(r for r in self._rpt_all_rows
                       if r["pos"] == pos_val and r["type_label"] == type_lbl)
        except (StopIteration, ValueError, IndexError):
            return

        f    = row["_f"]
        seq  = row["_seq"]
        pos0 = row["_pos0"]
        res  = row["residue"]

        # Wide context window
        win_start = max(0, pos0 - 11)
        win_end   = min(len(seq), pos0 + len(res) + 11)
        ctx       = seq[win_start:win_end]
        mot_s     = pos0 - win_start

        risk_tag = {"high": "high_t", "medium": "med_t",
                    "low":  "low_t",  "info":   "info_t"}.get(row["risk"], "bold")

        self.rpt_detail.configure(state="normal")
        self.rpt_detail.delete("1.0", tk.END)

        # Line 1: liability type + position metadata
        self.rpt_detail.insert(tk.END, row["type_label"], "bold")
        self.rpt_detail.insert(tk.END,
            f"  ·  Pos {row['pos']}  ·  Res: {res}  ·  Chain: {row['chain']}"
            f"  ·  Risk: ", "dim")
        self.rpt_detail.insert(tk.END, row["risk"].upper(), risk_tag)
        self.rpt_detail.insert(tk.END,
            f"  ·  {row['secstr']}  ·  Exposure: {row['exposure']}"
            f"  ·  RSA: {row['rsa']}%\n", "dim")

        # Line 2: sequence context with residue(s) highlighted in yellow
        self.rpt_detail.insert(tk.END, "  ...")
        self.rpt_detail.insert(tk.END, ctx[:mot_s])
        self.rpt_detail.insert(tk.END, ctx[mot_s:mot_s + len(res)], "motif")
        self.rpt_detail.insert(tk.END, ctx[mot_s + len(res):])
        self.rpt_detail.insert(tk.END, "...\n")

        # Line 3: scientific description from the liability catalogue
        desc = f.get("description", "")
        if desc:
            self.rpt_detail.insert(tk.END, f"  {desc}", "desc")

        self.rpt_detail.configure(state="disabled")

    def _export_csv(self):
        if not self._rpt_all_rows:
            return
        p = filedialog.asksaveasfilename(
            title="Export findings as CSV",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")],
            initialdir=self.v_outdir.get(),
            initialfile="liabilities_export.csv")
        if not p:
            return
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["Chain", "Pos", "Residue", "Category",
                        "Risk", "Sec.Str.", "Exposure", "RSA%",
                        "Liability Type", "Note"])
            for r in self._rpt_all_rows:
                w.writerow([r["chain"], r["pos"], r["residue"], r["category"],
                            r["risk"].upper(), r["secstr"], r["exposure"],
                            r["rsa"], r["type_label"], r["note"]])
        self.v_status.set(f"CSV exported -> {Path(p).name}")

    # ── Status bar ────────────────────────────────────────────────────────────
    def _build_statusbar(self):
        bar = tk.Frame(self.root, bg=DARK, height=22)
        bar.pack(fill="x", side="bottom")
        bar.pack_propagate(False)
        tk.Label(bar, textvariable=self.v_status, bg=DARK, fg="#7788aa",
                 font=("Segoe UI", 8), padx=12).pack(side="left", pady=2)
        tk.Label(bar, text="Local only - sequences never transmitted",
                 bg=DARK, fg="#445566", font=("Segoe UI", 8), padx=12).pack(side="right", pady=2)

    # ── UI callbacks ──────────────────────────────────────────────────────────
    def _toggle_mode(self):
        if self.v_mode.get() == "paste":
            self.file_frame.pack_forget()
            self.paste_frame.pack(fill="x", pady=(0,6))
        else:
            self.paste_frame.pack_forget()
            self.file_frame.pack(fill="x", pady=(0,6))

    def _update_len(self, _=None):
        import re as _re
        raw = self.seq_text.get("1.0", tk.END)
        count = sum(len(_re.sub(r"[^A-Za-z]", "", l))
                    for l in raw.splitlines() if not l.startswith(">"))
        self.v_seqlen.set(f"{count} residues" if count else "")

    def _pick_fasta(self):
        p = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA / text", "*.fasta *.fa *.fas *.txt"), ("All files", "*.*")],
            initialdir=self.v_outdir.get())
        if p:
            self.v_fasta.set(p); self.v_mode.set("file")
            self.v_outdir.set(str(Path(p).parent)); self._toggle_mode()

    def _pick_pdb(self):
        p = filedialog.askopenfilename(
            title="Select PDB file",
            filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
            initialdir=self.v_outdir.get())
        if p:
            self.v_pdb.set(p)

    def _pick_outdir(self):
        d = filedialog.askdirectory(title="Select output folder",
                                    initialdir=self.v_outdir.get())
        if d: self.v_outdir.set(d)

    # ── Log helpers ───────────────────────────────────────────────────────────
    def _log(self, msg: str, tag: str = "normal"):
        self.log.configure(state="normal")
        self.log.insert(tk.END, msg, tag)
        self.log.see(tk.END)
        self.log.configure(state="disabled")

    def _auto_tag(self, msg: str) -> str:
        if any(x in msg for x in ("Error:", "error:", "failed")):           return "red"
        if any(x in msg for x in ("[Warning]", "Warning")):        return "yellow"
        if any(x in msg for x in ("done", "saved ->", "complete", "Report saved")): return "green"
        if any(x in msg for x in ("Computing", "Predicting",
                                   "Estimating", "Parsing", "Analysing")): return "cyan"
        if "-" * 5 in msg:                                                return "header"
        if msg.strip().startswith("Pos "):                                return "bold"
        return "normal"

    def _drain_queue(self):
        try:
            while True:
                tag, msg = self._log_q.get_nowait()
                if tag == "normal": tag = self._auto_tag(msg)
                self._log(msg, tag)
        except queue.Empty:
            pass
        self.root.after(60, self._drain_queue)

    def _clear_log(self):
        self.log.configure(state="normal")
        self.log.delete("1.0", tk.END)
        self.log.configure(state="disabled")

    def _copy_log(self):
        self.root.clipboard_clear()
        self.root.clipboard_append(self.log.get("1.0", tk.END))
        self.v_status.set("Log copied to clipboard")

    def _open_html(self):
        if self._html_path and Path(self._html_path).exists():
            webbrowser.open(Path(self._html_path).as_uri())
        else:
            messagebox.showinfo("Not available", "HTML report not found.")

    def _open_folder(self):
        d = self.v_outdir.get()
        if not Path(d).exists():
            messagebox.showinfo("Not found", f"Folder not found:\n{d}"); return
        try:
            sys_name = platform.system()
            if   sys_name == "Darwin":  subprocess.Popen(["open",     d])
            elif sys_name == "Windows": os.startfile(d)
            else:                       subprocess.Popen(["xdg-open", d])
        except Exception as exc:
            messagebox.showerror("Error", str(exc))

    # ── Analysis ──────────────────────────────────────────────────────────────
    def _start(self):
        if self._running: return
        if not CORE_OK:
            messagebox.showerror("Error", "Analysis module not loaded."); return

        if self.v_mode.get() == "paste":
            if not self.seq_text.get("1.0", tk.END).strip():
                messagebox.showwarning("No input", "Please paste a sequence."); return
        else:
            fp = self.v_fasta.get().strip()
            if not fp:
                messagebox.showwarning("No input", "Please select a FASTA file."); return
            if not Path(fp).exists():
                messagebox.showwarning("File not found", f"FASTA file not found:\n{fp}"); return

        pdb = self.v_pdb.get().strip()
        if pdb and not Path(pdb).exists():
            messagebox.showwarning("PDB not found", f"PDB file not found:\n{pdb}"); return

        self.nb.select(1)  # Switch to log tab while running
        self._log(f"\n{chr(45)*62}\n  Started: {datetime.now():%Y-%m-%d %H:%M:%S}\n", "cyan")

        self.run_btn.configure(state="disabled", text="Running...", bg="#555")
        self.progress.pack(fill="x", pady=(0,8), before=self.run_btn)
        self.progress.start(12)
        self._running   = True
        self._html_path = None
        self.btn_html.configure(state="disabled", cursor="arrow")
        self.btn_folder.configure(state="disabled", cursor="arrow")
        self.v_status.set("Running analysis...")

        threading.Thread(target=self._worker, daemon=True).start()

    def _worker(self):
        q = self._log_q
        old_out, old_err = sys.stdout, sys.stderr
        sys.stdout = _QWriter(q, "normal")
        sys.stderr = _QWriter(q, "yellow")
        html_path = None
        try:
            # ── Input ──────────────────────────────────────────────────────
            if self.v_mode.get() == "paste":
                raw   = self.seq_text.get("1.0", tk.END).strip()
                pairs = pla.parse_fasta(raw) if ">" in raw else [("Input Sequence", raw)]
            else:
                raw   = Path(self.v_fasta.get()).read_text(encoding="utf-8", errors="replace")
                stem  = Path(self.v_fasta.get()).stem
                pairs = pla.parse_fasta(raw) if raw.lstrip().startswith(">") else [(stem, raw)]

            if not pairs:
                q.put(("red", "  No sequences found.\n"))
                self.root.after(0, self._fail, "No sequences found."); return

            # ── PDB ────────────────────────────────────────────────────────
            pdb_text = None
            pf = self.v_pdb.get().strip()
            if pf and Path(pf).exists():
                pdb_text = Path(pf).read_text(encoding="utf-8", errors="replace")
                q.put(("cyan", f"  PDB loaded: {Path(pf).name}\n"))
            self._3d_pdb_text = pdb_text   # stash for 3D tab

            chain_id = self.v_chain.get().strip() or None
            run_hos  = True   # always

            # ── Analysis ───────────────────────────────────────────────────
            results = pla.process_sequences(
                pairs, pdb_text=pdb_text, chain_id=chain_id,
                run_hos=run_hos, verbose=True)

            for entry in results:
                pla.print_console_report(
                    entry["name"], entry["seq"], entry["findings"], entry["summary"],
                    ss=entry.get("ss"), rsa=entry.get("rsa"), use_color=False)

            # ── HTML ───────────────────────────────────────────────────────
            if not self.v_nohtml.get():
                out_dir = Path(self.v_outdir.get())
                out_dir.mkdir(parents=True, exist_ok=True)
                stem      = self.v_outname.get().strip() or "protein_liabilities"
                suffix    = "_pdb" if pdb_text else "_hos"
                html_path = out_dir / f"{stem}{suffix}.html"
                title     = ("Protein Sequence Liability Analysis + PDB Structure"
                             if pdb_text else
                             "Protein Sequence Liability Analysis + HOS Prediction")
                html_content = pla.build_html_report(
                    results, title=title,
                    cdr_scheme=self.v_cdr_scheme.get())
                html_path.write_text(html_content, encoding="utf-8")
                q.put(("green", f"\n  Report saved -> {html_path}\n"))

            self._last_pdb_used = (pdb_text is not None)
            self.root.after(0, self._success, results, html_path)

        except Exception as exc:
            import traceback
            q.put(("red",  f"\n  Error: {exc}\n"))
            q.put(("dim",  traceback.format_exc() + "\n"))
            self.root.after(0, self._fail, str(exc))
        finally:
            sys.stdout, sys.stderr = old_out, old_err

    # ── 3D Structure tab ──────────────────────────────────────────────────────
    def _build_3d_tab(self, parent):
        """Inline 3D molecular viewer — pure tkinter, no external dependencies."""
        # ── Top status bar ────────────────────────────────────────────────────
        top = tk.Frame(parent, bg=ACCENT, padx=12, pady=6)
        top.pack(fill="x")
        tk.Label(top, text="3D Structure", bg=ACCENT, fg=WHITE,
                 font=("Segoe UI", 11, "bold")).pack(side="left")
        self._3d_status_var = tk.StringVar(value="Run an analysis with a PDB file to populate.")
        tk.Label(top, textvariable=self._3d_status_var,
                 bg=ACCENT, fg="#aabbff", font=F_SMALL).pack(side="right", padx=8)

        # ── Body: controls | canvas ───────────────────────────────────────────
        body = tk.Frame(parent, bg=DARK)
        body.pack(fill="both", expand=True)

        # Controls panel (left column)
        ctrl = tk.Frame(body, bg=DARK, width=200)
        ctrl.pack(side="left", fill="y")
        ctrl.pack_propagate(False)

        def _ch(txt, bold=False):
            tk.Label(ctrl, text=txt, bg=DARK,
                     fg="#aabbff" if bold else "#ccccdd",
                     font=("Segoe UI", 9, "bold") if bold else ("Segoe UI", 9),
                     anchor="w").pack(fill="x", padx=12, pady=(8 if bold else 1, 0))

        def _sep():
            tk.Frame(ctrl, bg="#2a3060", height=1).pack(fill="x", padx=8, pady=5)

        _ch("COLOR MODE", bold=True)
        self._3d_color_mode = tk.StringVar(value="risk")
        for val, lbl in [("risk",  "By Risk"),
                         ("chain", "By Chain"),
                         ("aa",    "By AA Type"),
                         ("ss",    "By Sec. Structure")]:
            tk.Radiobutton(
                ctrl, text=lbl, variable=self._3d_color_mode, value=val,
                bg=DARK, fg="#ddddee", selectcolor=ACCENT2,
                activebackground=DARK, activeforeground=WHITE,
                font=("Segoe UI", 9), anchor="w",
                command=self._on_3d_color_change,
            ).pack(fill="x", padx=16, pady=1)

        # Dynamic legend (populated by _update_3d_legend when mode is aa or ss)
        self._3d_legend_frame = tk.Frame(ctrl, bg=DARK)
        self._3d_legend_frame.pack(fill="x")

        _sep()
        _ch("RISK FILTER", bold=True)
        self._3d_risk_vars: dict = {}
        for risk, col in [("high",   "#e74c3c"), ("medium", "#f39c12"),
                          ("low",    "#3498db"), ("info",   "#95a5a6")]:
            v = tk.BooleanVar(value=True)
            self._3d_risk_vars[risk] = v
            tk.Checkbutton(
                ctrl, text=risk.capitalize(), variable=v,
                bg=DARK, fg=col, selectcolor=ACCENT,
                activebackground=DARK, activeforeground=col,
                font=("Segoe UI", 9, "bold"), anchor="w",
                command=self._update_risk_filter,
            ).pack(fill="x", padx=16, pady=1)

        _sep()
        _ch("DISPLAY", bold=True)
        self._3d_labels_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            ctrl, text="Show labels", variable=self._3d_labels_var,
            bg=DARK, fg="#ccccdd", selectcolor=ACCENT,
            activebackground=DARK, activeforeground=WHITE,
            font=("Segoe UI", 9), anchor="w",
            command=lambda: self._mol_canvas.set_show_labels(
                self._3d_labels_var.get()),
        ).pack(fill="x", padx=16, pady=1)

        _sep()
        _reset_btn = tk.Label(
            ctrl, text="Reset View", bg=ACCENT, fg=WHITE,
            font=("Segoe UI", 9, "bold"), padx=8, pady=5, cursor="hand2",
        )
        _reset_btn.pack(fill="x", padx=12, pady=4)
        _reset_btn.bind("<Button-1>", lambda _: self._mol_canvas.reset_view())
        _reset_btn.bind("<Enter>",    lambda _: _reset_btn.configure(bg=ACCENT2))
        _reset_btn.bind("<Leave>",    lambda _: _reset_btn.configure(bg=ACCENT))

        _sep()
        _ch("INTERACTION", bold=True)
        tk.Label(ctrl, text="Drag  → rotate\nScroll → zoom",
                 bg=DARK, fg="#667799", font=("Segoe UI", 8),
                 justify="left").pack(anchor="w", padx=16, pady=2)

        _sep()
        _ch("RISK LEGEND", bold=True)
        for risk, col in [("High",   "#e74c3c"), ("Medium", "#f39c12"),
                          ("Low",    "#3498db"), ("Info",   "#95a5a6")]:
            row = tk.Frame(ctrl, bg=DARK)
            row.pack(fill="x", padx=16, pady=1)
            dot = tk.Canvas(row, width=14, height=14, bg=DARK, highlightthickness=0)
            dot.pack(side="left")
            dot.create_oval(2, 2, 12, 12, fill=col, outline="")
            tk.Label(row, text=risk, bg=DARK, fg=col,
                     font=("Segoe UI", 9, "bold")).pack(side="left", padx=4)

        # Vertical separator
        tk.Frame(body, bg="#2a3060", width=1).pack(side="left", fill="y")

        # 3D canvas (fills remaining space)
        self._mol_canvas = Mol3DCanvas(body)
        self._mol_canvas.pack(side="left", fill="both", expand=True)

    def _populate_3d_tab(self, results, pdb_text):
        self._3d_results  = results
        self._3d_pdb_text = pdb_text
        total = sum(len(e["findings"]) for e in results)

        if not pdb_text:
            self._3d_status_var.set(
                f"No PDB provided — {total} liabilities detected (sequence only)")
            self._mol_canvas.clear()
            return

        # ── Parse Cα atoms from PDB ───────────────────────────────────────────
        _aa3 = {
            "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q",
            "GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K",
            "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
            "TYR":"Y","VAL":"V",
        }
        ca_by_key: dict = {}  # (chain, resnum) -> [x,y,z]
        aa_by_key: dict = {}  # (chain, resnum) -> 1-letter
        for line in pdb_text.splitlines():
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                chain  = line[21].strip() or "A"
                resnum = int(line[22:26].strip())
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except (ValueError, IndexError):
                continue
            key = (chain, resnum)
            ca_by_key[key] = [x, y, z]
            aa_by_key[key] = _aa3.get(line[17:20].strip(), "G")

        if not ca_by_key:
            self._3d_status_var.set("PDB parsed — no Cα atoms found.")
            self._mol_canvas.clear()
            return

        # Center coordinates
        xs = [v[0] for v in ca_by_key.values()]
        ys = [v[1] for v in ca_by_key.values()]
        zs = [v[2] for v in ca_by_key.values()]
        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        cz = sum(zs) / len(zs)
        for key in ca_by_key:
            v = ca_by_key[key]
            ca_by_key[key] = [v[0] - cx, v[1] - cy, v[2] - cz]

        # ── Build liability and CDR maps keyed by (chain, resnum) ─────────────
        liability_map: dict = {}  # (chain,resnum) -> {risk, label, size}
        cdr_name_map:  dict = {}  # (chain,resnum) -> cdr_name
        ss_map:        dict = {}  # (chain,resnum) -> ss_char

        cdr_scheme = self.v_cdr_scheme.get()
        _CDR_RISK_SIZE = {"high": 10, "medium": 8, "low": 6, "info": 5}

        for e in results:
            s2p = e.get("seq_to_pdb", {})
            seq = e.get("sequence", "")
            ss  = e.get("ss", "")

            # Secondary structure per residue
            for pos0, (ch, rn) in s2p.items():
                if pos0 < len(ss):
                    ss_map[(ch, rn)] = ss[pos0]

            # Findings
            for f in e["findings"]:
                pos0 = f.get("pos0")
                if pos0 is None:
                    pos0 = f.get("start", 1) - 1
                if pos0 not in s2p:
                    continue
                key  = s2p[pos0]
                risk = f.get("risk", "info")
                existing = liability_map.get(key)
                if existing is None or _risk_rank(risk) < _risk_rank(existing["risk"]):
                    raw = re.sub(r'\s*[–—-]\s*(High|Medium|Low|Info)\s*Risk\s*', ' ',
                                 f.get("label", "")).strip()
                    short = raw[:15] + "…" if len(raw) > 16 else raw
                    liability_map[key] = {
                        "risk":  risk,
                        "label": f"{pos0+1} {short}",
                        "size":  _CDR_RISK_SIZE.get(risk, 6),
                    }

            # CDR regions
            if hasattr(pla, "annotate_cdrs"):
                for cs, ce, cname in pla.annotate_cdrs(seq, scheme=cdr_scheme):
                    for p in range(cs, ce):
                        if p in s2p:
                            cdr_name_map[s2p[p]] = cname

        # ── Assemble ordered atom + bond lists ────────────────────────────────
        _CDR_COLORS_3D = Mol3DCanvas._CDR_COLORS_3D
        _CHAIN_COLORS  = Mol3DCanvas._CHAIN_COLORS
        _AA_COLORS     = Mol3DCanvas._AA_COLORS

        chains: dict = {}
        for (chain, resnum) in ca_by_key:
            chains.setdefault(chain, []).append(resnum)

        atoms: list = []
        bonds: list = []
        for chain in sorted(chains):
            ci = ord(chain) % len(_CHAIN_COLORS)
            prev_idx = None
            for resnum in sorted(chains[chain]):
                key = (chain, resnum)
                hl  = liability_map.get(key)
                cdr = cdr_name_map.get(key, "")
                ss_char = ss_map.get(key, "")
                aa = aa_by_key.get(key, "G")
                atom = {
                    "xyz":          ca_by_key[key],
                    "chain":        chain,
                    "resnum":       resnum,
                    "aa":           aa,
                    "color_chain":  _CHAIN_COLORS[ci],
                    "color_cdr":    _CDR_COLORS_3D.get(cdr, "") if cdr else "",
                    "color_aa":     _AA_COLORS.get(aa, "#aaaaaa"),
                    "ss_char":      ss_char,
                    "is_liability": hl is not None,
                    "hl_risk":      hl["risk"]  if hl else "",
                    "hl_label":     hl["label"] if hl else "",
                    "hl_size":      hl["size"]  if hl else 3,
                    "cdr_name":     cdr,
                }
                idx = len(atoms)
                atoms.append(atom)
                if prev_idx is not None:
                    bonds.append((prev_idx, idx))
                prev_idx = idx

        mapped = sum(1 for e in results for f in e["findings"]
                     if (f.get("pos0") or f.get("start", 1) - 1) in e.get("seq_to_pdb", {}))
        names  = ", ".join(e["name"] for e in results)
        self._3d_status_var.set(
            f"{len(atoms)} Cα  ·  {mapped}/{total} liabilities mapped  ·  {names}")
        self._mol_canvas.load(atoms, bonds)

    def _on_3d_color_change(self):
        self._mol_canvas.set_color_mode(self._3d_color_mode.get())
        self._update_3d_legend()

    def _update_3d_legend(self):
        for w in self._3d_legend_frame.winfo_children():
            w.destroy()
        mode = self._3d_color_mode.get()

        def _legend_row(parent, label, col):
            row = tk.Frame(parent, bg=DARK)
            row.pack(fill="x", padx=16, pady=1)
            dot = tk.Canvas(row, width=14, height=14, bg=DARK, highlightthickness=0)
            dot.pack(side="left")
            dot.create_oval(2, 2, 12, 12, fill=col, outline="")
            tk.Label(row, text=label, bg=DARK, fg="#ccccdd",
                     font=("Segoe UI", 8)).pack(side="left", padx=4)

        if mode == "aa":
            tk.Label(self._3d_legend_frame, text="AA TYPE",
                     bg=DARK, fg="#aabbff", font=("Segoe UI", 9, "bold"),
                     anchor="w").pack(fill="x", padx=12, pady=(8, 2))
            for lbl, col in [
                ("Hydrophobic",  "#e8a87c"),
                ("Aromatic",     "#d4956a"),
                ("Polar",        "#88d8a3"),
                ("Cysteine",     "#f5e642"),
                ("Acidic (D,E)", "#e07070"),
                ("Basic (K,R)",  "#70a0e0"),
                ("His (H)",      "#9b80d0"),
                ("Gly (G)",      "#cccccc"),
            ]:
                _legend_row(self._3d_legend_frame, lbl, col)

        elif mode == "ss":
            tk.Label(self._3d_legend_frame, text="SEC. STRUCTURE",
                     bg=DARK, fg="#aabbff", font=("Segoe UI", 9, "bold"),
                     anchor="w").pack(fill="x", padx=12, pady=(8, 2))
            for lbl, col in [
                ("Helix",  "#e74c3c"),
                ("Sheet",  "#3498db"),
                ("Coil",   "#666688"),
            ]:
                _legend_row(self._3d_legend_frame, lbl, col)

    def _update_risk_filter(self):
        visible = {r for r, v in self._3d_risk_vars.items() if v.get()}
        self._mol_canvas.set_risk_filter(visible)

    def _success(self, results, html_path):
        self._running   = False
        self._html_path = str(html_path) if html_path else None
        self.progress.stop(); self.progress.pack_forget()
        self.run_btn.configure(state="normal", text="Run Analysis", bg=ACCENT)

        total = sum(len(e["findings"]) for e in results)
        high  = sum(len(e["summary"]["by_risk"].get("high", [])) for e in results)
        self.v_status.set(
            f"Done - {len(results)} chain(s)  |  {total} liabilities  |  {high} high-risk")

        if html_path:
            self.btn_html.configure(state="normal", cursor="hand2")
        self.btn_folder.configure(state="normal", cursor="hand2")

        # Populate in-app report and 3D tab
        self._populate_report(results)
        self._populate_3d_tab(results, self._3d_pdb_text)
        self.nb.select(2)   # -> Report tab

        # Log summary
        self._log(f"\n  {chr(61)*58}\n", "dim")
        for e in results:
            h = len(e["summary"]["by_risk"].get("high",   []))
            m = len(e["summary"]["by_risk"].get("medium", []))
            l = len(e["summary"]["by_risk"].get("low",    []))
            i = len(e["summary"]["by_risk"].get("info",   []))
            self._log(f"  {e['name'][:55]}\n"
                      f"    High: {h}  ·  Medium: {m}  ·  Low: {l}  ·  Info: {i}\n",
                      "red" if h > 0 else "green")
        if html_path:
            self._log("\n  Switch to the Report tab - or click Open HTML Report.\n", "cyan")

    def _fail(self, msg):
        self._running = False
        self.progress.stop(); self.progress.pack_forget()
        self.run_btn.configure(state="normal", text="Run Analysis", bg=ACCENT)
        self.v_status.set(f"Failed: {msg[:100]}")
        messagebox.showerror("Analysis Error", f"Analysis failed:\n\n{msg}")


# ══════════════════════════════════════════════════════════════════════════════
def main():
    root = tk.Tk()
    try:
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass
    App(root)
    root.mainloop()

if __name__ == "__main__":
    main()
