#!/usr/bin/env python3
"""
Protein Liability Analyzer v2 — Desktop GUI
============================================
Graphical front-end for protein_liability_analyzer_v2.py.

Run:  python protein_liability_gui.py

Requires protein_liability_analyzer_v2.py in the same folder.
All computation is local — no internet required.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import sys
import os
import queue
import webbrowser
import platform
import subprocess
from pathlib import Path
from datetime import datetime

# ── Locate the analysis module ────────────────────────────────────────────────
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


# ══════════════════════════════════════════════════════════════════════════════
class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Protein Liability Analyzer v2 — Letarte Scientific")
        self.root.configure(bg=LIGHT)
        self.root.minsize(820, 640)
        self._centre(960, 800)

        # State
        self.v_mode     = tk.StringVar(value="paste")
        self.v_fasta    = tk.StringVar()
        self.v_pdb      = tk.StringVar()
        self.v_chain    = tk.StringVar(value="")
        self.v_outdir   = tk.StringVar(value=str(_HERE))
        self.v_outname  = tk.StringVar(value="protein_liabilities")
        self.v_hos      = tk.BooleanVar(value=False)
        self.v_nohtml   = tk.BooleanVar(value=False)
        self.v_seqlen   = tk.StringVar(value="")
        self.v_status   = tk.StringVar(value="Ready")

        self._html_path = None
        self._running   = False
        self._log_q: queue.Queue = queue.Queue()

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
        s.configure("TFrame",      background=LIGHT)
        s.configure("TLabel",      background=LIGHT, foreground=DARK, font=F_BODY)
        s.configure("TEntry",      fieldbackground=WHITE, foreground=DARK,
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

    # ── Header bar ────────────────────────────────────────────────────────────
    def _build_header(self):
        hdr = tk.Frame(self.root, bg=DARK, height=66)
        hdr.pack(fill="x")
        hdr.pack_propagate(False)
        inner = tk.Frame(hdr, bg=DARK)
        inner.pack(fill="both", expand=True, padx=20, pady=10)
        tk.Label(inner, text="🔬  Protein Liability Analyzer  v2",
                 bg=DARK, fg=WHITE, font=F_TITLE).pack(anchor="w")
        tk.Label(inner,
                 text="Higher-Order Structure Extension  ·  Letarte Scientific Consulting  ·  "
                      "All computation local — sequences never transmitted",
                 bg=DARK, fg="#7788aa", font=("Segoe UI", 8)).pack(anchor="w")

    # ── Notebook ──────────────────────────────────────────────────────────────
    def _build_notebook(self):
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill="both", expand=True, padx=10, pady=(6, 2))
        t1 = ttk.Frame(self.nb)
        t2 = ttk.Frame(self.nb)
        self.nb.add(t1, text="  ⚙  Analysis Setup  ")
        self.nb.add(t2, text="  📊  Results & Log  ")
        self._build_setup(t1)
        self._build_results(t2)

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
        # Scrollable canvas
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

        self._card(inner, "📥   Input",                         self._input_card)
        self._card(inner, "🏗️   Higher-Order Structure Analysis (optional)", self._hos_card)
        self._card(inner, "💾   Output",                        self._output_card)

        run_wrap = tk.Frame(inner, bg=LIGHT)
        run_wrap.pack(fill="x", padx=14, pady=(10,18))

        self.progress = ttk.Progressbar(run_wrap, mode="indeterminate")
        # hidden until running

        self.run_btn = tk.Button(
            run_wrap, text="▶   Run Analysis", command=self._start,
            bg=ACCENT, fg=WHITE, font=("Segoe UI", 13, "bold"),
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

        # Paste area
        self.paste_frame = tk.Frame(p, bg=WHITE)
        hdr_row = tk.Frame(self.paste_frame, bg=WHITE)
        hdr_row.pack(fill="x", pady=(0,4))
        tk.Label(hdr_row, text="Sequence — FASTA or plain amino-acid string:",
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

        # File area
        self.file_frame = tk.Frame(p, bg=WHITE)
        tk.Label(self.file_frame, text="FASTA file:",
                 bg=WHITE, fg=DARK, font=F_BODY, width=12, anchor="w").pack(side="left")
        tk.Entry(self.file_frame, textvariable=self.v_fasta, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left", fill="x", expand=True, padx=(0,8))
        self._flatbtn(self.file_frame, "Browse…", self._pick_fasta).pack(side="left")

        self._toggle_mode()

    def _hos_card(self, p):
        tk.Checkbutton(p,
                       text="Run HOS analysis  (Chou-Fasman secondary structure + sequence-based RSA estimation)",
                       variable=self.v_hos, bg=WHITE, fg=DARK, font=F_BOLD,
                       activebackground=WHITE).pack(anchor="w")
        tk.Label(p, text="      ↳  All prediction runs locally from sequence — no internet, "
                          "no external tools.  Enabled automatically when a PDB file is loaded.",
                 bg=WHITE, fg=GREY, font=F_SMALL, wraplength=760, justify="left").pack(anchor="w", pady=(0,10))

        pdb_row = tk.Frame(p, bg=WHITE)
        pdb_row.pack(fill="x", pady=(0,6))
        tk.Label(pdb_row, text="PDB file\n(optional):", bg=WHITE, fg=DARK, font=F_BODY,
                 width=13, justify="left", anchor="w").pack(side="left")
        tk.Entry(pdb_row, textvariable=self.v_pdb, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left", fill="x", expand=True, padx=(0,8))
        self._flatbtn(pdb_row, "Browse…", self._pick_pdb).pack(side="left")

        ch_row = tk.Frame(p, bg=WHITE)
        ch_row.pack(fill="x", pady=(0,8))
        tk.Label(ch_row, text="Chain ID:", bg=WHITE, fg=DARK, font=F_BODY,
                 width=13, anchor="w").pack(side="left")
        tk.Entry(ch_row, textvariable=self.v_chain, width=6, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left")
        tk.Label(ch_row, text="  Leave blank for automatic selection",
                 bg=WHITE, fg=GREY, font=F_SMALL).pack(side="left")

        info = tk.Frame(p, bg="#e8f4fd")
        info.pack(fill="x")
        tk.Label(info,
                 text="ℹ  With a PDB file: secondary structure comes from HELIX/SHEET records; "
                      "solvent accessibility is computed via Shrake-Rupley SASA (all-atom, 1.4 Å probe) — "
                      "entirely local, no BioPython required.",
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
        self._flatbtn(d_row, "Browse…", self._pick_outdir).pack(side="left")

        n_row = tk.Frame(p, bg=WHITE)
        n_row.pack(fill="x", pady=(0,4))
        tk.Label(n_row, text="Report name:", bg=WHITE, fg=DARK, font=F_BODY,
                 width=13, anchor="w").pack(side="left")
        tk.Entry(n_row, textvariable=self.v_outname, font=F_BODY,
                 relief="flat", bd=0, highlightthickness=1,
                 highlightbackground=BORDER).pack(side="left", fill="x", expand=True)
        tk.Label(n_row, text=".html", bg=WHITE, fg=GREY, font=F_BODY).pack(side="left", padx=4)

        tk.Checkbutton(p, text="Console output only — skip HTML report generation",
                       variable=self.v_nohtml, bg=WHITE, fg=GREY, font=F_SMALL,
                       activebackground=WHITE).pack(anchor="w", pady=(6,0))

    # ── Results / Log tab ─────────────────────────────────────────────────────
    def _build_results(self, parent):
        toolbar = tk.Frame(parent, bg=LIGHT, pady=5)
        toolbar.pack(fill="x", padx=10)
        tk.Label(toolbar, text="Analysis Log", bg=LIGHT, fg=DARK, font=F_BOLD).pack(side="left")

        right = tk.Frame(toolbar, bg=LIGHT)
        right.pack(side="right")
        self.btn_html = self._flatbtn(right, "🌐  Open HTML Report", self._open_html,
                                      bg=ACCENT, fg=WHITE, state="disabled")
        self.btn_html.pack(side="left", padx=(0,6))
        self.btn_folder = self._flatbtn(right, "📁  Open Folder", self._open_folder,
                                        state="disabled")
        self.btn_folder.pack(side="left", padx=(0,14))
        self._flatbtn(right, "Copy log", self._copy_log).pack(side="left", padx=(0,4))
        self._flatbtn(right, "Clear",    self._clear_log).pack(side="left")

        # Log widget with dark terminal look
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
            "\n  Welcome to Protein Liability Analyzer v2\n"
            "  ─────────────────────────────────────────────────────────\n"
            "  Set up your analysis in the Setup tab, then click Run Analysis.\n"
            "  Results and progress will stream here in real time.\n\n",
            "dim"
        )

    # ── Status bar ────────────────────────────────────────────────────────────
    def _build_statusbar(self):
        bar = tk.Frame(self.root, bg=DARK, height=22)
        bar.pack(fill="x", side="bottom")
        bar.pack_propagate(False)
        tk.Label(bar, textvariable=self.v_status, bg=DARK, fg="#7788aa",
                 font=("Segoe UI", 8), padx=12).pack(side="left", pady=2)
        tk.Label(bar, text="🔒  Local only — sequences never transmitted",
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
            self.v_pdb.set(p); self.v_hos.set(True)

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
        if any(x in msg for x in ("❌", "Error:", "error:")):          return "red"
        if any(x in msg for x in ("⚠", "[Warning]", "Warning")):       return "yellow"
        if any(x in msg for x in ("✓", "done", "saved →", "complete", "📄")): return "green"
        if any(x in msg for x in ("▸", "Computing", "Predicting",
                                   "Estimating", "Parsing", "Analysing")): return "cyan"
        if "─" * 5 in msg:                                               return "header"
        if msg.strip().startswith("Pos "):                               return "bold"
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

        # Validate
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

        self.nb.select(1)
        self._log(f"\n{'─'*62}\n  🔬  Started: {datetime.now():%Y-%m-%d %H:%M:%S}\n", "cyan")

        self.run_btn.configure(state="disabled", text="⏳  Running…", bg="#555")
        self.progress.pack(fill="x", pady=(0,8), before=self.run_btn)
        self.progress.start(12)
        self._running = True
        self._html_path = None
        self.btn_html.configure(state="disabled", cursor="arrow")
        self.btn_folder.configure(state="disabled", cursor="arrow")
        self.v_status.set("Running analysis…")

        threading.Thread(target=self._worker, daemon=True).start()

    def _worker(self):
        q = self._log_q
        old_out, old_err = sys.stdout, sys.stderr
        sys.stdout = _QWriter(q, "normal")
        sys.stderr = _QWriter(q, "yellow")
        html_path  = None
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
                q.put(("red", "  ❌  No sequences found.\n"))
                self.root.after(0, self._fail, "No sequences found."); return

            # ── PDB ────────────────────────────────────────────────────────
            pdb_text = None
            pf = self.v_pdb.get().strip()
            if pf and Path(pf).exists():
                pdb_text = Path(pf).read_text(encoding="utf-8", errors="replace")
                q.put(("cyan", f"  📂  PDB loaded: {Path(pf).name}\n"))

            chain_id = self.v_chain.get().strip() or None
            run_hos  = self.v_hos.get() or (pdb_text is not None)

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
                stem    = self.v_outname.get().strip() or "protein_liabilities"
                suffix  = "_hos" if run_hos else ""
                html_path = out_dir / f"{stem}{suffix}.html"
                title   = "Protein Sequence Liability Analysis"
                if run_hos:
                    title += " + PDB Structure" if pdb_text else " + HOS Prediction"
                pla.build_html_report(results, title=title).encode  # call for side-effect check
                html_content = pla.build_html_report(results, title=title)
                html_path.write_text(html_content, encoding="utf-8")
                q.put(("green", f"\n  📄  Report saved → {html_path}\n"))

            self.root.after(0, self._success, results, html_path)

        except Exception as exc:
            import traceback
            q.put(("red", f"\n  ❌  Error: {exc}\n"))
            q.put(("dim", traceback.format_exc() + "\n"))
            self.root.after(0, self._fail, str(exc))
        finally:
            sys.stdout, sys.stderr = old_out, old_err

    def _success(self, results, html_path):
        self._running   = False
        self._html_path = str(html_path) if html_path else None
        self.progress.stop(); self.progress.pack_forget()
        self.run_btn.configure(state="normal", text="▶   Run Analysis", bg=ACCENT)
        total = sum(len(e["findings"]) for e in results)
        high  = sum(len(e["summary"]["by_risk"].get("high", [])) for e in results)
        self.v_status.set(f"Done — {len(results)} chain(s) · {total} liabilities · {high} high-risk")
        if html_path:
            self.btn_html.configure(state="normal", cursor="hand2")
        self.btn_folder.configure(state="normal", cursor="hand2")
        # Summary in log
        self._log(f"\n  {'═'*58}\n", "dim")
        for e in results:
            h = len(e["summary"]["by_risk"].get("high",   []))
            m = len(e["summary"]["by_risk"].get("medium", []))
            l = len(e["summary"]["by_risk"].get("low",    []))
            i = len(e["summary"]["by_risk"].get("info",   []))
            self._log(f"  {e['name'][:55]}\n"
                      f"    High: {h}  ·  Medium: {m}  ·  Low: {l}  ·  Info: {i}\n",
                      "red" if h > 0 else "green")
        if html_path:
            self._log("\n  👉  Click 'Open HTML Report' to view the full interactive report.\n", "cyan")

    def _fail(self, msg):
        self._running = False
        self.progress.stop(); self.progress.pack_forget()
        self.run_btn.configure(state="normal", text="▶   Run Analysis", bg=ACCENT)
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
