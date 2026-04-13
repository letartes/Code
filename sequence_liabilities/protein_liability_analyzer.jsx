/**
 * Protein Sequence Liability Analyzer — React Component
 * ======================================================
 * Drop-in React component (no external CSS required).
 * For use in a build pipeline: import React and ReactDOM separately.
 * For standalone browser use, see protein_liability_analyzer.html.
 *
 * Detected liabilities:
 *   Deamidation · Oxidation · Isomerization · Glycosylation
 *   Pyroglutamate · Truncation/Clipping
 *
 * CDR annotation (VH / VL) is performed automatically when antibody
 * variable-domain sequences are detected (Chothia-like approximation).
 *
 * Letarte Scientific Consulting
 */

// ─── CONSTANTS ───────────────────────────────────────────────────────────────

const RISK_ORDER = { high: 0, medium: 1, low: 2, info: 3 };

const LIABILITIES = {
  deamid_high: {
    label: "Deamidation – High Risk",
    display_label: "Deamidation (NG / NS)",
    short: "Deam-H", category: "Deamidation", risk: "high",
    color_text: "#c0392b", color_bg: "#fde8e8", color_border: "#e74c3c",
    description: "NG and NS motifs. Asn–Gly is the most susceptible sequence to deamidation, forming succinimide intermediates rapidly at physiological pH.",
    span_length: 2,
  },
  deamid_medium: {
    label: "Deamidation – Medium Risk",
    display_label: "Deamidation (NA / NT / NQ / NR / QG / QS)",
    short: "Deam-M", category: "Deamidation", risk: "medium",
    color_text: "#d35400", color_bg: "#fef0e6", color_border: "#e67e22",
    description: "NA, NT, NQ, NR, NH motifs and Gln-based equivalents (QG, QS). Glutamine can deamidate, though more slowly than asparagine.",
    span_length: 2,
  },
  oxid_high: {
    label: "Oxidation – High Risk",
    display_label: "Oxidation (Met, Trp)",
    short: "Ox-H", category: "Oxidation", risk: "high",
    color_text: "#6c3483", color_bg: "#f4ecf7", color_border: "#9b59b6",
    description: "Met (M) and Trp (W). Methionine sulfoxide formation is a primary oxidative liability; tryptophan oxidation is common under photo-oxidative stress.",
    span_length: 1,
  },
  oxid_medium: {
    label: "Oxidation – Medium Risk",
    display_label: "Oxidation (His, Cys)",
    short: "Ox-M", category: "Oxidation", risk: "medium",
    color_text: "#1a5276", color_bg: "#eaf4fb", color_border: "#2980b9",
    description: "His (H) and Cys (C). Histidine oxidised to 2-oxo-histidine under metal-catalysed conditions; free cysteines prone to disulfide scrambling.",
    span_length: 1,
  },
  oxid_low: {
    label: "Oxidation – Low Risk",
    display_label: "Oxidation (Phe, Tyr)",
    short: "Ox-L", category: "Oxidation", risk: "low",
    color_text: "#0e6655", color_bg: "#e8f8f5", color_border: "#1abc9c",
    description: "Phe (F) and Tyr (Y). Susceptible under severe photo-oxidative or radiolytic stress conditions.",
    span_length: 1,
  },
  isom_high: {
    label: "Isomerization – High Risk",
    display_label: "Isomerization (DG / DS)",
    short: "Isom-H", category: "Isomerization", risk: "high",
    color_text: "#7d6608", color_bg: "#fefde7", color_border: "#f1c40f",
    description: "DG and DS motifs. Asp–Gly is the most prone to succinimide-mediated racemization and isomerization.",
    span_length: 2,
  },
  isom_medium: {
    label: "Isomerization – Medium Risk",
    display_label: "Isomerization (DT / DA / DN / DD / DH)",
    short: "Isom-M", category: "Isomerization", risk: "medium",
    color_text: "#6e2f0e", color_bg: "#fdf0e7", color_border: "#e08040",
    description: "DT, DA, DN, DD, DH motifs. Moderate succinimide intermediate formation rate.",
    span_length: 2,
  },
  n_glycan: {
    label: "N-Glycosylation Sequon",
    display_label: "N-Glycan sequon (N-X-S/T, X≠P)",
    short: "N-Gly", category: "Glycosylation", risk: "medium",
    color_text: "#1d6d2e", color_bg: "#e8f8e8", color_border: "#27ae60",
    description: "N-X-S/T sequon (X ≠ Pro). Potential site for complex or high-mannose N-glycan addition in the ER/Golgi.",
    span_length: 3,
  },
  o_glycan: {
    label: "O-Glycosylation Site",
    display_label: "O-Glycan (S/T–P)",
    short: "O-Gly", category: "Glycosylation", risk: "low",
    color_text: "#117a65", color_bg: "#d1f2eb", color_border: "#1abc9c",
    description: "S/T–Pro motifs. Potential for mucin-type O-glycan addition.",
    span_length: 2,
  },
  pyroglu: {
    label: "Pyroglutamate Formation",
    display_label: "Pyroglutamate (N-term Q/E cyclisation)",
    short: "pyroGlu", category: "Pyroglutamate", risk: "medium",
    color_text: "#7e5109", color_bg: "#fef5e7", color_border: "#f0a500",
    description: "N-terminal Gln (–17 Da, loss of NH₃) or Glu (–18 Da, loss of H₂O) spontaneously cyclises to pyroglutamate.",
    span_length: 1,
  },
  trunc_high: {
    label: "Truncation – High Risk",
    display_label: "Truncation (Asp–Pro acid-labile bond)",
    short: "Trunc-H", category: "Truncation", risk: "high",
    color_text: "#922b21", color_bg: "#f9ebea", color_border: "#cb4335",
    description: "Asp–Pro (DP) peptide bond is unusually labile under mildly acidic conditions (pH 3–4), causing chain fragmentation.",
    span_length: 2,
  },
  trunc_medium: {
    label: "Truncation – Medium Risk",
    display_label: "Truncation (C-terminal Lys / N-terminal Met)",
    short: "Trunc-M", category: "Truncation", risk: "medium",
    color_text: "#7b241c", color_bg: "#fdedec", color_border: "#e74c3c",
    description: "C-terminal Lys clipping by carboxypeptidase B (CHO expression) or N-terminal Met removal by Met aminopeptidase.",
    span_length: 1,
  },
};

const CAT_ORDER = [
  "Deamidation", "Oxidation", "Isomerization",
  "Glycosylation", "Pyroglutamate", "Truncation",
];
const CAT_ICONS = {
  Deamidation: "🧪", Oxidation: "⚡", Isomerization: "🔄",
  Glycosylation: "🍬", Pyroglutamate: "🔵", Truncation: "✂️",
};

const RISK_COLORS = {
  high:   { bg: "#fde8e8", text: "#c0392b", border: "#e74c3c" },
  medium: { bg: "#fef0e6", text: "#d35400", border: "#e67e22" },
  low:    { bg: "#e8f8f5", text: "#0e6655", border: "#1abc9c" },
  info:   { bg: "#d6eaf8", text: "#154360", border: "#2e86c1" },
};

const THEME = {
  bgDark:   "#1a1a2e",
  bgMid:    "#16213e",
  bgAccent: "#0f3460",
  btnBlue:  "#0d2d6b",
  bgLight:  "#f0f2f5",
  bgWhite:  "#ffffff",
  bgPanel:  "#f8f9fa",
  fgDark:   "#1a1a2e",
  fgMuted:  "#6b7280",
  fgSubtle: "#9ba5b0",
  border:   "#e5e7eb",
};


// ─── ANALYSIS ENGINE ──────────────────────────────────────────────────────────

function parseFasta(text) {
  const sequences = [];
  let currentName = null;
  let currentSeq = [];
  for (const rawLine of text.split("\n")) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith(">")) {
      if (currentSeq.length > 0 && currentName !== null) {
        sequences.push([currentName, currentSeq.join("").toUpperCase()]);
      }
      currentName = line.slice(1).trim() || "Sequence";
      currentSeq = [];
    } else {
      currentSeq.push(line.replace(/[^A-Za-z]/g, ""));
    }
  }
  if (currentSeq.length > 0 && currentName !== null) {
    sequences.push([currentName, currentSeq.join("").toUpperCase()]);
  }
  return sequences;
}

function cleanSequence(seq) {
  return seq.replace(/[^A-Za-z]/g, "").toUpperCase();
}

function findCdrRegions(seq) {
  const s = seq.toUpperCase();
  const cdrs = [];

  // ── VH ───────────────────────────────────────────────────────────────────
  const fr4VhMatch = s.match(/WG[QKHR]G[TS]/);
  if (fr4VhMatch) {
    const fr4 = fr4VhMatch.index;
    let h3Found = null;

    // CDR-H3
    const win3 = s.slice(Math.max(0, fr4 - 30), fr4);
    const rc3 = win3.lastIndexOf("C");
    if (rc3 >= 0) {
      const h3s = Math.max(0, fr4 - 30) + rc3 + 1;
      const h3e = fr4 - 1;
      const len = h3e - h3s + 1;
      if (len >= 2 && len <= 28) { cdrs.push([h3s, h3e, "CDR-H3"]); h3Found = h3s; }
    }

    // FR2 anchor → CDR-H1, CDR-H2
    const fr2Seg = s.slice(0, Math.max(0, fr4 - 35));
    let fr2Match = null;
    for (const m of [...fr2Seg.matchAll(/W[VILMAF][RQHK]/g)]) fr2Match = m;
    if (fr2Match) {
      const fr2 = fr2Match.index;
      const win1 = s.slice(Math.max(0, fr2 - 22), fr2);
      const rc1 = win1.lastIndexOf("C");
      if (rc1 >= 0) {
        const cAbs = Math.max(0, fr2 - 22) + rc1;
        const h1s = cAbs + 3; const h1e = fr2 - 1;
        const len = h1e - h1s + 1;
        if (len >= 4 && len <= 15) cdrs.push([h1s, h1e, "CDR-H1"]);
      }
      let h2s = fr2 + 14;
      let h2e = h2s + 16;
      if (h3Found !== null) h2e = Math.min(h2e, h3Found - 29);
      const len2 = h2e - h2s + 1;
      if (len2 >= 5 && len2 <= 22) cdrs.push([h2s, h2e, "CDR-H2"]);
    }
    return cdrs.sort((a, b) => a[0] - b[0]);
  }

  // ── VL ───────────────────────────────────────────────────────────────────
  const fr4VlMatch = s.match(/FG[QNRST]G[TS]/);
  if (fr4VlMatch) {
    const fr4 = fr4VlMatch.index;
    let l3Found = null;

    // CDR-L3
    const win3 = s.slice(Math.max(0, fr4 - 18), fr4);
    const rc3 = win3.lastIndexOf("C");
    if (rc3 >= 0) {
      const l3s = Math.max(0, fr4 - 18) + rc3 + 1;
      const l3e = fr4 - 1;
      const len = l3e - l3s + 1;
      if (len >= 3 && len <= 14) { cdrs.push([l3s, l3e, "CDR-L3"]); l3Found = l3s; }
    }

    const fr2Seg = s.slice(0, Math.max(0, fr4 - 35));
    let fr2Match = null;
    for (const m of [...fr2Seg.matchAll(/W[YFH][LQP]/g)]) fr2Match = m;
    if (fr2Match) {
      const fr2 = fr2Match.index;
      const win1 = s.slice(Math.max(0, fr2 - 28), fr2);
      const rc1 = win1.lastIndexOf("C");
      if (rc1 >= 0) {
        const cAbs = Math.max(0, fr2 - 28) + rc1;
        const l1s = cAbs + 1; const l1e = fr2 - 1;
        const len = l1e - l1s + 1;
        if (len >= 6 && len <= 18) cdrs.push([l1s, l1e, "CDR-L1"]);
      }
      let l2s = fr2 + 13;
      let l2e = l2s + 6;
      if (l3Found !== null) l2e = Math.min(l2e, l3Found - 29);
      const len2 = l2e - l2s + 1;
      if (len2 >= 5 && len2 <= 12) cdrs.push([l2s, l2e, "CDR-L2"]);
    }
    return cdrs.sort((a, b) => a[0] - b[0]);
  }

  return [];
}

function cdrNameAt(pos0, cdrs) {
  for (const [start, end, name] of cdrs) {
    if (pos0 >= start && pos0 <= end) return name;
  }
  return null;
}

function findLiabilities(seq, cdrs = []) {
  const findings = [];
  const n = seq.length;

  function add(key, pos0, span, note = "") {
    const lib = LIABILITIES[key];
    findings.push({
      type: key, pos0, span, pos1: pos0 + 1,
      residues: seq.slice(pos0, pos0 + span),
      note, cdr: cdrNameAt(pos0, cdrs),
      label: lib.label, display_label: lib.display_label, short: lib.short,
      category: lib.category, risk: lib.risk,
      color_text: lib.color_text, color_bg: lib.color_bg, color_border: lib.color_border,
    });
  }

  // Deamidation
  const HIGH_DEAMID = new Set(["NG", "NS"]);
  const MED_DEAMID  = new Set(["NA", "NT", "NQ", "NR", "NH", "QG", "QS"]);
  for (let i = 0; i < n - 1; i++) {
    const pair = seq.slice(i, i + 2);
    if (HIGH_DEAMID.has(pair)) add("deamid_high", i, 2);
    else if (MED_DEAMID.has(pair)) add("deamid_medium", i, 2);
  }

  // Oxidation
  for (let i = 0; i < n; i++) {
    const aa = seq[i];
    if ("MW".includes(aa)) add("oxid_high", i, 1);
    else if ("HC".includes(aa)) add("oxid_medium", i, 1);
    else if ("FY".includes(aa)) add("oxid_low", i, 1);
  }

  // Isomerization
  const HIGH_ISOM = new Set(["DG", "DS"]);
  const MED_ISOM  = new Set(["DT", "DA", "DN", "DD", "DH"]);
  for (let i = 0; i < n - 1; i++) {
    const pair = seq.slice(i, i + 2);
    if (HIGH_ISOM.has(pair)) add("isom_high", i, 2);
    else if (MED_ISOM.has(pair)) add("isom_medium", i, 2);
  }

  // N-Glycosylation sequon
  for (let i = 0; i < n - 2; i++) {
    if (seq[i] === "N" && seq[i + 1] !== "P" && "ST".includes(seq[i + 2]))
      add("n_glycan", i, 3);
  }

  // O-Glycosylation
  for (let i = 0; i < n - 1; i++) {
    if ("ST".includes(seq[i]) && seq[i + 1] === "P") add("o_glycan", i, 2);
  }

  // Pyroglutamate
  if (n > 0 && "QE".includes(seq[0])) {
    const note = seq[0] === "Q"
      ? "N-term Gln → pyroglutamate (–17 Da, loss of NH₃)"
      : "N-term Glu → pyroglutamate (–18 Da, loss of H₂O)";
    add("pyroglu", 0, 1, note);
  }

  // Truncation: DP acid-labile
  for (let i = 0; i < n - 1; i++) {
    if (seq.slice(i, i + 2) === "DP") add("trunc_high", i, 2, "Acid-labile Asp–Pro bond");
  }

  // C-terminal Lys clipping
  if (n > 0 && seq[n - 1] === "K")
    add("trunc_medium", n - 1, 1, "C-terminal Lys clipping (CHO expression)");

  // N-terminal Met removal
  if (n > 0 && seq[0] === "M")
    add("trunc_medium", 0, 1, "N-terminal Met — Met aminopeptidase removal");

  return findings;
}

function summarize(findings, seqLen) {
  const byCategory = {}, byRisk = {};
  for (const f of findings) {
    if (!byCategory[f.category]) byCategory[f.category] = [];
    byCategory[f.category].push(f);
    if (!byRisk[f.risk]) byRisk[f.risk] = [];
    byRisk[f.risk].push(f);
  }
  return { total: findings.length, byCategory, byRisk, seqLen };
}

function processSequences(pairs) {
  return pairs.map(([name, raw]) => {
    const seq = cleanSequence(raw);
    const cdrs = findCdrRegions(seq);
    const findings = findLiabilities(seq, cdrs);
    const summary = summarize(findings, seq.length);
    return { name, seq, findings, summary, cdrs };
  });
}


// ─── UTILITY ──────────────────────────────────────────────────────────────────

function RiskBadge({ risk, small = false }) {
  const c = RISK_COLORS[risk] || RISK_COLORS.info;
  return (
    <span style={{
      display: "inline-block",
      padding: small ? "1px 6px" : "2px 8px",
      borderRadius: 20,
      fontSize: small ? "0.72em" : "0.8em",
      fontWeight: 700,
      letterSpacing: "0.4px",
      background: c.bg, color: c.text,
      border: `1px solid ${c.border}`,
      whiteSpace: "nowrap",
    }}>
      {risk.toUpperCase()}
    </span>
  );
}

function SortIcon({ col, sortCol, sortAsc }) {
  if (sortCol !== col) return <span style={{ color: "#ccc", marginLeft: 4 }}>⇅</span>;
  return <span style={{ marginLeft: 4 }}>{sortAsc ? "↑" : "↓"}</span>;
}


// ─── SUMMARY TAB ─────────────────────────────────────────────────────────────

function SummaryTab({ result }) {
  if (!result) return (
    <div style={{ display: "flex", alignItems: "center", justifyContent: "center",
                  height: "100%", color: THEME.fgMuted, fontSize: 16 }}>
      🔬&nbsp; Paste a sequence and click <strong style={{ margin: "0 6px" }}>ANALYZE</strong>
    </div>
  );

  const { findings, summary, cdrs } = result;
  const highCt = (summary.byRisk.high || []).length;

  return (
    <div style={{ padding: "16px 20px", overflowY: "auto", height: "100%" }}>
      {/* Alert banner */}
      <div style={{
        padding: "10px 16px", borderRadius: 8, marginBottom: 12, fontWeight: 600,
        fontSize: "0.9em",
        background: highCt ? "#fde8e8" : "#d6eaf8",
        color:      highCt ? "#c0392b" : "#154360",
      }}>
        {highCt
          ? `⚠️  ${highCt} high-risk liabilit${highCt === 1 ? "y" : "ies"} detected`
          : "✅  No high-risk liabilities detected"}
      </div>

      {/* CDR banner */}
      {cdrs.length > 0 && (
        <div style={{ background: "#e3f2fd", color: "#1565c0", padding: "8px 14px",
                      borderRadius: 6, marginBottom: 12, fontSize: "0.88em" }}>
          📐&nbsp;<strong>CDRs detected (Chothia approx.):</strong>&nbsp;
          {cdrs.map(([s, e, n]) => (
            <span key={n} style={{ marginRight: 12 }}>
              <strong>{n}</strong> ({s + 1}–{e + 1})
            </span>
          ))}
        </div>
      )}

      {/* Category cards */}
      <div style={{ display: "grid", gridTemplateColumns: "repeat(6, 1fr)", gap: 8, marginBottom: 14 }}>
        {CAT_ORDER.map(cat => {
          const items = summary.byCategory[cat] || [];
          const high   = items.filter(f => f.risk === "high").length;
          const medium = items.filter(f => f.risk === "medium").length;
          const low    = items.filter(f => f.risk === "low").length;
          return (
            <div key={cat} style={{
              background: THEME.bgWhite, border: `1px solid ${THEME.border}`,
              borderRadius: 8, padding: "10px 8px", textAlign: "center",
            }}>
              <div style={{ fontSize: 20 }}>{CAT_ICONS[cat]}</div>
              <div style={{ fontSize: 22, fontWeight: 700, color: THEME.fgDark, lineHeight: 1.2 }}>
                {items.length}
              </div>
              <div style={{ fontSize: "0.72em", fontWeight: 700, color: THEME.fgMuted, marginBottom: 4 }}>
                {cat}
              </div>
              <div style={{ display: "flex", flexWrap: "wrap", gap: 2, justifyContent: "center" }}>
                {high   > 0 && <RiskBadge risk="high"   small />}
                {medium > 0 && <RiskBadge risk="medium" small />}
                {low    > 0 && <RiskBadge risk="low"    small />}
              </div>
            </div>
          );
        })}
      </div>

      {/* Stats row */}
      <div style={{ fontSize: "0.88em", color: THEME.fgMuted, marginBottom: 16 }}>
        Length: <strong>{summary.seqLen}</strong> aa &nbsp;·&nbsp;
        Total: <strong>{summary.total}</strong> &nbsp;·&nbsp;
        High: <strong style={{ color: "#c0392b" }}>{(summary.byRisk.high||[]).length}</strong> &nbsp;·&nbsp;
        Med: <strong style={{ color: "#d35400" }}>{(summary.byRisk.medium||[]).length}</strong> &nbsp;·&nbsp;
        Low: <strong style={{ color: "#0e6655" }}>{(summary.byRisk.low||[]).length}</strong>
      </div>

      {/* Legend */}
      <div style={{ background: THEME.bgPanel, borderRadius: 8, padding: "12px 16px",
                    border: `1px solid ${THEME.border}` }}>
        <div style={{ fontWeight: 700, fontSize: "0.88em", marginBottom: 8, color: THEME.fgMuted }}>
          Liability Legend
        </div>
        <div style={{ display: "flex", flexWrap: "wrap", gap: 6 }}>
          {Object.entries(LIABILITIES).map(([key, lib]) => (
            <span key={key} style={{
              padding: "2px 7px", borderRadius: 4, fontSize: "0.78em", fontWeight: 700,
              background: lib.color_bg, color: lib.color_text,
              border: `1px solid ${lib.color_border}`, fontFamily: "monospace",
              cursor: "default", title: lib.description,
            }}>
              {lib.short}
            </span>
          ))}
        </div>
        <div style={{ fontSize: "0.76em", color: "#aaa", marginTop: 8 }}>
          CDR boundaries are approximate (Chothia-like anchor residues).
          Tryptic K/R sites are excluded from truncation reporting.
        </div>
      </div>
    </div>
  );
}


// ─── ANNOTATED SEQUENCE TAB ───────────────────────────────────────────────────

function AnnotatedSequenceTab({ results }) {
  const WRAP = 60;

  if (!results || results.length === 0) return (
    <div style={{ display: "flex", alignItems: "center", justifyContent: "center",
                  height: "100%", color: THEME.fgMuted, fontSize: 16 }}>
      Run an analysis to see the annotated sequence.
    </div>
  );

  return (
    <div style={{ padding: "12px 16px", overflowY: "auto", height: "100%" }}>
      {/* Color key */}
      <div style={{ background: THEME.bgPanel, padding: "6px 12px", borderRadius: 6,
                    marginBottom: 12, display: "flex", flexWrap: "wrap", gap: 6,
                    alignItems: "center", fontSize: "0.82em", color: THEME.fgMuted }}>
        <strong>Key:</strong>
        {Object.entries(LIABILITIES).map(([key, lib]) => (
          <span key={key} style={{
            padding: "1px 5px", borderRadius: 3, fontFamily: "monospace", fontWeight: 700,
            background: lib.color_bg, color: lib.color_text,
            border: `1px solid ${lib.color_border}`,
          }}>{lib.short}</span>
        ))}
        <span style={{ marginLeft: 8 }}>·</span>
        <span style={{ textDecoration: "underline", textDecorationThickness: 2,
                       fontFamily: "monospace", color: "#1565c0", fontWeight: 700 }}>
          CDR
        </span>
        <span style={{ color: THEME.fgMuted, fontStyle: "italic" }}>(underlined)</span>
      </div>

      {results.map(({ name, seq, findings, cdrs }, chainIdx) => {
        const n = seq.length;

        // per-position liability map (best = highest risk)
        const posF = Array.from({ length: n }, () => []);
        for (const f of findings) {
          for (let j = f.pos0; j < Math.min(f.pos0 + f.span, n); j++) posF[j].push(f);
        }
        const bestAt = i => {
          const fl = posF[i];
          if (!fl.length) return null;
          return fl.reduce((a, b) => RISK_ORDER[a.risk] <= RISK_ORDER[b.risk] ? a : b);
        };

        // CDR position set
        const cdrPos = {};
        for (const [cs, ce, cn] of cdrs) {
          for (let p = cs; p <= ce; p++) cdrPos[p] = cn;
        }

        // Build wrapped lines
        const lines = [];
        for (let start = 0; start < n; start += WRAP) {
          const end = Math.min(start + WRAP, n);
          const chars = [];
          for (let i = start; i < end; i++) {
            const aa = seq[i];
            const best = bestAt(i);
            const inCdr = cdrPos[i] !== undefined;
            let style = {
              fontFamily: "monospace", fontSize: 13, display: "inline-block",
              lineHeight: "22px",
            };
            let title = `Pos ${i + 1}`;
            if (best) {
              style = {
                ...style,
                background: best.color_bg, color: best.color_text,
                borderBottom: `2px solid ${best.color_border}`,
                fontWeight: 700,
                ...(inCdr ? { textDecoration: "underline", textDecorationThickness: 2 } : {}),
              };
              title = `Pos ${i + 1}: ${best.short}${inCdr ? ` [${cdrPos[i]}]` : ""}`;
              if (posF[i].length > 1)
                title += ` + ${posF[i].length - 1} more`;
            } else if (inCdr) {
              style = { ...style, textDecoration: "underline", textDecorationThickness: 2,
                        textUnderlineOffset: 2, color: "#1565c0" };
              title = `Pos ${i + 1}: ${cdrPos[i]}`;
            }
            chars.push(
              <span key={i} style={style} title={title}>{aa}</span>
            );
          }

          // CDR annotation for this line
          const cdrTags = cdrs
            .filter(([cs, ce]) => cs < end && ce >= start)
            .map(([cs, ce, cn]) => {
              const ls = Math.max(cs, start);
              const le = Math.min(ce, end - 1);
              return (
                <span key={cn} style={{ fontSize: "0.75em", color: "#1565c0",
                                        fontWeight: 700, marginLeft: 6 }}>
                  {cn}:{ls + 1}–{le + 1}
                </span>
              );
            });

          lines.push(
            <div key={start} style={{ lineHeight: "22px", marginBottom: 2 }}>
              <span style={{ color: THEME.fgSubtle, fontSize: 11,
                             fontFamily: "monospace", minWidth: 40,
                             display: "inline-block", userSelect: "none" }}>
                {String(start + 1).padStart(5)}
              </span>
              {" "}{chars}{" "}
              <span style={{ color: THEME.fgSubtle, fontSize: 11,
                             fontFamily: "monospace", userSelect: "none" }}>
                {end}
              </span>
              {cdrTags}
            </div>
          );
        }

        // Chain header
        const cdrHeader = cdrs.length > 0
          ? "  [" + cdrs.map(([s, e, n]) => `${n}:${s + 1}-${e + 1}`).join("  ") + "]"
          : "";

        return (
          <div key={chainIdx} style={{ marginBottom: 24 }}>
            <div style={{ fontFamily: "monospace", fontWeight: 700, fontSize: 13,
                          color: "#1565c0", marginBottom: 4 }}>
              ▶ {name}{cdrHeader}
            </div>
            <div style={{ background: THEME.bgPanel, borderRadius: 8, padding: "12px 14px",
                          border: `1px solid ${THEME.border}` }}>
              {lines}
            </div>
          </div>
        );
      })}
    </div>
  );
}


// ─── ALL FINDINGS TAB ────────────────────────────────────────────────────────

function AllFindingsTab({ results }) {
  const [riskFilter, setRiskFilter] = React.useState("All");
  const [catFilter,  setCatFilter]  = React.useState("All");
  const [cdrOnly,    setCdrOnly]    = React.useState(false);
  const [sortCol,    setSortCol]    = React.useState("pos");
  const [sortAsc,    setSortAsc]    = React.useState(true);

  const allRows = React.useMemo(() => {
    if (!results) return [];
    const rows = [];
    for (const { name, findings } of results) {
      for (const f of findings) {
        rows.push({ ...f, chainName: name });
      }
    }
    return rows;
  }, [results]);

  const filtered = React.useMemo(() => {
    let rows = allRows;
    if (riskFilter !== "All") rows = rows.filter(r => r.risk === riskFilter.toLowerCase());
    if (catFilter  !== "All") rows = rows.filter(r => r.category === catFilter);
    if (cdrOnly)               rows = rows.filter(r => r.cdr);

    const colKey = {
      pos:      r => r.pos1,
      residue:  r => r.residues[0],
      category: r => r.category,
      risk:     r => RISK_ORDER[r.risk],
      type:     r => r.display_label,
      cdr:      r => r.cdr || "",
      note:     r => r.note || "",
    }[sortCol] || (r => r.pos1);

    return [...rows].sort((a, b) => {
      const av = colKey(a), bv = colKey(b);
      const cmp = typeof av === "number"
        ? av - bv
        : String(av).localeCompare(String(bv));
      return sortAsc ? cmp : -cmp;
    });
  }, [allRows, riskFilter, catFilter, cdrOnly, sortCol, sortAsc]);

  const handleSort = col => {
    if (sortCol === col) setSortAsc(v => !v);
    else { setSortCol(col); setSortAsc(true); }
  };

  const thStyle = {
    padding: "8px 10px", textAlign: "left", fontWeight: 700,
    fontSize: "0.83em", color: "#444", background: THEME.bgPanel,
    borderBottom: `2px solid ${THEME.border}`, cursor: "pointer",
    userSelect: "none", whiteSpace: "nowrap",
  };
  const tdStyle = (risk) => ({
    padding: "7px 10px", fontSize: "0.85em", verticalAlign: "middle",
    borderBottom: `1px solid ${THEME.border}`,
    background: RISK_COLORS[risk]
      ? RISK_COLORS[risk].bg + "55"
      : "transparent",
  });

  if (!results || results.length === 0) return (
    <div style={{ display: "flex", alignItems: "center", justifyContent: "center",
                  height: "100%", color: THEME.fgMuted, fontSize: 16 }}>
      Run an analysis to see findings.
    </div>
  );

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
      {/* Filter bar */}
      <div style={{ background: THEME.bgPanel, padding: "8px 14px", display: "flex",
                    flexWrap: "wrap", gap: 12, alignItems: "center", borderBottom: `1px solid ${THEME.border}` }}>
        <span style={{ fontSize: "0.85em", color: THEME.fgMuted, fontWeight: 600 }}>Risk:</span>
        {["All", "High", "Medium", "Low"].map(r => (
          <label key={r} style={{ display: "flex", alignItems: "center", gap: 4,
                                   fontSize: "0.85em", cursor: "pointer" }}>
            <input type="radio" name="risk" value={r}
                   checked={riskFilter === r}
                   onChange={() => setRiskFilter(r)} />
            {r}
          </label>
        ))}
        <span style={{ fontSize: "0.85em", color: THEME.fgMuted, fontWeight: 600, marginLeft: 8 }}>Category:</span>
        <select value={catFilter} onChange={e => setCatFilter(e.target.value)}
                style={{ fontSize: "0.85em", padding: "2px 6px", borderRadius: 4,
                         border: `1px solid ${THEME.border}` }}>
          {["All", ...CAT_ORDER].map(c => <option key={c} value={c}>{c}</option>)}
        </select>
        <label style={{ display: "flex", alignItems: "center", gap: 4,
                         fontSize: "0.85em", cursor: "pointer", marginLeft: 8 }}>
          <input type="checkbox" checked={cdrOnly} onChange={e => setCdrOnly(e.target.checked)} />
          CDR only
        </label>
        <span style={{ marginLeft: "auto", fontSize: "0.82em", color: THEME.fgMuted }}>
          {filtered.length} finding{filtered.length !== 1 ? "s" : ""}
        </span>
      </div>

      {/* Table */}
      <div style={{ overflowY: "auto", flex: 1 }}>
        <table style={{ width: "100%", borderCollapse: "collapse", fontSize: "0.88em" }}>
          <thead>
            <tr>
              {[
                ["pos",      "Pos."],
                ["residue",  "Residue"],
                ["category", "Category"],
                ["risk",     "Risk"],
                ["type",     "Liability Type"],
                ["cdr",      "CDR"],
                ["note",     "Note"],
              ].map(([col, label]) => (
                <th key={col} style={thStyle} onClick={() => handleSort(col)}>
                  {label}<SortIcon col={col} sortCol={sortCol} sortAsc={sortAsc} />
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {filtered.map((f, i) => (
              <tr key={i}>
                <td style={{ ...tdStyle(f.risk), fontFamily: "monospace", color: THEME.fgMuted,
                              fontWeight: 700, textAlign: "center", width: 55 }}>
                  {f.pos1}
                </td>
                <td style={{ ...tdStyle(f.risk), textAlign: "center", width: 60 }}>
                  <span style={{
                    background: f.color_bg, color: f.color_text,
                    border: `1px solid ${f.color_border}`,
                    padding: "1px 6px", borderRadius: 3,
                    fontFamily: "monospace", fontWeight: 700,
                  }}>
                    {f.residues[0]}
                  </span>
                </td>
                <td style={{ ...tdStyle(f.risk), color: THEME.fgDark }}>{f.category}</td>
                <td style={{ ...tdStyle(f.risk), width: 80 }}>
                  <RiskBadge risk={f.risk} small />
                </td>
                <td style={{ ...tdStyle(f.risk), fontSize: "0.82em", color: "#555" }}>
                  {f.display_label}
                </td>
                <td style={{ ...tdStyle(f.risk), width: 80, textAlign: "center" }}>
                  {f.cdr ? (
                    <span style={{ background: "#e3f2fd", color: "#1565c0", padding: "1px 6px",
                                   borderRadius: 10, fontSize: "0.78em", fontWeight: 700 }}>
                      {f.cdr}
                    </span>
                  ) : "—"}
                </td>
                <td style={{ ...tdStyle(f.risk), fontSize: "0.8em", color: "#555", maxWidth: 260 }}>
                  {f.note || "—"}
                </td>
              </tr>
            ))}
            {filtered.length === 0 && (
              <tr>
                <td colSpan={7} style={{ textAlign: "center", padding: 24,
                                         color: THEME.fgMuted, fontStyle: "italic" }}>
                  No findings match the current filters.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
}


// ─── MAIN APP ─────────────────────────────────────────────────────────────────

function App() {
  const [inputText,       setInputText]       = React.useState("");
  const [results,         setResults]         = React.useState(null);
  const [selectedChain,   setSelectedChain]   = React.useState(0);
  const [activeTab,       setActiveTab]       = React.useState(0);
  const [analyzing,       setAnalyzing]       = React.useState(false);
  const [status,          setStatus]          = React.useState("Ready.");
  const [placeholder,     setPlaceholder]     = React.useState(true);
  const fileInputRef = React.useRef(null);

  const PLACEHOLDER = ">MyProtein_HC\nEVQLVESGGGLVQPGGSLRLSCAAGFNIKDTYIH...\n\nAccepts:\n  • FASTA (single or multi-chain)\n  • Raw amino-acid sequence\n  • Multiple FASTA entries in one paste";

  const handleFocus = () => { if (placeholder) { setInputText(""); setPlaceholder(false); } };
  const handleBlur  = () => { if (!inputText.trim()) { setInputText(PLACEHOLDER); setPlaceholder(true); } };

  const handleFileLoad = e => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = ev => {
      setInputText(ev.target.result.trim());
      setPlaceholder(false);
      setStatus(`Loaded: ${file.name}`);
    };
    reader.readAsText(file);
    e.target.value = "";
  };

  const handleAnalyze = () => {
    const raw = placeholder ? "" : inputText.trim();
    if (!raw) { alert("Please paste or load a protein sequence first."); return; }
    setAnalyzing(true);
    setStatus("Analyzing…");
    // Use setTimeout to allow UI to update before heavy computation
    setTimeout(() => {
      try {
        const pairs = raw.includes(">") ? parseFasta(raw) : [["Sequence", raw]];
        if (!pairs.length) { alert("No valid sequences found."); setAnalyzing(false); return; }
        const res = processSequences(pairs);
        setResults(res);
        setSelectedChain(0);
        setActiveTab(0);
        const total = res.reduce((s, r) => s + r.summary.total, 0);
        const highs = res.reduce((s, r) => s + (r.summary.byRisk.high || []).length, 0);
        setStatus(`Done — ${res.length} chain(s), ${total} liabilities flagged.`);
      } catch (err) {
        alert("Analysis error: " + err.message);
      }
      setAnalyzing(false);
    }, 20);
  };

  const handleClear = () => {
    setInputText(PLACEHOLDER);
    setPlaceholder(true);
    setResults(null);
    setStatus("Ready.");
  };

  // Download HTML report
  const handleDownloadReport = () => {
    if (!results) return;
    const html = buildHtmlReport(results);
    const blob = new Blob([html], { type: "text/html" });
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement("a");
    a.href     = url;
    a.download = "protein_liabilities.html";
    a.click();
    URL.revokeObjectURL(url);
    setStatus("HTML report downloaded.");
  };

  const currentResult = results ? results[selectedChain] : null;
  const totalCount    = results ? results.reduce((s, r) => s + r.summary.total, 0) : 0;
  const highCount     = results ? results.reduce((s, r) => s + (r.summary.byRisk.high || []).length, 0) : 0;

  const tabs = ["Summary", "Annotated Sequence", "All Findings"];

  return (
    <div style={{ height: "100vh", display: "flex", flexDirection: "column",
                  background: THEME.bgDark, fontFamily: "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif" }}>

      {/* ── Header ─────────────────────────────────────────────────────────── */}
      <div style={{ background: THEME.bgAccent, padding: "14px 22px",
                    display: "flex", justifyContent: "space-between", alignItems: "center" }}>
        <div>
          <div style={{ color: "#fff", fontWeight: 700, fontSize: 17 }}>
            🔬&nbsp; Protein Sequence Liability Analyzer
          </div>
          <div style={{ color: "#a8b2c8", fontSize: 10, marginTop: 3 }}>
            deamidation · oxidation · isomerization · glycosylation · pyroglutamate · truncation · CDR annotation
          </div>
        </div>
        <div style={{ color: "#7b8fb5", fontSize: 11, fontStyle: "italic" }}>
          Letarte Scientific Consulting
        </div>
      </div>

      {/* ── Body ───────────────────────────────────────────────────────────── */}
      <div style={{ flex: 1, display: "flex", overflow: "hidden",
                    background: THEME.bgLight, padding: "14px 14px 0" }}>

        {/* ── Left panel ─────────────────────────────────────────────────── */}
        <div style={{ width: 420, flexShrink: 0, display: "flex", flexDirection: "column",
                      marginRight: 12 }}>
          <div style={{ background: THEME.bgWhite, border: `1px solid ${THEME.border}`,
                        borderRadius: 8, flex: 1, display: "flex", flexDirection: "column",
                        overflow: "hidden" }}>

            {/* Card header */}
            <div style={{ background: THEME.bgPanel, padding: "10px 14px",
                          borderBottom: `1px solid ${THEME.border}`, fontWeight: 700, fontSize: 12 }}>
              Sequence Input
            </div>

            {/* Toolbar */}
            <div style={{ padding: "8px 10px", display: "flex", gap: 6,
                          borderBottom: `1px solid ${THEME.border}` }}>
              <input ref={fileInputRef} type="file"
                     accept=".fasta,.fa,.faa,.txt"
                     style={{ display: "none" }}
                     onChange={handleFileLoad} />
              <button
                onClick={() => fileInputRef.current?.click()}
                style={{ padding: "5px 12px", fontSize: 12, borderRadius: 5,
                         border: `1px solid ${THEME.border}`, background: THEME.bgPanel,
                         cursor: "pointer", fontWeight: 600 }}>
                📁 Load FASTA
              </button>
              <button
                onClick={handleClear}
                style={{ padding: "5px 12px", fontSize: 12, borderRadius: 5,
                         border: `1px solid ${THEME.border}`, background: THEME.bgPanel,
                         cursor: "pointer", fontWeight: 600 }}>
                🗑 Clear
              </button>
            </div>

            {/* Text area */}
            <textarea
              value={placeholder ? PLACEHOLDER : inputText}
              onChange={e => { setInputText(e.target.value); setPlaceholder(false); }}
              onFocus={handleFocus}
              onBlur={handleBlur}
              style={{
                flex: 1, resize: "none", border: "none", outline: "none",
                fontFamily: "'Courier New', monospace", fontSize: 11,
                padding: "10px 12px", background: THEME.bgPanel,
                color: placeholder ? THEME.fgMuted : THEME.fgDark,
                lineHeight: 1.5,
              }}
              spellCheck={false}
            />

            {/* Info strip */}
            <div style={{ background: "#f0f7ff", padding: "6px 12px", fontSize: 10,
                          color: "#2563eb", borderTop: `1px solid ${THEME.border}` }}>
              ℹ FASTA or raw sequence · single or multi-chain · CDR auto-detected
            </div>
          </div>

          {/* Analyze button */}
          <button
            onClick={handleAnalyze}
            disabled={analyzing}
            style={{
              marginTop: 8, padding: "14px", fontSize: 14, fontWeight: 700,
              background: analyzing ? "#6b7280" : THEME.btnBlue,
              color: "#ffffff", border: "none", borderRadius: 6, cursor: analyzing ? "wait" : "pointer",
              letterSpacing: "0.5px", transition: "background 0.15s",
            }}
            onMouseEnter={e => { if (!analyzing) e.target.style.background = "#1a4080"; }}
            onMouseLeave={e => { if (!analyzing) e.target.style.background = THEME.btnBlue; }}
          >
            {analyzing ? "  Analyzing…  " : "  ANALYZE  "}
          </button>
        </div>

        {/* ── Right panel ────────────────────────────────────────────────── */}
        <div style={{ flex: 1, display: "flex", flexDirection: "column",
                      background: THEME.bgWhite, border: `1px solid ${THEME.border}`,
                      borderRadius: 8, overflow: "hidden" }}>

          {/* Tab bar + chain selector */}
          <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between",
                        borderBottom: `2px solid ${THEME.border}`, background: THEME.bgPanel,
                        padding: "0 14px", flexShrink: 0 }}>
            <div style={{ display: "flex" }}>
              {tabs.map((tab, i) => (
                <button key={i} onClick={() => setActiveTab(i)}
                  style={{
                    padding: "12px 16px", fontSize: 13, fontWeight: activeTab === i ? 700 : 400,
                    color: activeTab === i ? THEME.fgDark : THEME.fgMuted,
                    background: "transparent", border: "none", cursor: "pointer",
                    borderBottom: activeTab === i ? `2px solid ${THEME.btnBlue}` : "2px solid transparent",
                    marginBottom: -2,
                  }}>
                  {tab}
                </button>
              ))}
            </div>
            {/* Chain selector (only when multiple chains) */}
            {results && results.length > 1 && (
              <select value={selectedChain}
                      onChange={e => setSelectedChain(Number(e.target.value))}
                      style={{ fontSize: 12, padding: "3px 8px", borderRadius: 4,
                               border: `1px solid ${THEME.border}` }}>
                {results.map((r, i) => (
                  <option key={i} value={i}>{r.name}</option>
                ))}
              </select>
            )}
          </div>

          {/* Tab content */}
          <div style={{ flex: 1, overflow: "hidden", display: "flex", flexDirection: "column" }}>
            {activeTab === 0 && <SummaryTab result={currentResult} />}
            {activeTab === 1 && (
              <AnnotatedSequenceTab
                results={results ? (results.length > 1 ? results : results) : null}
              />
            )}
            {activeTab === 2 && <AllFindingsTab results={results} />}
          </div>
        </div>
      </div>

      {/* ── Footer ─────────────────────────────────────────────────────────── */}
      <div style={{
        background: THEME.bgMid, padding: "8px 14px",
        display: "flex", alignItems: "center", justifyContent: "space-between",
      }}>
        <span style={{ color: "#7b8fb5", fontSize: 11 }}>{status}</span>
        <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
          {results && (
            <span style={{ color: "#7b8fb5", fontSize: 11, fontWeight: 700, marginRight: 8 }}>
              {totalCount} liabilit{totalCount !== 1 ? "ies" : "y"} · {highCount} high-risk
            </span>
          )}
          <button onClick={handleDownloadReport} disabled={!results}
            style={{
              padding: "5px 12px", fontSize: 11, borderRadius: 4,
              background: results ? "#1e293b" : "#374151",
              color: results ? "#e2e8f0" : "#6b7280",
              border: "1px solid #374151", cursor: results ? "pointer" : "default",
              fontWeight: 600,
            }}>
            💾 Download HTML Report
          </button>
        </div>
      </div>
    </div>
  );
}


// ─── HTML REPORT GENERATOR ────────────────────────────────────────────────────

function buildHtmlReport(results) {
  const now = new Date().toISOString().slice(0, 16).replace("T", " ");
  const catIcons = { Deamidation:"🧪", Oxidation:"⚡", Isomerization:"🔄",
                     Glycosylation:"🍬", Pyroglutamate:"🔵", Truncation:"✂️" };

  const sections = results.map(({ name, seq, findings, summary, cdrs }) => {
    const highCt = (summary.byRisk.high || []).length;

    const cards = CAT_ORDER.map(cat => {
      const items = summary.byCategory[cat] || [];
      const h = items.filter(f => f.risk === "high").length;
      const m = items.filter(f => f.risk === "medium").length;
      const l = items.filter(f => f.risk === "low").length;
      const bars = [
        h ? `<span class="rb risk-high">${h} High</span>` : "",
        m ? `<span class="rb risk-medium">${m} Med</span>` : "",
        l ? `<span class="rb risk-low">${l} Low</span>` : "",
      ].filter(Boolean).join(" ");
      return `<div class="scard"><div class="si">${catIcons[cat]||"•"}</div>
        <div class="sc">${items.length}</div><div class="st">${cat}</div>
        <div class="sb">${bars||"—"}</div></div>`;
    }).join("");

    const cdrInfo = cdrs.length > 0
      ? `<div class="cdr-info">📐 <b>CDRs detected:</b> ` +
        cdrs.map(([s,e,n]) => `<b>${n}</b> (${s+1}–${e+1})`).join(", ") +
        ` <em>(Chothia-like approx.)</em></div>`
      : "";

    // Annotated sequence
    const n = seq.length;
    const posF = Array.from({ length: n }, () => []);
    for (const f of findings) {
      for (let j = f.pos0; j < Math.min(f.pos0 + f.span, n); j++) posF[j].push(f);
    }
    const cdrPos = {};
    for (const [cs, ce, cn] of cdrs) for (let p = cs; p <= ce; p++) cdrPos[p] = cn;

    const WRAP = 60;
    const seqLines = [];
    for (let start = 0; start < n; start += WRAP) {
      const end = Math.min(start + WRAP, n);
      let line = `<span class="pl">${String(start + 1).padStart(5)} </span>`;
      for (let i = start; i < end; i++) {
        const aa = seq[i];
        const fl = posF[i];
        const inCdr = cdrPos[i] !== undefined;
        if (fl.length > 0) {
          const best = fl.reduce((a, b) => RISK_ORDER[a.risk] <= RISK_ORDER[b.risk] ? a : b);
          const tt = `Pos ${i+1}: ${best.short}${inCdr?" ["+cdrPos[i]+"]":""}`;
          const ud = inCdr ? "text-decoration:underline;text-decoration-thickness:2px;" : "";
          line += `<span style="background:${best.color_bg};color:${best.color_text};border-bottom:2px solid ${best.color_border};${ud}font-weight:700" title="${tt}">${aa}</span>`;
        } else if (inCdr) {
          line += `<span class="cdr-aa" title="Pos ${i+1}: ${cdrPos[i]}">${aa}</span>`;
        } else {
          line += aa;
        }
      }
      const cdrTags = cdrs
        .filter(([cs,ce]) => cs < end && ce >= start)
        .map(([cs,ce,cn]) => `<span class="cdr-lbl">${cn}:${Math.max(cs,start)+1}–${Math.min(ce,end-1)+1}</span>`)
        .join("");
      line += ` <span class="pl">${end}</span>` + (cdrTags ? "  " + cdrTags : "");
      seqLines.push(line);
    }

    const sortedF = [...findings].sort((a, b) =>
      (RISK_ORDER[a.risk] - RISK_ORDER[b.risk]) || (a.pos0 - b.pos0));
    const rows = sortedF.map(f => {
      const rb = `<span class="rb risk-${f.risk}">${f.risk.toUpperCase()}</span>`;
      const res = `<span style="background:${f.color_bg};color:${f.color_text};border:1px solid ${f.color_border};padding:1px 5px;border-radius:3px;font-family:monospace;font-weight:700">${f.residues[0]}</span>`;
      const cdrBadge = f.cdr ? `<span class="cdr-badge">${f.cdr}</span>` : "";
      const note = (f.note || "—") + (cdrBadge ? " " + cdrBadge : "");
      return `<tr><td class="pc">${f.pos1}</td><td>${res}</td><td>${f.category}</td>` +
             `<td>${rb}</td><td>${f.display_label}</td><td>${note}</td></tr>`;
    }).join("");

    return `<div class="sec">
      <div class="sh">${name}
        <span class="sm">${summary.seqLen} residues · ${summary.total} findings${cdrs.length ? " · " + cdrs.length + " CDR(s)" : ""}</span>
      </div>
      <div class="sb2">
        <div class="alert" style="background:${highCt?"#fde8e8":"#d6eaf8"};color:${highCt?"#c0392b":"#154360"}">
          ${highCt ? "⚠️ " + highCt + " high-risk liabilit" + (highCt===1?"y":"ies") + " detected." : "✅ No high-risk liabilities detected."}
        </div>
        ${cdrInfo}
        <div class="sgrid">${cards}</div>
        <div class="sl">Annotated Sequence <span class="sl-note">hover for details · underlined = CDR</span></div>
        <div class="seq-blk">${seqLines.join("\n")}</div>
        <div class="sl" style="margin-top:18px">Liabilities Detail</div>
        <div style="overflow-x:auto">
          <table><thead><tr>
            <th>Pos.</th><th>Residue</th><th>Category</th>
            <th>Risk</th><th>Liability Type</th><th>Note / CDR</th>
          </tr></thead><tbody>${rows}</tbody></table>
        </div>
      </div>
    </div>`;
  }).join("");

  const legRows = Object.entries(LIABILITIES).map(([, lib]) =>
    `<tr><td><span style="background:${lib.color_bg};color:${lib.color_text};border:1px solid ${lib.color_border};padding:2px 6px;border-radius:3px;font-family:monospace;font-weight:600">${lib.short}</span></td>` +
    `<td>${lib.label}</td><td><span class="rb risk-${lib.risk}">${lib.risk.toUpperCase()}</span></td>` +
    `<td style="color:#555;font-size:.85em">${lib.description}</td></tr>`
  ).join("");

  return `<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Protein Sequence Liability Analysis</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;background:#f0f2f5;color:#1a1a2e;font-size:15px}
.pw{max-width:1200px;margin:0 auto;padding:24px 16px 60px}
.hdr{background:linear-gradient(135deg,#1a1a2e,#0f3460);color:#fff;border-radius:12px;padding:32px 36px;margin-bottom:28px}
.hdr h1{font-size:1.9em;font-weight:700;margin-bottom:6px}
.hdr .sub{color:#a8b2c8;font-size:.95em}
.hdr .meta{margin-top:16px;font-size:.85em;color:#c0c8dc;display:flex;gap:24px;flex-wrap:wrap}
.sec{background:#fff;border-radius:10px;box-shadow:0 1px 6px rgba(0,0,0,.08);margin-bottom:24px;overflow:hidden}
.sh{padding:16px 24px;border-bottom:1px solid #eee;font-weight:700;font-size:1.05em;background:#fafbfc}
.sm{font-weight:400;font-size:.88em;color:#888;margin-left:10px}
.sb2{padding:20px 24px}
.alert{padding:10px 16px;border-radius:6px;margin-bottom:14px;font-weight:600;font-size:.9em}
.cdr-info{background:#e3f2fd;color:#1565c0;padding:8px 14px;border-radius:6px;margin-bottom:12px;font-size:.88em}
.sgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:12px;margin-bottom:20px}
.scard{background:#fff;border:1px solid #e5e7eb;border-radius:8px;padding:12px 10px;text-align:center}
.si{font-size:1.6em;margin-bottom:4px}.sc{font-size:1.8em;font-weight:700;color:#1a1a2e}
.st{font-size:.72em;font-weight:700;color:#6b7280;margin:2px 0 6px}
.sb{display:flex;flex-wrap:wrap;gap:3px;justify-content:center}
.rb{display:inline-block;padding:2px 7px;border-radius:20px;font-size:.72em;font-weight:700;letter-spacing:.4px;white-space:nowrap}
.risk-high{background:#fde8e8;color:#c0392b;border:1px solid #e74c3c}
.risk-medium{background:#fef0e6;color:#d35400;border:1px solid #e67e22}
.risk-low{background:#e8f8f5;color:#0e6655;border:1px solid #1abc9c}
.risk-info{background:#d6eaf8;color:#154360;border:1px solid #2e86c1}
.sl{font-weight:700;margin-bottom:8px;font-size:.95em}.sl-note{font-weight:400;font-size:.78em;color:#888;margin-left:8px}
.seq-blk{font-family:"Courier New",monospace;font-size:.9em;line-height:2.1;background:#f8f9fa;border-radius:8px;padding:14px;overflow-x:auto;white-space:pre-wrap;word-break:break-all}
.pl{color:#9ba5b0;user-select:none;font-size:.86em}
.cdr-aa{text-decoration:underline;text-decoration-thickness:2px;text-underline-offset:3px;color:#1565c0}
.cdr-lbl{color:#1a5276;font-size:.76em;font-weight:600;margin-left:6px}
.cdr-badge{display:inline-block;padding:1px 6px;border-radius:10px;background:#e3f2fd;color:#1565c0;font-size:.75em;font-weight:700;margin-left:4px}
table{width:100%;border-collapse:collapse;font-size:.88em}
thead th{background:#f4f5f7;padding:10px 12px;text-align:left;font-weight:700;color:#444;border-bottom:2px solid #e0e3e8}
tbody tr{border-bottom:1px solid #f0f0f0}tbody tr:hover{background:#fafbff}
tbody td{padding:8px 12px;vertical-align:middle}.pc{font-family:monospace;color:#7f8c8d;font-weight:600}
.footer{text-align:center;margin-top:40px;font-size:.82em;color:#aaa}
</style></head><body>
<div class="pw">
  <div class="hdr">
    <h1>🔬 Protein Sequence Liability Analysis</h1>
    <div class="sub">PTM risk · deamidation · oxidation · isomerization · glycosylation · pyroglutamate · truncation</div>
    <div class="meta"><span>📅 ${now}</span><span>🔢 ${results.length} chain(s)</span><span>⚗️ Letarte Scientific Consulting</span></div>
  </div>
  ${sections}
  <div class="sec">
    <div class="sh">📖 Liability Legend</div>
    <div class="sb2">
      <table><thead><tr><th>Code</th><th>Liability</th><th>Risk</th><th>Description</th></tr></thead>
      <tbody>${legRows}</tbody></table>
      <p style="margin-top:14px;font-size:.82em;color:#888">CDR boundaries are approximate (Chothia-like). Tryptic K/R sites excluded from truncation.</p>
    </div>
  </div>
  <div class="footer">Protein Liability Analyzer · Letarte Scientific Consulting · For R&amp;D use only</div>
</div></body></html>`;
}


export default App;
