import { useState, useRef } from "react";

// ─── Monoisotopic residue masses (Da) ───────────────────────────────────────
const AA_MASSES = {
  G: 57.02146, A: 71.03711, V: 99.06841, L: 113.08406, I: 113.08406,
  P: 97.05276,  F: 147.06841, W: 186.07931, M: 131.04049, S: 87.03203,
  T: 101.04768, C: 103.00919, Y: 163.06333, H: 137.05891, D: 115.02694,
  E: 129.04259, N: 114.04293, Q: 128.05858, K: 128.09496, R: 156.10111,
};
const WATER = 18.01056;

// ─── Cysteine modifications ──────────────────────────────────────────────────
const CYS_MODS = {
  unmodified:      { label: "Unmodified",                  delta: 0,         info: "Free thiol (SH)" },
  carbamidomethyl: { label: "Carbamidomethyl (CAM)",       delta: 57.02146,  info: "+C₂H₃NO  (iodoacetamide)" },
  carboxymethyl:   { label: "Carboxymethyl (CM)",          delta: 58.00548,  info: "+C₂H₂O₂  (iodoacetic acid)" },
};

// ─── Enzyme definitions ──────────────────────────────────────────────────────
const ENZYMES = {
  trypsin: {
    name: "Trypsin",
    desc: "Cleaves C-terminal to K, R  (not before P)",
    sites: (seq) => {
      const s = [];
      for (let i = 0; i < seq.length - 1; i++)
        if ((seq[i] === "K" || seq[i] === "R") && seq[i + 1] !== "P") s.push(i);
      return s;
    },
  },
  chymotrypsin: {
    name: "Chymotrypsin",
    desc: "Cleaves C-terminal to F, Y, W  (not before P)",
    sites: (seq) => {
      const s = [];
      for (let i = 0; i < seq.length - 1; i++)
        if (["F", "Y", "W"].includes(seq[i]) && seq[i + 1] !== "P") s.push(i);
      return s;
    },
  },
  lysc: {
    name: "Lys-C",
    desc: "Cleaves C-terminal to K",
    sites: (seq) => {
      const s = [];
      for (let i = 0; i < seq.length - 1; i++)
        if (seq[i] === "K") s.push(i);
      return s;
    },
  },
  aspn: {
    name: "Asp-N",
    desc: "Cleaves N-terminal to D  (before D)",
    sites: (seq) => {
      const s = [];
      for (let i = 1; i < seq.length; i++)
        if (seq[i] === "D") s.push(i - 1);
      return s;
    },
  },
  thermolysin: {
    name: "Thermolysin",
    desc: "Cleaves N-terminal to L, I, V, F, A, M",
    sites: (seq) => {
      const s = [];
      for (let i = 1; i < seq.length; i++)
        if (["L", "I", "V", "F", "A", "M"].includes(seq[i])) s.push(i - 1);
      return s;
    },
  },
  ides: {
    name: "IdeS",
    desc: "Cleaves IgG1 hinge  ELLGG↓PSVF",
    sites: (seq) => {
      const s = [];
      let idx = seq.indexOf("ELLGGP");
      while (idx !== -1) {
        s.push(idx + 4); // after 2nd G, before P
        idx = seq.indexOf("ELLGGP", idx + 1);
      }
      return s;
    },
  },
};

// ─── Mass & digestion logic ──────────────────────────────────────────────────
function peptideMass(pep, cysMod) {
  const d = CYS_MODS[cysMod].delta;
  let m = WATER;
  for (const aa of pep) {
    const rm = AA_MASSES[aa];
    if (rm !== undefined) { m += rm; if (aa === "C") m += d; }
  }
  return m;
}

function digest(seq, enzymeKey, minLen, cysMod) {
  const cuts = ENZYMES[enzymeKey].sites(seq);
  const bounds = [0, ...cuts.map((c) => c + 1), seq.length];
  const out = [];
  for (let i = 0; i < bounds.length - 1; i++) {
    const pep = seq.slice(bounds[i], bounds[i + 1]);
    if (pep.length < minLen) continue;
    out.push({
      idx:    i + 1,
      seq:    pep,
      start:  bounds[i] + 1,
      end:    bounds[i + 1],
      len:    pep.length,
      nCys:   (pep.match(/C/g) || []).length,
      mass:   peptideMass(pep, cysMod),
    });
  }
  return out;
}

function parseInput(raw) {
  if (!raw.trim()) return { clean: "", invalid: [] };
  let seq = raw.trim().startsWith(">")
    ? raw.split("\n").filter((l) => !l.startsWith(">")).join("")
    : raw;
  seq = seq.replace(/[^A-Za-z]/g, "").toUpperCase();
  const invalid = [...new Set(seq.split("").filter((aa) => !AA_MASSES[aa]))];
  return { clean: seq, invalid };
}

// ─── Icons / helpers ─────────────────────────────────────────────────────────
const SortArrow = ({ col, sortKey, dir }) =>
  sortKey !== col ? (
    <span style={{ color: "#94a3b8", marginLeft: 4 }}>↕</span>
  ) : (
    <span style={{ color: "#60a5fa", marginLeft: 4 }}>{dir === "asc" ? "↑" : "↓"}</span>
  );

// ─── Main component ──────────────────────────────────────────────────────────
export default function ProteinDigestionTool() {
  const [rawSeq, setRawSeq] = useState("");
  const [enzyme, setEnzyme] = useState("trypsin");
  const [minLen, setMinLen] = useState(6);
  const [cysMod, setCysMod] = useState("carbamidomethyl");
  const [results, setResults] = useState(null);
  const [error, setError] = useState("");
  const [sortKey, setSortKey] = useState("start");
  const [sortDir, setSortDir] = useState("asc");
  const [filter, setFilter] = useState("");
  const fileRef = useRef(null);

  // ── handlers ──
  const loadFile = (e) => {
    const f = e.target.files[0];
    if (!f) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
      setRawSeq(ev.target.result);
      setError("");
    };
    reader.readAsText(f);
  };

  const runDigest = () => {
    const { clean, invalid } = parseInput(rawSeq);
    if (!clean) { setError("Please enter or load a protein sequence."); return; }
    if (invalid.length) {
      setError(`Unknown characters found: ${invalid.join(", ")}. Please check your sequence.`);
      return;
    }
    setError("");
    setFilter("");
    const peptides = digest(clean, enzyme, minLen, cysMod);
    setResults({ peptides, seqLen: clean.length, enzyme, cysMod, minLen });
  };

  const handleSort = (k) => {
    if (k === sortKey) setSortDir((d) => (d === "asc" ? "desc" : "asc"));
    else { setSortKey(k); setSortDir("asc"); }
  };

  const displayRows = () => {
    if (!results) return [];
    let rows = results.peptides;
    if (filter) {
      const q = filter.toUpperCase();
      rows = rows.filter((r) => r.seq.includes(q));
    }
    return [...rows].sort((a, b) => {
      const cmp = typeof a[sortKey] === "string"
        ? a[sortKey].localeCompare(b[sortKey])
        : a[sortKey] - b[sortKey];
      return sortDir === "asc" ? cmp : -cmp;
    });
  };

  const exportCSV = () => {
    const rows = displayRows();
    const header = ["#", "Sequence", "Start", "End", "Length", "Cys", "Monoisotopic Mass (Da)"];
    const lines = [header, ...rows.map((r) =>
      [r.idx, r.seq, r.start, r.end, r.len, r.nCys, r.mass.toFixed(5)]
    )].map((r) => r.join(",")).join("\n");
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines], { type: "text/csv" }));
    a.download = `digestion_${enzyme}_${new Date().toISOString().slice(0,10)}.csv`;
    a.click();
  };

  const { clean: seqClean } = parseInput(rawSeq);
  const rows = displayRows();

  // ── styles ──
  const card = {
    background: "white", borderRadius: 12, padding: 20,
    boxShadow: "0 1px 4px rgba(0,0,0,0.09)", marginBottom: 16,
  };
  const btn = (bg = "#1e3a5f") => ({
    background: bg, color: "white", border: "none", borderRadius: 8,
    padding: "10px 18px", cursor: "pointer", fontSize: 13, fontWeight: 600,
  });
  const label14 = { fontSize: 14, fontWeight: 600, color: "#1e3a5f", marginBottom: 10, display: "block" };
  const radioRow = { display: "flex", alignItems: "flex-start", marginBottom: 10, cursor: "pointer", gap: 8 };

  return (
    <div style={{ fontFamily: "'Inter', system-ui, sans-serif", maxWidth: 1100, margin: "0 auto", padding: 24, background: "#f1f5f9", minHeight: "100vh" }}>

      {/* ── Header ── */}
      <div style={{ background: "linear-gradient(135deg,#1e3a5f 0%,#1a5fa8 100%)", borderRadius: 14, padding: "22px 28px", marginBottom: 20, color: "white" }}>
        <div style={{ display: "flex", alignItems: "center", gap: 12 }}>
          <span style={{ fontSize: 36 }}>🧬</span>
          <div>
            <h1 style={{ margin: 0, fontSize: 22, fontWeight: 700, letterSpacing: "-0.3px" }}>Protein Digestion Tool</h1>
            <p style={{ margin: "5px 0 0", opacity: 0.8, fontSize: 13 }}>
              In silico enzymatic digestion · Exact monoisotopic masses · Cysteine alkylation support
            </p>
          </div>
        </div>
      </div>

      {/* ── Sequence input ── */}
      <div style={card}>
        <span style={label14}>① Protein Sequence</span>
        <div style={{ display: "flex", gap: 10, marginBottom: 10, flexWrap: "wrap" }}>
          <input ref={fileRef} type="file" accept=".fasta,.fa,.txt" onChange={loadFile} style={{ display: "none" }} />
          <button style={btn()} onClick={() => fileRef.current?.click()}>📂 Load FASTA / .txt file</button>
          <button style={btn("#64748b")} onClick={() => { setRawSeq(""); setResults(null); setError(""); }}>✕ Clear</button>
        </div>
        <textarea
          value={rawSeq}
          onChange={(e) => setRawSeq(e.target.value)}
          placeholder={"Paste your sequence here (plain AA or FASTA format):\n\n>sp|P01857|IGHG1_HUMAN Immunoglobulin heavy constant gamma 1\nAST KGPSVFPLAP SSKSTSGGTA ALGCLVKDYF PEPVTVSWNSG ALTSGVHTFP..."}
          style={{
            width: "100%", height: 130, padding: 12, border: "1.5px solid #e2e8f0",
            borderRadius: 8, fontSize: 13, fontFamily: "monospace", resize: "vertical",
            boxSizing: "border-box", outline: "none", lineHeight: 1.6,
          }}
        />
        {seqClean.length > 0 && (
          <p style={{ margin: "6px 0 0", fontSize: 12, color: "#64748b" }}>
            ✅ Sequence loaded: <strong>{seqClean.length.toLocaleString()} residues</strong>
          </p>
        )}
      </div>

      {/* ── Settings row ── */}
      <div style={{ display: "grid", gridTemplateColumns: "1.2fr 1fr 0.9fr", gap: 16, marginBottom: 16 }}>

        {/* Enzyme */}
        <div style={card}>
          <span style={label14}>② Digestion Enzyme</span>
          {Object.entries(ENZYMES).map(([key, val]) => (
            <label key={key} style={radioRow}>
              <input type="radio" name="enzyme" value={key}
                checked={enzyme === key} onChange={() => setEnzyme(key)}
                style={{ marginTop: 2, accentColor: "#1e3a5f" }} />
              <span>
                <span style={{ fontWeight: 600, fontSize: 13, color: "#1e293b" }}>{val.name}</span>
                <span style={{ display: "block", fontSize: 11, color: "#64748b", marginTop: 1 }}>{val.desc}</span>
              </span>
            </label>
          ))}
        </div>

        {/* Cys mod */}
        <div style={card}>
          <span style={label14}>③ Cysteine Modification</span>
          {Object.entries(CYS_MODS).map(([key, val]) => (
            <label key={key} style={radioRow}>
              <input type="radio" name="cys" value={key}
                checked={cysMod === key} onChange={() => setCysMod(key)}
                style={{ marginTop: 2, accentColor: "#1e3a5f" }} />
              <span>
                <span style={{ fontWeight: 600, fontSize: 13, color: "#1e293b" }}>{val.label}</span>
                <span style={{ display: "block", fontSize: 11, color: "#64748b", marginTop: 1 }}>
                  {val.delta > 0 ? `Δ +${val.delta.toFixed(3)} Da · ` : ""}{val.info}
                </span>
              </span>
            </label>
          ))}
        </div>

        {/* Parameters + Run */}
        <div style={card}>
          <span style={label14}>④ Parameters</span>
          <div style={{ marginBottom: 6, fontSize: 13, fontWeight: 500, color: "#374151" }}>
            Minimum peptide length
          </div>
          <div style={{ display: "flex", alignItems: "center", gap: 12, marginBottom: 6 }}>
            <input
              type="range" min={1} max={25} value={minLen}
              onChange={(e) => setMinLen(+e.target.value)}
              style={{ flex: 1, accentColor: "#1e3a5f" }}
            />
            <span style={{
              background: "#1e3a5f", color: "white", borderRadius: 6,
              padding: "3px 10px", fontSize: 15, fontWeight: 700, minWidth: 34, textAlign: "center"
            }}>
              {minLen}
            </span>
          </div>
          <p style={{ fontSize: 11, color: "#64748b", margin: "0 0 20px" }}>
            Peptides &lt; {minLen} AA will be excluded
          </p>

          <button
            onClick={runDigest}
            style={{
              ...btn("linear-gradient(135deg,#1e3a5f,#1a5fa8)"),
              width: "100%", padding: "13px", fontSize: 14, borderRadius: 10,
              boxShadow: "0 2px 8px rgba(30,58,95,0.25)"
            }}
          >
            ⚗️ Run Digestion
          </button>
        </div>
      </div>

      {/* ── Error ── */}
      {error && (
        <div style={{ background: "#fef2f2", border: "1.5px solid #fca5a5", borderRadius: 8, padding: "11px 16px", marginBottom: 16, color: "#b91c1c", fontSize: 13 }}>
          ⚠️ {error}
        </div>
      )}

      {/* ── Results ── */}
      {results && (
        <div style={card}>
          {/* Results header */}
          <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 14, flexWrap: "wrap", gap: 10 }}>
            <div>
              <h2 style={{ margin: 0, fontSize: 16, fontWeight: 700, color: "#1e3a5f" }}>Digestion Results</h2>
              <p style={{ margin: "3px 0 0", fontSize: 12, color: "#64748b" }}>
                {ENZYMES[results.enzyme].name} &nbsp;·&nbsp;
                {CYS_MODS[results.cysMod].label} &nbsp;·&nbsp;
                Min length {results.minLen} AA &nbsp;·&nbsp;
                <strong style={{ color: "#1e3a5f" }}>{results.peptides.length}</strong> peptides generated
                {rows.length !== results.peptides.length &&
                  <span style={{ color: "#e05c00" }}> ({rows.length} shown after filter)</span>}
              </p>
            </div>
            <div style={{ display: "flex", gap: 9, flexWrap: "wrap" }}>
              <input
                type="text"
                placeholder="🔍 Filter by sequence…"
                value={filter}
                onChange={(e) => setFilter(e.target.value)}
                style={{ padding: "8px 12px", border: "1.5px solid #e2e8f0", borderRadius: 8, fontSize: 13, width: 200 }}
              />
              <button style={btn("#16a34a")} onClick={exportCSV}>📥 Export CSV</button>
            </div>
          </div>

          {/* Summary stats */}
          <div style={{ display: "grid", gridTemplateColumns: "repeat(4,1fr)", gap: 10, marginBottom: 16 }}>
            {[
              { label: "Protein Length",  value: `${results.seqLen.toLocaleString()} AA` },
              { label: "Peptides",        value: results.peptides.length },
              { label: "Cys-containing",  value: results.peptides.filter((p) => p.nCys > 0).length },
              {
                label: "Mass range",
                value: results.peptides.length
                  ? `${Math.min(...results.peptides.map((p) => p.mass)).toFixed(0)} – ${Math.max(...results.peptides.map((p) => p.mass)).toFixed(0)} Da`
                  : "—"
              },
            ].map((s) => (
              <div key={s.label} style={{ background: "#f0f6ff", borderRadius: 9, padding: "12px 14px", textAlign: "center" }}>
                <div style={{ fontSize: 17, fontWeight: 800, color: "#1e3a5f" }}>{s.value}</div>
                <div style={{ fontSize: 11, color: "#64748b", marginTop: 2 }}>{s.label}</div>
              </div>
            ))}
          </div>

          {/* Table */}
          <div style={{ overflowX: "auto" }}>
            <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 13 }}>
              <thead>
                <tr style={{ background: "#1e3a5f", color: "white" }}>
                  {[
                    { k: "idx",   label: "#",                      align: "center" },
                    { k: "seq",   label: "Peptide Sequence",        align: "left" },
                    { k: "start", label: "Start",                   align: "center" },
                    { k: "end",   label: "End",                     align: "center" },
                    { k: "len",   label: "Length",                  align: "center" },
                    { k: "nCys",  label: "Cys",                     align: "center" },
                    { k: "mass",  label: "Monoisotopic Mass (Da)",   align: "right" },
                  ].map((col) => (
                    <th
                      key={col.k}
                      onClick={() => handleSort(col.k)}
                      style={{
                        padding: "10px 12px", textAlign: col.align, cursor: "pointer",
                        userSelect: "none", whiteSpace: "nowrap", fontWeight: 600,
                      }}
                    >
                      {col.label}
                      <SortArrow col={col.k} sortKey={sortKey} dir={sortDir} />
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {rows.length === 0 ? (
                  <tr>
                    <td colSpan={7} style={{ padding: 30, textAlign: "center", color: "#94a3b8" }}>
                      No peptides match the current filter
                    </td>
                  </tr>
                ) : (
                  rows.map((r, i) => (
                    <tr
                      key={r.idx}
                      style={{
                        background: i % 2 === 0 ? "white" : "#f8fafc",
                        borderBottom: "1px solid #e2e8f0",
                      }}
                    >
                      <td style={{ padding: "8px 12px", textAlign: "center", color: "#94a3b8", fontSize: 11 }}>{r.idx}</td>
                      <td style={{ padding: "8px 12px", fontFamily: "monospace", color: "#1e293b", letterSpacing: "0.4px" }}>{r.seq}</td>
                      <td style={{ padding: "8px 12px", textAlign: "center", color: "#475569" }}>{r.start}</td>
                      <td style={{ padding: "8px 12px", textAlign: "center", color: "#475569" }}>{r.end}</td>
                      <td style={{ padding: "8px 12px", textAlign: "center", color: "#475569" }}>{r.len}</td>
                      <td style={{
                        padding: "8px 12px", textAlign: "center",
                        color: r.nCys > 0 ? "#b91c1c" : "#475569",
                        fontWeight: r.nCys > 0 ? 700 : 400,
                      }}>
                        {r.nCys > 0 ? `⚠ ${r.nCys}` : "0"}
                      </td>
                      <td style={{ padding: "8px 12px", textAlign: "right", fontFamily: "monospace", fontWeight: 600, color: "#1e3a5f" }}>
                        {r.mass.toFixed(5)}
                      </td>
                    </tr>
                  ))
                )}
              </tbody>
            </table>
          </div>

          <p style={{ margin: "10px 0 0", fontSize: 11, color: "#94a3b8" }}>
            † Masses are monoisotopic and include water (+18.01056 Da). Cysteine modification applied: {CYS_MODS[results.cysMod].label}
            {CYS_MODS[results.cysMod].delta > 0 && ` (Δ +${CYS_MODS[results.cysMod].delta.toFixed(5)} Da per Cys)`}.
          </p>
        </div>
      )}
    </div>
  );
}
