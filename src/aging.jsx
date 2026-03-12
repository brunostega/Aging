import { useState, useEffect, useRef, useCallback, useMemo } from "react";

const DEFAULT_N = 80;
const STATES = { M: 0, D: 1, F: 2 };
const STATE_NAMES = ["Monomer", "Disordered", "Fibril"];
const STATE_COLORS = ["#4ade80", "#fb923c", "#f87171"];

// Quenched disorder matrix: jMatrix[i][j] in [0,1], symmetric
function buildJMatrix(n) {
  const m = Array.from({ length: n }, () => new Float32Array(n));
  for (let i = 0; i < n; i++)
    for (let j = i + 1; j < n; j++) {
      const v = Math.random();
      m[i][j] = v;
      m[j][i] = v;
    }
  return m;
}

function initChain(n) { return new Array(n).fill(STATES.M); }

// Length of the contiguous fibril run containing site i
function fibrilRunLength(chain, i) {
  if (chain[i] !== STATES.F) return 0;
  const N = chain.length;
  let lo = i, hi = i;
  while (lo > 0 && chain[lo - 1] === STATES.F) lo--;
  while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;
  return hi - lo + 1;
}

function computeEnergy(chain, params, jMatrix) {
  const { eM, eD, eF, jF, jFF, jD, rD, minRun } = params;
  const N = chain.length;
  const baseE = [eM, eD, eF];
  let E = 0;
  for (let i = 0; i < N; i++) {
    E += baseE[chain[i]];
    if (chain[i] === STATES.F && fibrilRunLength(chain, i) >= minRun) {
      // Short-range: direct neighbors
      if (i > 0     && chain[i-1] === STATES.F && fibrilRunLength(chain, i-1) >= minRun) E -= jF / 2;
      if (i < N - 1 && chain[i+1] === STATES.F && fibrilRunLength(chain, i+1) >= minRun) E -= jF / 2;
      // Long-range: any other active fibril site (|d| > 1)
      for (let j = 0; j < N; j++) {
        if (Math.abs(j - i) > 1 && chain[j] === STATES.F && fibrilRunLength(chain, j) >= minRun) E -= jFF / 2;
      }
    }
    if (chain[i] === STATES.D) {
      for (let d = 1; d <= rD; d++) {
        if (i - d >= 0 && chain[i-d] === STATES.D) E -= (jD * jMatrix[i][i-d]) / 2;
        if (i + d <  N && chain[i+d] === STATES.D) E -= (jD * jMatrix[i][i+d]) / 2;
      }
    }
  }
  return E;
}

// Local energy at idx (chain already contains the trial state)
function localEnergy(chain, idx, params, jMatrix) {
  const { eM, eD, eF, jF, jFF, jD, rD, minRun } = params;
  const N = chain.length;
  const baseE = [eM, eD, eF];
  let E = baseE[chain[idx]];
  if (chain[idx] === STATES.F && fibrilRunLength(chain, idx) >= minRun) {
    // Short-range: direct neighbors
    if (idx > 0     && chain[idx-1] === STATES.F && fibrilRunLength(chain, idx-1) >= minRun) E -= jF;
    if (idx < N - 1 && chain[idx+1] === STATES.F && fibrilRunLength(chain, idx+1) >= minRun) E -= jF;
    // Long-range: any other active fibril site (|d| > 1)
    for (let j = 0; j < N; j++) {
      if (Math.abs(j - idx) > 1 && chain[j] === STATES.F && fibrilRunLength(chain, j) >= minRun) E -= jFF;
    }
  }
  if (chain[idx] === STATES.D) {
    for (let d = 1; d <= rD; d++) {
      if (idx - d >= 0 && chain[idx-d] === STATES.D) E -= jD * jMatrix[idx][idx-d];
      if (idx + d <  N && chain[idx+d] === STATES.D) E -= jD * jMatrix[idx][idx+d];
    }
  }
  return E;
}

function mcStep(chain, params, T, jMatrix, locked) {
  const N = chain.length;
  const c = [...chain];
  for (let _ = 0; _ < N; _++) {
    const idx = Math.floor(Math.random() * N);
    if (locked && locked[idx]) continue;          // irreversible: skip locked sites
    const oldState = c[idx];
    const newState = Math.floor(Math.random() * 3);
    if (newState === oldState) continue;
    const oldE = localEnergy(c, idx, params, jMatrix);
    c[idx] = newState;
    const newE = localEnergy(c, idx, params, jMatrix);
    const dE = newE - oldE;
    if (dE > 0 && Math.random() >= Math.exp(-dE / T)) c[idx] = oldState;
  }

  // After the sweep, lock any newly active fibril sites
  const newLocked = locked ? [...locked] : new Array(N).fill(false);
  if (locked !== null) {
    for (let i = 0; i < N; i++) {
      if (!newLocked[i] && c[i] === STATES.F && fibrilRunLength(c, i) >= params.minRun) {
        newLocked[i] = true;
      }
    }
  }

  return { chain: c, locked: newLocked };
}

function countStates(chain) {
  const counts = [0, 0, 0];
  chain.forEach(s => counts[s]++);
  return counts;
}

export default function App() {
  const defaultParams = { eM: 0, eD: 1.5, eF: 3.0, jF: 2.5, jFF: 0.5, jD: 1.2, rD: 3, minRun: 3 };
  const [nMonomers, setNMonomers] = useState(DEFAULT_N);
  const [params, setParams]        = useState(defaultParams);
  const [T, setT]                = useState(1.0);
  const [chain, setChain]        = useState(() => initChain(DEFAULT_N));
  const [jMatrix, setJMatrix]    = useState(() => buildJMatrix(DEFAULT_N));
  const [irreversible, setIrrev] = useState(false);
  const [locked, setLocked]      = useState(() => new Array(DEFAULT_N).fill(false));
  const [running, setRunning]    = useState(false);
  const [step, setStep]          = useState(0);
  const [history, setHistory]    = useState({ M: [], D: [], F: [], E: [] });
  const MAX_HIST = 200;

  const refs = {
    chain:        useRef(chain),
    params:       useRef(params),
    T:            useRef(T),
    running:      useRef(running),
    jMatrix:      useRef(jMatrix),
    irreversible: useRef(irreversible),
    locked:       useRef(locked),
  };
  refs.chain.current        = chain;
  refs.params.current       = params;
  refs.T.current            = T;
  refs.running.current      = running;
  refs.jMatrix.current      = jMatrix;
  refs.irreversible.current = irreversible;
  refs.locked.current       = locked;

  const animRef = useRef(null);
  const stepRef = useRef(0);

  const pushHistory = useCallback((c, p, jm) => {
    const ct = countStates(c);
    const E  = computeEnergy(c, p, jm);
    setHistory(prev => {
      const trim = a => a.length >= MAX_HIST ? a.slice(1) : a;
      return {
        M: [...trim(prev.M), ct[0] / c.length],
        D: [...trim(prev.D), ct[1] / c.length],
        F: [...trim(prev.F), ct[2] / c.length],
        E: [...trim(prev.E), E],
      };
    });
  }, []);

  const tick = useCallback(() => {
    if (!refs.running.current) return;
    const activeLocked = refs.irreversible.current ? refs.locked.current : null;
    const { chain: nc, locked: nl } = mcStep(refs.chain.current, refs.params.current, refs.T.current, refs.jMatrix.current, activeLocked);
    stepRef.current += 1;
    setChain(nc);
    if (refs.irreversible.current) setLocked(nl);
    setStep(stepRef.current);
    pushHistory(nc, refs.params.current, refs.jMatrix.current);
    animRef.current = requestAnimationFrame(tick);
  }, [pushHistory]);

  useEffect(() => {
    if (running) animRef.current = requestAnimationFrame(tick);
    else cancelAnimationFrame(animRef.current);
    return () => cancelAnimationFrame(animRef.current);
  }, [running, tick]);

  const reset = (n) => {
    const nn = n ?? nMonomers;
    setRunning(false);
    cancelAnimationFrame(animRef.current);
    const c = initChain(nn), jm = buildJMatrix(nn);
    setChain(c); setJMatrix(jm);
    setLocked(new Array(nn).fill(false));
    stepRef.current = 0; setStep(0);
    setHistory({ M: [], D: [], F: [], E: [] });
  };

  const doStep = () => {
    if (running) return;
    const activeLocked = irreversible ? refs.locked.current : null;
    const { chain: nc, locked: nl } = mcStep(refs.chain.current, refs.params.current, refs.T.current, refs.jMatrix.current, activeLocked);
    if (irreversible) setLocked(nl);
    stepRef.current += 1;
    setChain(nc); setStep(stepRef.current);
    pushHistory(nc, refs.params.current, refs.jMatrix.current);
  };

  const counts = countStates(chain);
  const E = computeEnergy(chain, params, jMatrix);

  const fibrilRuns = useMemo(() => {
    const n = chain.length;
    const runs = []; let i = 0;
    while (i < n) {
      if (chain[i] === STATES.F) {
        let j = i; while (j < n && chain[j] === STATES.F) j++;
        runs.push(j - i); i = j;
      } else i++;
    }
    return runs;
  }, [chain]);

  const activeFibrilCount = fibrilRuns.filter(r => r >= params.minRun).reduce((a, r) => a + r, 0);
  const activeFibrilFrac  = counts[2] > 0 ? activeFibrilCount / counts[2] : 0;

  const nField = { val: nMonomers, setVal: (v) => { setNMonomers(v); reset(v); } };

  const paramConfig = [
    { key: "eM",     label: "E_Monomer",      min: -2, max: 4,  step: 0.1, col: STATE_COLORS[0], desc: "Intrinsic cost — monomer" },
    { key: "eD",     label: "E_Disordered",   min: -2, max: 6,  step: 0.1, col: STATE_COLORS[1], desc: "Intrinsic cost — disordered" },
    { key: "eF",     label: "E_Fibril",       min: -2, max: 8,  step: 0.1, col: STATE_COLORS[2], desc: "Intrinsic cost — fibril" },
    { key: "jF",     label: "J_Fibril",       min: 0,  max: 8,  step: 0.1, col: "#f87171",       desc: "Fibril coupling (run≥minRun required)" },
    { key: "jD",     label: "J_D max",        min: 0,  max: 4,  step: 0.1, col: "#fb923c",       desc: "Max disordered coupling (quenched)" },
    { key: "rD",     label: "r_Disordered",   min: 1,  max: 10, step: 1,   col: "#fbbf24",       desc: "Range of D-D interactions" },
    { key: "minRun", label: "Min Fibril Run", min: 2,  max: 8,  step: 1,   col: "#c084fc",       desc: "Min run length to activate J_F" },
    { key: "jFF",    label: "J_Fibril long-range", min: 0, max: 4, step: 0.05, col: "#38bdf8",   desc: "Active fibril–fibril coupling (|d|>1, any distance)" },
  ];

  const Sparkline = ({ data, color, height = 38 }) => {
    if (data.length < 2) return <svg width="100%" height={height} />;
    const w = 260, h = height;
    const mn = Math.min(...data), mx = Math.max(...data), range = mx - mn || 1;
    const pts = data.map((v, i) => `${(i/(data.length-1))*w},${h - ((v-mn)/range)*(h-6) - 3}`).join(" ");
    return (
      <svg width="100%" viewBox={`0 0 ${w} ${h}`} preserveAspectRatio="none" style={{ display: "block" }}>
        <polyline points={pts} fill="none" stroke={color} strokeWidth="1.5" opacity="0.85" />
      </svg>
    );
  };

  return (
    <div style={{ background: "#0a0e1a", minHeight: "100vh", color: "#e2e8f0",
      fontFamily: "'IBM Plex Mono','Courier New',monospace", padding: 24, boxSizing: "border-box" }}>
      <style>{`
        @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@300;400;500;600&family=Space+Grotesk:wght@300;400;600&display=swap');
        input[type=range]{cursor:pointer;width:100%}
        ::-webkit-scrollbar{width:6px}::-webkit-scrollbar-track{background:#1e2536}
        ::-webkit-scrollbar-thumb{background:#4a5568;border-radius:3px}
      `}</style>

      <div style={{ marginBottom: 18 }}>
        <div style={{ fontSize: 10, letterSpacing: 4, color: "#7c3aed", textTransform: "uppercase", marginBottom: 3 }}>
          Monte Carlo · Metropolis Algorithm
        </div>
        <h1 style={{ margin: 0, fontSize: 20, fontWeight: 600, fontFamily: "'Space Grotesk',sans-serif", letterSpacing: -0.5 }}>
          1D Protein Aggregation — Ising Model
        </h1>
        <div style={{ fontSize: 10, color: "#475569", marginTop: 3 }}>
          Quenched random D-D couplings &nbsp;·&nbsp; Cooperative fibril threshold &nbsp;·&nbsp; N={chain.length}
        </div>
      </div>

      <div style={{ display: "grid", gridTemplateColumns: "1fr 268px", gap: 18 }}>

        {/* LEFT */}
        <div>
          {/* Chain visualization */}
          <div style={{ background: "#111827", border: "1px solid #1e2d4a", borderRadius: 8, padding: 14, marginBottom: 14 }}>
            <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 10 }}>
              <span style={{ fontSize: 10, color: "#64748b", letterSpacing: 2, textTransform: "uppercase" }}>
                Chain — step {step}
              </span>
              <span style={{ fontSize: 10, color: "#64748b" }}>
                {fibrilRuns.filter(r => r >= params.minRun).length} active fibril cluster(s)
              </span>
            </div>
            <div style={{ display: "flex", flexWrap: "wrap", gap: 2 }}>
              {chain.map((s, i) => {
                const isActiveF = s === STATES.F && fibrilRunLength(chain, i) >= params.minRun;
                const isLocked  = locked[i];
                return (
                  <div key={i} title={`${i}: ${STATE_NAMES[s]}${isLocked ? " [locked]" : ""}`} style={{
                    width: 10, height: 26, borderRadius: 2, flexShrink: 0,
                    background: STATE_COLORS[s],
                    opacity: s === STATES.F ? (isActiveF ? 1 : 0.28) : 0.82,
                    boxShadow: isLocked
                      ? `0 0 6px #a855f7, inset 0 0 3px #a855f744`
                      : isActiveF ? `0 0 5px ${STATE_COLORS[2]}99` : "none",
                    outline: isLocked ? "1px solid #a855f788" : "none",
                  }} />
                );
              })}
            </div>
            <div style={{ display: "flex", gap: 16, marginTop: 10, flexWrap: "wrap" }}>
              {STATE_NAMES.map((name, i) => (
                <div key={i} style={{ display: "flex", alignItems: "center", gap: 5, fontSize: 10 }}>
                  <div style={{ width: 9, height: 9, borderRadius: 2, background: STATE_COLORS[i] }} />
                  <span style={{ color: "#94a3b8" }}>{name}</span>
                </div>
              ))}
              <div style={{ display: "flex", alignItems: "center", gap: 5, fontSize: 10 }}>
                <div style={{ width: 9, height: 9, borderRadius: 2, background: STATE_COLORS[2], opacity: 0.28 }} />
                <span style={{ color: "#475569" }}>Fibril (sub-threshold)</span>
              </div>
              {irreversible && (
                <div style={{ display: "flex", alignItems: "center", gap: 5, fontSize: 10 }}>
                  <div style={{ width: 9, height: 9, borderRadius: "50%", background: "#a855f7", boxShadow: "0 0 4px #a855f7" }} />
                  <span style={{ color: "#c084fc" }}>Locked ({locked.filter(Boolean).length})</span>
                </div>
              )}
            </div>
          </div>

          {/* Stats */}
          <div style={{ display: "grid", gridTemplateColumns: "repeat(5,1fr)", gap: 8, marginBottom: 14 }}>
            {[
              { label: "Monomer",    val: `${(counts[0]/chain.length*100).toFixed(0)}%`, col: STATE_COLORS[0] },
              { label: "Disordered", val: `${(counts[1]/chain.length*100).toFixed(0)}%`, col: STATE_COLORS[1] },
              { label: "Fibril",     val: `${(counts[2]/chain.length*100).toFixed(0)}%`, col: STATE_COLORS[2] },
              { label: "Active F",   val: `${(activeFibrilFrac*100).toFixed(0)}%`, col: "#c084fc" },
              { label: "Energy",     val: E.toFixed(1), col: "#818cf8" },
            ].map(({ label, val, col }) => (
              <div key={label} style={{ background: "#111827", border: `1px solid ${col}33`, borderRadius: 8, padding: "8px 10px" }}>
                <div style={{ fontSize: 9, color: "#64748b", textTransform: "uppercase", letterSpacing: 1 }}>{label}</div>
                <div style={{ fontSize: 18, fontWeight: 600, color: col, marginTop: 1 }}>{val}</div>
              </div>
            ))}
          </div>

          {/* Sparklines */}
          <div style={{ background: "#111827", border: "1px solid #1e2d4a", borderRadius: 8, padding: 14, marginBottom: 14 }}>
            <div style={{ fontSize: 10, color: "#64748b", letterSpacing: 2, textTransform: "uppercase", marginBottom: 8 }}>
              Population Dynamics
            </div>
            {[
              { key: "M", color: STATE_COLORS[0], label: "Monomer" },
              { key: "D", color: STATE_COLORS[1], label: "Disordered" },
              { key: "F", color: STATE_COLORS[2], label: "Fibril" },
              { key: "E", color: "#818cf8",        label: "Energy" },
            ].map(({ key, color, label }) => (
              <div key={key} style={{ marginBottom: 5 }}>
                <div style={{ fontSize: 9, color: color, marginBottom: 1 }}>{label}</div>
                <Sparkline data={history[key]} color={color} />
              </div>
            ))}
          </div>

          {/* Buttons */}
          <div style={{ display: "flex", gap: 8, alignItems: "center", flexWrap: "wrap" }}>
            {[
              { label: running ? "⏸ Pause" : "▶ Run", fn: () => setRunning(r => !r),
                bg: running ? "#7f1d1d" : "#1e3a5f", bdr: running ? "#ef4444" : "#3b82f6", col: running ? "#fca5a5" : "#93c5fd" },
              { label: "Step",  fn: doStep, bg: "#1a1f35", bdr: "#374151", col: "#94a3b8", dis: running },
              { label: "Reset", fn: reset,  bg: "#1a1f35", bdr: "#374151", col: "#94a3b8" },
            ].map(({ label, fn, bg, bdr, col, dis }) => (
              <button key={label} onClick={fn} disabled={dis} style={{
                background: bg, border: `1px solid ${bdr}`, color: col,
                padding: "9px 20px", borderRadius: 6, cursor: dis ? "not-allowed" : "pointer",
                fontFamily: "inherit", fontSize: 11, letterSpacing: 1, textTransform: "uppercase",
                opacity: dis ? 0.4 : 1,
              }}>{label}</button>
            ))}

            {/* Irreversible fibril toggle */}
            <button
              onClick={() => setIrrev(v => !v)}
              title="When ON, fibril sites that reach the active run threshold can never revert"
              style={{
                background: irreversible ? "#2d1a4a" : "#1a1f35",
                border: `1px solid ${irreversible ? "#a855f7" : "#374151"}`,
                color: irreversible ? "#d8b4fe" : "#64748b",
                padding: "9px 14px", borderRadius: 6, cursor: "pointer",
                fontFamily: "inherit", fontSize: 11, letterSpacing: 1, textTransform: "uppercase",
                display: "flex", alignItems: "center", gap: 6,
              }}>
              <span style={{
                display: "inline-block", width: 10, height: 10, borderRadius: "50%",
                background: irreversible ? "#a855f7" : "#374151",
                boxShadow: irreversible ? "0 0 6px #a855f7" : "none",
                flexShrink: 0,
              }} />
              Irreversible F
            </button>
          </div>
        </div>

        {/* RIGHT */}
        <div>
          {/* N Monomers */}
          <div style={{ background: "#111827", border: "1px solid #1e2d4a", borderRadius: 8, padding: 14, marginBottom: 14 }}>
            <div style={{ fontSize: 10, color: "#64748b", letterSpacing: 2, textTransform: "uppercase", marginBottom: 8 }}>Chain Length</div>
            <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
              <input
                type="number" min={10} max={300} step={1}
                value={nMonomers}
                onChange={e => {
                  const v = Math.max(10, Math.min(300, parseInt(e.target.value) || DEFAULT_N));
                  setNMonomers(v);
                }}
                onBlur={() => reset(nMonomers)}
                onKeyDown={e => { if (e.key === "Enter") reset(nMonomers); }}
                style={{
                  background: "#0d1117", border: "1px solid #374151", color: "#e2e8f0",
                  borderRadius: 6, padding: "6px 10px", fontFamily: "inherit",
                  fontSize: 18, fontWeight: 600, width: 80, textAlign: "center",
                  outline: "none",
                }}
              />
              <span style={{ fontSize: 10, color: "#475569" }}>monomers<br/>(Enter or blur to apply)</span>
            </div>
            <input type="range" min={10} max={300} step={5} value={nMonomers}
              onChange={e => setNMonomers(parseInt(e.target.value))}
              onMouseUp={e => reset(parseInt(e.target.value))}
              onTouchEnd={e => reset(parseInt(e.target.value))}
              style={{ accentColor: "#64748b", marginTop: 8 }} />
            <div style={{ display: "flex", justifyContent: "space-between", fontSize: 9, color: "#374151", marginTop: 2 }}>
              <span>10</span><span>300</span>
            </div>
          </div>

          {/* Temperature */}
          <div style={{ background: "#111827", border: "1px solid #1e2d4a", borderRadius: 8, padding: 14, marginBottom: 14 }}>
            <div style={{ fontSize: 10, color: "#64748b", letterSpacing: 2, textTransform: "uppercase", marginBottom: 8 }}>Temperature</div>
            <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
              <span style={{ fontSize: 11, color: "#94a3b8" }}>k_B T</span>
              <span style={{ fontSize: 14, color: "#f0abfc", fontWeight: 600 }}>{T.toFixed(2)}</span>
            </div>
            <input type="range" min={0.1} max={5} step={0.05} value={T}
              onChange={e => setT(parseFloat(e.target.value))} style={{ accentColor: "#c084fc" }} />
            <div style={{ display: "flex", justifyContent: "space-between", fontSize: 9, color: "#475569", marginTop: 2 }}>
              <span>ordered</span><span>disordered</span>
            </div>
          </div>

          {/* Parameters */}
          <div style={{ background: "#111827", border: "1px solid #1e2d4a", borderRadius: 8, padding: 14 }}>
            <div style={{ fontSize: 10, color: "#64748b", letterSpacing: 2, textTransform: "uppercase", marginBottom: 12 }}>Parameters</div>
            {paramConfig.map(({ key, label, min, max, step: s, col, desc }) => (
              <div key={key} style={{ marginBottom: 12 }}>
                <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 2 }}>
                  <span style={{ fontSize: 11, color: col, fontWeight: 500 }}>{label}</span>
                  <span style={{ fontSize: 12, color: "#e2e8f0", fontWeight: 600 }}>
                    {(key === "rD" || key === "minRun") ? params[key] : params[key].toFixed(1)}
                  </span>
                </div>
                <input type="range" min={min} max={max} step={s} value={params[key]}
                  onChange={e => setParams(p => ({ ...p, [key]: parseFloat(e.target.value) }))}
                  style={{ accentColor: col }} />
                <div style={{ fontSize: 9, color: "#374151", marginTop: 1 }}>{desc}</div>
              </div>
            ))}
          </div>

          {/* Hamiltonian */}
          <div style={{ background: "#0d1117", border: "1px solid #1e2d4a", borderRadius: 8, padding: 12, marginTop: 14, fontSize: 10, lineHeight: 1.9 }}>
            <div style={{ color: "#64748b", marginBottom: 4, letterSpacing: 1, textTransform: "uppercase", fontSize: 9 }}>Hamiltonian</div>
            <div style={{ color: "#7dd3fc" }}>H = Σᵢ εₛᵢ</div>
            <div style={{ color: "#f87171", marginLeft: 8 }}>− J_F Σ⟨i,j⟩ δ(F,F)·𝟙[run≥{params.minRun}]</div>
            <div style={{ color: "#fb923c", marginLeft: 8 }}>− Σ|i−j|≤r J̃ᵢⱼ δ(D,D)</div>
            <div style={{ color: "#38bdf8", marginLeft: 8 }}>− J_FF Σ|i−j|&gt;1 δ(F*,F*)</div>
            <div style={{ marginTop: 8, borderTop: "1px solid #1e2d4a", paddingTop: 8, color: "#374151", fontSize: 9 }}>
              J̃ᵢⱼ ~ U[0, J_D max] — quenched, resampled on Reset
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
