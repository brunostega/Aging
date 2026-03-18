import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import {
  STATES, STATE_NAMES, STATE_COLORS,
  DEFAULT_N, DEFAULT_PARAMS, MIN_N, MAX_N,
  initChain, parseChain, countStates, fibrilRunLengths,
} from "./simulation.js";
import { computeEnergy, fibrilRunLength } from "./model.js";

// ── Sparkline ───────────────────────────────────────────────────────────────

const Y_TICKS = 3;   // number of y-axis tick lines (including min and max)
const X_TICKS = 4;   // number of x-axis tick marks on the bottom plot
const PAD_L   = 36;  // left padding for y-axis labels
const PAD_B   = 18;  // bottom padding for x-axis labels
const PAD_T   = 4;   // top padding
const PAD_R   = 4;   // right padding

function niceNum(v) {
  if (Math.abs(v) >= 1000) return v.toFixed(0);
  if (Math.abs(v) >= 10)   return v.toFixed(1);
  return v.toFixed(2);
}

function Sparkline({ data, color, height = 48, showXAxis = false, steps = null }) {
  if (data.length < 2) return <svg width="100%" height={height + (showXAxis ? PAD_B : 0)} />;

  const totalH = height + (showXAxis ? PAD_B : 0);
  const W = 260;
  const px0 = PAD_L, px1 = W - PAD_R;
  const py0 = PAD_T, py1 = height - PAD_T;
  const pw = px1 - px0, ph = py1 - py0;

  const mn = Math.min(...data), mx = Math.max(...data), range = mx - mn || 1;

  const xOf = i => px0 + (i / (data.length - 1)) * pw;
  const yOf = v => py1 - ((v - mn) / range) * ph;

  const pts = data.map((v, i) => `${xOf(i)},${yOf(v)}`).join(" ");

  // Y ticks
  const yTickVals = Array.from({ length: Y_TICKS }, (_, i) => mn + (i / (Y_TICKS - 1)) * range);

  // X ticks: evenly spaced sample indices, labelled with actual step numbers
  const xTickSampleIdxs = Array.from({ length: X_TICKS }, (_, i) =>
    Math.round(i * (data.length - 1) / (X_TICKS - 1))
  );

  return (
    <svg width="100%" viewBox={`0 0 ${W} ${totalH}`} preserveAspectRatio="none" style={{ display: "block" }}>
      {/* Y grid lines + labels */}
      {yTickVals.map((v, i) => {
        const y = yOf(v);
        return (
          <g key={i}>
            <line x1={px0} y1={y} x2={px1} y2={y} stroke="#1e2d4a" strokeWidth="0.5" />
            <text x={px0 - 3} y={y + 3} textAnchor="end" fontSize="7" fill="#475569">{niceNum(v)}</text>
          </g>
        );
      })}

      {/* X ticks + labels using actual step numbers */}
      {showXAxis && xTickSampleIdxs.map((si, i) => {
        const x = xOf(si);
        const label = steps ? steps[si] : si;
        return (
          <g key={i}>
            <line x1={x} y1={py1} x2={x} y2={py1 + 3} stroke="#475569" strokeWidth="0.5" />
            <text x={x} y={py1 + PAD_B - 3} textAnchor="middle" fontSize="7" fill="#475569">{label}</text>
          </g>
        );
      })}

      {/* Line */}
      <polyline points={pts} fill="none" stroke={color} strokeWidth="1.5" opacity="0.85" />
    </svg>
  );
}

// ── Parameter slider config ─────────────────────────────────────────────────

const PARAM_CONFIG = [
  { key: "eM",     label: "E_Monomer",           min: 0, max: 8, step: 0.1, col: STATE_COLORS[0], desc: "Intrinsic cost — monomer" },
  { key: "eD",     label: "E_Disordered",        min: 0, max: 8, step: 0.1, col: STATE_COLORS[1], desc: "Intrinsic cost — disordered" },
  { key: "eF",     label: "E_Fibril",            min: 0, max: 8, step: 0.1, col: STATE_COLORS[2], desc: "Intrinsic cost — fibril" },
  { key: "jD",     label: "J_Disordered",        min: 0, max: 8, step: 0.1, col: "#E69F00",       desc: "Disordered–disordered nearest-neighbour coupling" },
  { key: "jF",     label: "J_Fibril",            min: 0, max: 8, step: 0.1, col: "#CC79A7",       desc: "Fibril coupling (run≥minRun required)" },
  { key: "hFF",    label: "h_FF background",     min: 0, max: 8, step: 0.1, col: "#38bdf8",       desc: "Background field on all F sites when an active run exists" },
  { key: "minRun", label: "Min Fibril Run",      min: 2, max: 8, step: 1,   col: "#c084fc",       desc: "Min run length to activate J_F" },
];


const MAX_HIST     = 400;  // max samples kept per series
const PRUNE_TARGET = 200;  // downsample to this size when MAX_HIST is reached

// Uniformly downsample array a to targetLen by picking evenly spaced indices.
function downsample(a, targetLen) {
  if (a.length <= targetLen) return a;
  return Array.from({ length: targetLen }, (_, i) =>
    a[Math.round(i * (a.length - 1) / (targetLen - 1))]
  );
}

// ── App ─────────────────────────────────────────────────────────────────────

export default function App() {
  const [params, setParams]               = useState(DEFAULT_PARAMS);
  const [T, setT]                         = useState(1.0);
  const [chain, setChain]                 = useState(() => initChain(DEFAULT_N));
  const [irreversible, setIrrev]          = useState(false);
  const [stopOnFibril, setStopOnFibril]   = useState(true);
  const [locked, setLocked]               = useState(() => new Array(DEFAULT_N).fill(false));
  const [running, setRunning]             = useState(false);
  const [step, setStep]                   = useState(0);
  const [history, setHistory]             = useState({ M: [], D: [], F: [], E: [], steps: [] });
  const [snapCount, setSnapCount]         = useState(0);
  const [seqInput, setSeqInput]           = useState("");
  const [seqError, setSeqError]           = useState(null);
  const [showChain, setShowChain]         = useState(true);
  const SWEEPS_PER_BATCH = 10;

  // Refs for values needed inside the rAF loop without stale closures
  const refs = {
    chain:        useRef(chain),
    params:       useRef(params),
    T:            useRef(T),
    running:      useRef(running),
    irreversible: useRef(irreversible),
    locked:       useRef(locked),
    stopOnFibril: useRef(stopOnFibril),
  };
  refs.chain.current        = chain;
  refs.params.current       = params;
  refs.T.current            = T;
  refs.running.current      = running;
  refs.irreversible.current = irreversible;
  refs.locked.current       = locked;
  refs.stopOnFibril.current = stopOnFibril;

  // Single source of truth for chain length (uncontrolled input pattern)
  const sizeRef        = useRef(DEFAULT_N);
  const inputRef       = useRef(null);
  const energyRef      = useRef(computeEnergy(initChain(DEFAULT_N), DEFAULT_PARAMS));
  const workerRef      = useRef(null);
  const stepRef        = useRef(0);
  const energySumRef   = useRef(0);
  const mSumRef        = useRef(0);
  const dSumRef        = useRef(0);
  const fSumRef        = useRef(0);
  const generationRef  = useRef(0);  // incremented on reset to discard stale batches
  const trajectoryRef  = useRef([]);

  // ── history / trajectory ──────────────────────────────────────────────────

  const pushHistory = useCallback((c, E) => {
    const ct = countStates(c);
    // averages accumulated in onmessage for per-snapshot accuracy
    const avgE = energySumRef.current / stepRef.current;
    trajectoryRef.current.push({ step: stepRef.current, chain: [...c], E });
    setSnapCount(trajectoryRef.current.length);
    setHistory(prev => {
      const appendAndPrune = (a, v) => {
        const next = [...a, v];
        return next.length > MAX_HIST ? downsample(next, PRUNE_TARGET) : next;
      };
      return {
        M: appendAndPrune(prev.M, ct[0] / c.length),
        D: appendAndPrune(prev.D, ct[1] / c.length),
        F: appendAndPrune(prev.F, ct[2] / c.length),
        E: appendAndPrune(prev.E, E),
        steps: appendAndPrune(prev.steps, stepRef.current),
      };
    });
  }, []);

  // ── worker setup ──────────────────────────────────────────────────────────

  // Create worker once on mount
  useEffect(() => {
    workerRef.current = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
    const gen = generationRef; // capture ref
    workerRef.current.onmessage = (e) => {
      const msg = e.data;
      if (msg.type === "batch") {
        // Discard batches from a previous generation (i.e. before last reset)
        if (msg.generation !== gen.current) return;
        const { chain: nc, locked: nl, E: newE, step, snapshots } = msg;
        energyRef.current = newE;
        stepRef.current = step;
        setChain(nc);
        setLocked(nl);
        setStep(step);
        if (snapshots.length > 0) {
          const last = snapshots[snapshots.length - 1];
          // Accumulate averages locally — don't trust worker's energySum
          snapshots.forEach(snap => {
            energySumRef.current += snap.E;
            const ct = snap.chain.reduce((a, s) => { a[s]++; return a; }, [0,0,0]);
            const N = snap.chain.length;
            mSumRef.current += ct[0] / N;
            dSumRef.current += ct[1] / N;
            fSumRef.current += ct[2] / N;
          });
          pushHistory(last.chain, last.E);
          for (const snap of snapshots) {
            trajectoryRef.current.push(snap);
          }
          setSnapCount(trajectoryRef.current.length);
        }
      }
      if (msg.type === "done") {
        if (msg.generation !== gen.current) return;
        setRunning(false);
      }
    };
    return () => workerRef.current.terminate();
  }, [pushHistory]);

  // Start worker when running becomes true.
  // Stop is sent directly by the pause button, reset, and applySequence.
  useEffect(() => {
    if (!workerRef.current) return;
    if (running) {
      workerRef.current.postMessage({
        type: "start",
        chain:         refs.chain.current,
        locked:        refs.locked.current,
        params:        refs.params.current,
        T:             refs.T.current,
        irreversible:  refs.irreversible.current,
        stopOnFibril:  refs.stopOnFibril.current,
        E:             energyRef.current,
        step:          stepRef.current,
        energySum:     energySumRef.current,
        sweepsPerBatch: SWEEPS_PER_BATCH,
        generation:     generationRef.current,
      });
    }
  }, [running]);

  // When params or T change, resync energy and inform worker
  useEffect(() => {
    energyRef.current = computeEnergy(refs.chain.current, refs.params.current);
    // Push updated params to the worker so changes take effect immediately
    if (workerRef.current) {
      workerRef.current.postMessage({
        type: "params",
        params: refs.params.current,
        T:      refs.T.current,
        E:      energyRef.current,
      });
    }
  }, [params, T]);



  // ── reset ─────────────────────────────────────────────────────────────────

  const reset = (n) => {
    const nn = (n !== undefined && !isNaN(n))
      ? Math.max(MIN_N, Math.min(MAX_N, n))
      : sizeRef.current;
    sizeRef.current = nn;
    if (inputRef.current) inputRef.current.value = nn;
    generationRef.current += 1;  // invalidate any in-flight batches
    setRunning(false);
    if (workerRef.current) workerRef.current.postMessage({ type: "stop" });
    const freshChain = initChain(nn);
    energyRef.current = computeEnergy(freshChain, params);
    setChain(freshChain);
    setLocked(new Array(nn).fill(false));
    stepRef.current = 0; setStep(0);
    energySumRef.current = 0;
    mSumRef.current = 0; dSumRef.current = 0; fSumRef.current = 0;
    trajectoryRef.current = [];
    setSnapCount(0);
    setHistory({ M: [], D: [], F: [], E: [], steps: [] });
  };

  // ── sequence apply ───────────────────────────────────────────────────────

  const applySequence = () => {
    const { chain: parsed, error } = parseChain(seqInput);
    if (error) { setSeqError(error); return; }
    setSeqError(null);
    const nn = parsed.length;
    sizeRef.current = nn;
    if (inputRef.current) inputRef.current.value = nn;
    generationRef.current += 1;
    setRunning(false);
    if (workerRef.current) workerRef.current.postMessage({ type: "stop" });
    const freshChain = parsed;
    energyRef.current = computeEnergy(freshChain, params);
    setChain(freshChain);
    setLocked(new Array(nn).fill(false));
    stepRef.current = 0; setStep(0);
    energySumRef.current = 0;
    mSumRef.current = 0; dSumRef.current = 0; fSumRef.current = 0;
    trajectoryRef.current = [];
    setSnapCount(0);
    setHistory({ M: [], D: [], F: [], E: [], steps: [] });
  };

  // ── single step ───────────────────────────────────────────────────────────

  const doStep = () => {
    if (running) return;
    workerRef.current.postMessage({
      type:         "step",
      chain:        refs.chain.current,
      locked:       refs.locked.current,
      params:       refs.params.current,
      T:            refs.T.current,
      irreversible: refs.irreversible.current,
      E:            energyRef.current,
      step:         stepRef.current,
      energySum:    energySumRef.current,
    });
  };

  // ── trajectory export ─────────────────────────────────────────────────────

  const saveTrajectory = () => {
    const payload = {
      metadata: {
        created:    new Date().toISOString(),
        software:   "Aging",
        n:          chain.length,
        totalSteps: stepRef.current,
        snapshots:  trajectoryRef.current.length,
        params:     { ...params, T, irreversible },
        states:     { 0: "Monomer", 1: "Disordered", 2: "Fibril" },
      },
      trajectory: trajectoryRef.current,
    };
    const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement("a");
    a.href     = url;
    a.download = `aging_trajectory_${Date.now()}.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  // ── derived display values ────────────────────────────────────────────────

  const counts          = countStates(chain);
  const E               = computeEnergy(chain, params);
  const fibrilRuns      = useMemo(() => fibrilRunLengths(chain), [chain]);
  const activeFibrilCount = fibrilRuns.filter(r => r >= params.minRun).reduce((a, r) => a + r, 0);
  const activeFibrilFrac  = counts[2] > 0 ? activeFibrilCount / counts[2] : 0;

  // ── render ────────────────────────────────────────────────────────────────

  return (
    <div style={{
      background: "#0a0e1a", minHeight: "100vh", color: "#e2e8f0",
      fontFamily: "'IBM Plex Mono','Courier New',monospace", padding: 24, boxSizing: "border-box",
    }}>
      <style>{`
        @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@300;400;500;600&family=Space+Grotesk:wght@300;400;600&display=swap');
        input[type=range]{cursor:pointer;width:100%}
        ::-webkit-scrollbar{width:6px}::-webkit-scrollbar-track{background:#1e2536}
        ::-webkit-scrollbar-thumb{background:#4a5568;border-radius:3px}
      `}</style>

      {/* Header */}
      <div style={{ marginBottom: 18, display: "flex", alignItems: "center", gap: 16 }}>
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 220 160" style={{ height: 52, flexShrink: 0 }}>
          <line x1="20" y1="80" x2="200" y2="80" stroke="#475569" strokeWidth="1.5" opacity="0.35"/>
          <circle cx="40" cy="80" r="14" fill="#CC79A7" opacity="0.9"/>
          <line x1="40" y1="66" x2="40" y2="42" stroke="#CC79A7" strokeWidth="2.5"/>
          <polygon points="40,36 35,46 45,46" fill="#CC79A7"/>
          <circle cx="88" cy="80" r="14" fill="#CC79A7" opacity="0.9"/>
          <line x1="88" y1="66" x2="88" y2="42" stroke="#CC79A7" strokeWidth="2.5"/>
          <polygon points="88,36 83,46 93,46" fill="#CC79A7"/>
          <circle cx="136" cy="80" r="18" fill="none" stroke="#CC79A7" strokeWidth="1" opacity="0.35"/>
          <circle cx="136" cy="80" r="14" fill="#CC79A7" opacity="0.9"/>
          <line x1="136" y1="66" x2="136" y2="42" stroke="#CC79A7" strokeWidth="2.5"/>
          <polygon points="136,36 131,46 141,46" fill="#CC79A7"/>
          <path d="M40,58 Q88,28 136,58" fill="none" stroke="#CC79A7" strokeWidth="1" strokeDasharray="3 3" opacity="0.55"/>
          <text x="88" y="22" textAnchor="middle" fontFamily="'IBM Plex Mono', monospace" fontSize="11" fill="#CC79A7" opacity="0.9">J_F</text>
          <circle cx="176" cy="80" r="14" fill="#E69F00" opacity="0.85"/>
          <line x1="169" y1="88" x2="155" y2="62" stroke="#E69F00" strokeWidth="2.5"/>
          <polygon points="152,57 159,65 165,57" fill="#E69F00" transform="rotate(-30,159,61)"/>
          <text x="40"  y="112" textAnchor="middle" fontFamily="'IBM Plex Mono', monospace" fontSize="11" fill="#CC79A7">F</text>
          <text x="88"  y="112" textAnchor="middle" fontFamily="'IBM Plex Mono', monospace" fontSize="11" fill="#CC79A7">F</text>
          <text x="136" y="112" textAnchor="middle" fontFamily="'IBM Plex Mono', monospace" fontSize="11" fill="#CC79A7">F</text>
          <text x="176" y="112" textAnchor="middle" fontFamily="'IBM Plex Mono', monospace" fontSize="11" fill="#E69F00">D</text>
        </svg>
        <div>
          <h1 style={{ margin: 0, fontSize: 22, fontWeight: 600, fontFamily: "'Space Grotesk',sans-serif", letterSpacing: -0.5 }}>
            <span style={{ color: "#e2e8f0" }}>Ag</span><span style={{ color: "#7c3aed" }}>I</span><span style={{ color: "#e2e8f0" }}>ng</span>
          </h1>
          <div style={{ fontSize: 11, color: "#94a3b8", marginTop: 3 }}>
            Modelling protein aggregation using a 1D Ising model
          </div>
          <div style={{ fontSize: 10, color: "#475569", marginTop: 2 }}>
            Monte Carlo · Metropolis &nbsp;·&nbsp; N={chain.length}
          </div>
        </div>
      </div>

      <div style={{ display: "grid", gridTemplateColumns: "1fr 268px", gap: 18 }}>

        {/* LEFT */}
        <div>
          {/* Chain visualization */}
          <div style={{ background: "#111827", border: "1px solid #1e2d4a", borderRadius: 8, padding: 14, marginBottom: 14 }}>
            <div style={{ display: "flex", justifyContent: "space-between", marginBottom: showChain ? 10 : 0 }}>
              <span style={{ fontSize: 10, color: "#64748b", letterSpacing: 2, textTransform: "uppercase" }}>
                Chain — step {step}
              </span>
              <div style={{ display: "flex", gap: 12, alignItems: "center" }}>
                <span style={{ fontSize: 10, color: "#64748b" }}>
                  {fibrilRuns.filter(r => r >= params.minRun).length} active fibril cluster(s)
                </span>
                <button onClick={() => setShowChain(v => !v)} style={{
                  background: "none", border: "none", color: "#475569", cursor: "pointer",
                  fontSize: 10, padding: 0, fontFamily: "inherit", letterSpacing: 1,
                  textTransform: "uppercase",
                }}>{showChain ? "▲ hide" : "▼ show"}</button>
              </div>
            </div>
            {showChain && <div style={{ display: "flex", flexWrap: "wrap", gap: 2 }}>
              {chain.map((s, i) => {
                const isActiveF = s === STATES.F && fibrilRunLength(chain, i) >= params.minRun;
                const isLocked  = locked[i];
                return (
                  <div key={i} title={`${i}: ${STATE_NAMES[s]}${isLocked ? " [locked]" : ""}`} style={{
                    width: 10, height: 26, borderRadius: 2, flexShrink: 0,
                    background: STATE_COLORS[s],
                    opacity: s === STATES.F ? (isActiveF ? 1 : 0.28) : 0.82,
                    boxShadow: isLocked
                      ? "0 0 6px #a855f7, inset 0 0 3px #a855f744"
                      : isActiveF ? `0 0 5px ${STATE_COLORS[2]}99` : "none",
                    outline: isLocked ? "1px solid #a855f788" : "none",
                  }} />
                );
              })}
            </div>}
            {showChain && <div style={{ display: "flex", gap: 16, marginTop: 10, flexWrap: "wrap" }}>
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
            </div>}
          </div>

          {/* Stats */}
          <div style={{ display: "grid", gridTemplateColumns: "repeat(5,1fr)", gap: 8, marginBottom: 14 }}>
            {[
              { label: "⟨M⟩",        val: step > 0 ? `${(mSumRef.current / step * 100).toFixed(0)}%` : "—", col: STATE_COLORS[0] },
              { label: "⟨D⟩",        val: step > 0 ? `${(dSumRef.current / step * 100).toFixed(0)}%` : "—", col: STATE_COLORS[1] },
              { label: "⟨F⟩",        val: step > 0 ? `${(fSumRef.current / step * 100).toFixed(0)}%` : "—", col: STATE_COLORS[2] },
              { label: "Active F",   val: `${(activeFibrilFrac * 100).toFixed(0)}%`,          col: "#c084fc" },
              { label: "Avg Energy", val: step > 0 ? (energySumRef.current / step).toFixed(1) : "—",  col: "#818cf8" },
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
              { key: "M", color: STATE_COLORS[0], label: "Monomer",         xAxis: false },
              { key: "D", color: STATE_COLORS[1], label: "Disordered",      xAxis: false },
              { key: "F", color: STATE_COLORS[2], label: "Fibril",          xAxis: false },
              { key: "E", color: "#818cf8",       label: "Energy",          xAxis: true  },
            ].map(({ key, color, label, xAxis }) => (
              <div key={key} style={{ marginBottom: 5 }}>
                <div style={{ fontSize: 9, color: color, marginBottom: 1 }}>{label}</div>
                <Sparkline data={history[key]} color={color} showXAxis={xAxis} steps={xAxis ? history.steps : null} />
              </div>
            ))}
          </div>

          {/* Buttons */}
          <div style={{ display: "flex", gap: 8, alignItems: "center", flexWrap: "wrap" }}>
            {[
              { label: running ? "⏸ Pause" : "▶ Run", fn: () => {
                  if (running && workerRef.current) workerRef.current.postMessage({ type: "stop" });
                  setRunning(r => !r);
                },
                bg: running ? "#7f1d1d" : "#1e3a5f", bdr: running ? "#ef4444" : "#3b82f6", col: running ? "#fca5a5" : "#93c5fd" },
              { label: "Step",  fn: doStep,         bg: "#1a1f35", bdr: "#374151", col: "#94a3b8", dis: running },
              { label: "Reset", fn: reset,          bg: "#1a1f35", bdr: "#374151", col: "#94a3b8" },
              { label: `Save (${snapCount})`, fn: saveTrajectory, bg: "#1a2a1a", bdr: "#56B4E9", col: "#93d6f5", dis: snapCount === 0 },
            ].map(({ label, fn, bg, bdr, col, dis }) => (
              <button key={label} onClick={fn} disabled={dis} style={{
                background: bg, border: `1px solid ${bdr}`, color: col,
                padding: "9px 20px", borderRadius: 6, cursor: dis ? "not-allowed" : "pointer",
                fontFamily: "inherit", fontSize: 11, letterSpacing: 1, textTransform: "uppercase",
                opacity: dis ? 0.4 : 1,
              }}>{label}</button>
            ))}

            {/* Stop on full fibril toggle */}
            <button
              onClick={() => setStopOnFibril(v => !v)}
              title="Stop simulation when chain reaches 100% fibril"
              style={{
                background: stopOnFibril ? "#1a2a1a" : "#1a1f35",
                border: `1px solid ${stopOnFibril ? "#56B4E9" : "#374151"}`,
                color: stopOnFibril ? "#93d6f5" : "#64748b",
                padding: "9px 14px", borderRadius: 6, cursor: "pointer",
                fontFamily: "inherit", fontSize: 11, letterSpacing: 1, textTransform: "uppercase",
                display: "flex", alignItems: "center", gap: 6,
              }}>
              <span style={{
                display: "inline-block", width: 10, height: 10, borderRadius: "50%",
                background: stopOnFibril ? "#56B4E9" : "#374151",
                boxShadow: stopOnFibril ? "0 0 6px #56B4E9" : "none",
                flexShrink: 0,
              }} />
              Stop at 100% F
            </button>

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
          {/* Chain Length */}
          <div style={{ background: "#111827", border: "1px solid #1e2d4a", borderRadius: 8, padding: 14, marginBottom: 14 }}>
            <div style={{ fontSize: 10, color: "#64748b", letterSpacing: 2, textTransform: "uppercase", marginBottom: 8 }}>Chain Length</div>
            <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
              <input
                ref={inputRef}
                type="number" min={MIN_N} max={MAX_N} step={1}
                defaultValue={DEFAULT_N}
                onKeyDown={e => {
                  if (e.key === "Enter") {
                    const v = parseInt(e.target.value);
                    reset(isNaN(v) ? sizeRef.current : v);
                    e.target.blur();
                  } else if (e.key === "Escape") {
                    e.target.value = sizeRef.current;
                    e.target.blur();
                  }
                }}
                style={{
                  background: "#0d1117", border: "1px solid #374151", color: "#e2e8f0",
                  borderRadius: 6, padding: "6px 10px", fontFamily: "inherit",
                  fontSize: 18, fontWeight: 600, width: 80, textAlign: "center", outline: "none",
                }}
              />
              <span style={{ fontSize: 10, color: "#475569" }}>sites<br/>(Enter to apply)</span>
            </div>
            <input type="range" min={MIN_N} max={MAX_N} step={5} defaultValue={DEFAULT_N}
              onChange={e => { if (inputRef.current) inputRef.current.value = e.target.value; }}
              onMouseUp={e => reset(parseInt(e.target.value))}
              onTouchEnd={e => reset(parseInt(e.target.value))}
              style={{ accentColor: "#64748b", marginTop: 8 }} />
            <div style={{ display: "flex", justifyContent: "space-between", fontSize: 9, color: "#374151", marginTop: 2 }}>
              <span>{MIN_N}</span><span>{MAX_N}</span>
            </div>
          </div>

          {/* Sequence input */}
          <div style={{ background: "#111827", border: `1px solid ${seqError ? "#ef444466" : "#1e2d4a"}`, borderRadius: 8, padding: 14, marginBottom: 14 }}>
            <div style={{ fontSize: 10, color: "#64748b", letterSpacing: 2, textTransform: "uppercase", marginBottom: 8 }}>
              Starting conformation
            </div>
            <div style={{ display: "flex", gap: 6 }}>
              <input
                type="text"
                value={seqInput}
                onChange={e => { setSeqInput(e.target.value); setSeqError(null); }}
                onKeyDown={e => { if (e.key === "Enter") applySequence(); }}
                placeholder="e.g. MMMFFFDDD"
                spellCheck={false}
                style={{
                  flex: 1, background: "#0d1117", border: `1px solid ${seqError ? "#ef4444" : "#374151"}`,
                  color: "#e2e8f0", borderRadius: 6, padding: "6px 8px",
                  fontFamily: "inherit", fontSize: 11, outline: "none",
                  letterSpacing: 1, textTransform: "uppercase",
                }}
              />
              <button onClick={applySequence} style={{
                background: "#1a1f35", border: "1px solid #374151", color: "#94a3b8",
                borderRadius: 6, padding: "6px 12px", cursor: "pointer",
                fontFamily: "inherit", fontSize: 11, letterSpacing: 1, textTransform: "uppercase",
                flexShrink: 0,
              }}>Apply</button>
            </div>
            {seqError
              ? <div style={{ fontSize: 9, color: "#f87171", marginTop: 5 }}>{seqError}</div>
              : <div style={{ fontSize: 9, color: "#374151", marginTop: 5 }}>M · D · F — Enter or Apply</div>
            }
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
            {PARAM_CONFIG.map(({ key, label, min, max, step: s, col, desc }) => (
              <div key={key} style={{ marginBottom: 12 }}>
                <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 2 }}>
                  <span style={{ fontSize: 11, color: col, fontWeight: 500 }}>{label}</span>
                  <span style={{ fontSize: 12, color: "#e2e8f0", fontWeight: 600 }}>
                    {key === "minRun" ? params[key] : params[key].toFixed(1)}
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
            <div style={{ color: "#CC79A7", marginLeft: 8 }}>− J_F Σ⟨i,j⟩ δ(F,F)·𝟙[run≥{params.minRun}]</div>
            <div style={{ color: "#E69F00", marginLeft: 8 }}>− J_D Σ⟨i,j⟩ δ(D,D)</div>
            <div style={{ color: "#38bdf8", marginLeft: 8 }}>− h_FF · N_F · 𝟙[∃ active run]</div>
          </div>
        </div>

      </div>
    </div>
  );
}
