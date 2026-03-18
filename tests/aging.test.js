import { describe, it, expect } from "vitest";
import { STATES, fibrilRunLength, computeEnergy,
         buildRunLen, countActiveRuns, buildSimState, updateRunLen,
         deltaEnergyFast } from "../src/model.js";
import { initChain, countStates, fibrilRunLengths, mcStep, DEFAULT_PARAMS } from "../src/simulation.js";

const { M, D, F } = STATES;

// ── shared parameter set ────────────────────────────────────────────────────
const P = { eM: 1, eD: 2, eF: 4, jF: 3, hFF: 1, jD: 2, minRun: 3 };

// ── fibrilRunLength ─────────────────────────────────────────────────────────
describe("fibrilRunLength", () => {
  it("returns 0 for a monomer site", () => {
    expect(fibrilRunLength([M, M, M], 1)).toBe(0);
  });

  it("returns 0 for a disordered site", () => {
    expect(fibrilRunLength([D, D, D], 1)).toBe(0);
  });

  it("returns 1 for an isolated fibril", () => {
    expect(fibrilRunLength([M, F, M], 1)).toBe(1);
  });

  it("returns full run length for a contiguous block", () => {
    // chain: M F F F M  — run of 3 at positions 1,2,3
    const chain = [M, F, F, F, M];
    expect(fibrilRunLength(chain, 1)).toBe(3);
    expect(fibrilRunLength(chain, 2)).toBe(3);
    expect(fibrilRunLength(chain, 3)).toBe(3);
  });

  it("does not bleed across non-fibril sites", () => {
    // chain: F F M F F — two runs of 2
    const chain = [F, F, M, F, F];
    expect(fibrilRunLength(chain, 0)).toBe(2);
    expect(fibrilRunLength(chain, 4)).toBe(2);
  });
});

// ── computeEnergy: pure intrinsic (no interactions) ────────────────────────
describe("computeEnergy — intrinsic only", () => {
  it("all monomers: E = N × eM", () => {
    const N = 5;
    const chain = new Array(N).fill(M);
    // No interactions for monomers
    expect(computeEnergy(chain, P)).toBeCloseTo(N * P.eM);
  });

  it("single monomer: E = eM", () => {
    expect(computeEnergy([M], P)).toBeCloseTo(P.eM);
  });

  it("single disordered: E = eD (no neighbours)", () => {
    expect(computeEnergy([D], P)).toBeCloseTo(P.eD);
  });

  it("single fibril: E = eF, run=1 < minRun so no coupling", () => {
    expect(computeEnergy([F], P)).toBeCloseTo(P.eF);
  });
});

// ── computeEnergy: disordered interactions ─────────────────────────────────
describe("computeEnergy — disordered coupling", () => {
  it("two adjacent D sites: E = 2×eD − jD", () => {
    // one D-D pair → −jD (counted once with /2 per site × 2 sites)
    const expected = 2 * P.eD - P.jD;
    expect(computeEnergy([D, D], P)).toBeCloseTo(expected);
  });

  it("three adjacent D sites: E = 3×eD − 2×jD", () => {
    // two D-D pairs
    const expected = 3 * P.eD - 2 * P.jD;
    expect(computeEnergy([D, D, D], P)).toBeCloseTo(expected);
  });

  it("two D sites separated by M: no interaction", () => {
    const expected = 2 * P.eD + P.eM;
    expect(computeEnergy([D, M, D], P)).toBeCloseTo(expected);
  });
});

// ── computeEnergy: fibril coupling ─────────────────────────────────────────
describe("computeEnergy — fibril coupling", () => {
  it("fibril run < minRun: no coupling, no hFF (no active run)", () => {
    // run of 2, minRun=3 → no active run exists, no jF, no hFF
    const chain = [F, F];
    const expected = 2 * P.eF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("fibril run exactly minRun=3: 2 nn pairs + hFF background on all 3 F sites", () => {
    // [F,F,F]: hasActive=true, nF=3 → −hFF×3
    const chain = [F, F, F];
    const expected = 3 * P.eF - 2 * P.jF - P.hFF * 3;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("fibril run of 4: 3 nn pairs + hFF background on all 4 F sites", () => {
    // [F,F,F,F]: hasActive=true, nF=4
    const chain = [F, F, F, F];
    const expected = 4 * P.eF - 3 * P.jF - P.hFF * 4;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("fibril run of 5: 4 nn pairs + hFF background on all 5 F sites", () => {
    // [F,F,F,F,F]: hasActive=true, nF=5
    const chain = [F, F, F, F, F];
    const expected = 5 * P.eF - 4 * P.jF - P.hFF * 5;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("two active runs: nn jF within each run + hFF background on all 6 F sites", () => {
    // [F,F,F,M,F,F,F]: run1=0,1,2  run2=4,5,6  nF=6
    // nn: (0,1),(1,2),(4,5),(5,6) = 4×jF; hFF acts on all 6 F sites
    const chain = [F, F, F, M, F, F, F];
    const expected = 6 * P.eF + P.eM - 4 * P.jF - P.hFF * 6;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });
});

// ── computeEnergy: mixed states ────────────────────────────────────────────
describe("computeEnergy — mixed states", () => {
  it("[F,F,F,M,M] with minRun=3: active run, hFF on all 3 F sites", () => {
    // fibril block [0,1,2]: active, nF=3
    const chain = [F, F, F, M, M];
    const expected = 3 * P.eF + 2 * P.eM - 2 * P.jF - P.hFF * 3;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("[M,D,D,F,F,F] with minRun=3: D-D, F-F nn, and hFF background", () => {
    // D: 1 nn pair → −jD; F: run=3, nn=2 → −2×jF; hFF on 3 F sites
    const chain = [M, D, D, F, F, F];
    const expected = P.eM + 2*P.eD + 3*P.eF - P.jD - 2*P.jF - P.hFF * 3;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });
});

describe("computeEnergy — hFF background field", () => {
  it("hFF does not fire when no active run exists", () => {
    // [F,F,M,F]: runs of 2 and 1, minRun=3 — no active run, no hFF
    const chain = [F, F, M, F];
    const expected = 3 * P.eF + P.eM;  
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("hFF fires on ALL F sites including sub-threshold ones once an active run exists", () => {
    // [F,F,F,M,F,F]: run1=3 (active), run2=2 (sub-threshold), nF=5
    // nn within run1: (0,1),(1,2)=2×jF; no nn in run2 (below minRun, no jF)
    // hFF acts on all 5 F sites
    const chain = [F, F, F, M, F, F];
    const expected = 5 * P.eF + P.eM - 2 * P.jF - P.hFF * 5;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });
});

// ── deltaEnergy consistency ───────────────────────────────────────────────
// For every flip, deltaEnergy must equal computeEnergy(new) - computeEnergy(old).

function checkDelta(chain, idx, newState, params) {
  const oldState = chain[idx];
  const { nF, nActiveRuns, runLen } = buildSimState(chain, params.minRun);
  const oldE = computeEnergy(chain, params);
  chain[idx] = newState;
  updateRunLen(chain, runLen, idx);
  const dE = deltaEnergyFast(chain, idx, oldState, newState, params, nF, nActiveRuns, runLen);
  const newE = computeEnergy(chain, params);
  chain[idx] = oldState; // restore
  return { dE, expected: newE - oldE };
}

describe("deltaEnergy — exact consistency with computeEnergy", () => {
  const baseCases = [
    { desc: "M→F inside existing active run",  chain: [F,F,F,M,F,F,F], idx: 3, to: F },
    { desc: "F→M breaking an active run",      chain: [F,F,F,F,F],     idx: 2, to: M },
    { desc: "M→D with no interactions",        chain: [M,M,M,M,M],     idx: 2, to: D },
    { desc: "D→M inside D block",              chain: [D,D,D,M,F,F,F], idx: 1, to: M },
    { desc: "M→F creating first active run",   chain: [F,F,M,M,M],     idx: 2, to: F },
    { desc: "F→M removing only active run",    chain: [F,F,F,M,M],     idx: 1, to: M },
    { desc: "M→F next to sub-threshold run",   chain: [F,F,M,M,M],     idx: 2, to: F },
    { desc: "F→D inside active run",           chain: [F,F,F,F,F],     idx: 3, to: D },
    { desc: "F→D inside active run",           chain: [M,M,M,F,F],     idx: 2, to: F },
    { desc: "F→D inside active run",           chain: [F,F,M,F,F],     idx: 2, to: F },
  ];

  const cases = [
    { desc: "M→D next to D neighbour (left)",         chain: [D,M,M,M,M],     idx: 1, to: D },
    { desc: "M→D next to D neighbour (right)",        chain: [M,M,M,M,D],     idx: 3, to: D },
    { desc: "M→D between two D neighbours",           chain: [D,M,D,M,M],     idx: 1, to: D },
    { desc: "D→M breaking D-D pair",                  chain: [D,D,M,M,M],     idx: 0, to: M },
    { desc: "D→M in middle of D block",               chain: [D,D,D,M,M],     idx: 1, to: M },
    { desc: "D→F next to D (loses D-D, gains nothing since run<minRun)", chain: [D,D,M,M,M], idx: 1, to: F },
    ...baseCases
  ];

  cases.forEach(({ desc, chain: c, idx, to }) => {
    it(desc, () => {
      const chain = [...c];
      const { dE, expected } = checkDelta(chain, idx, to, P);
      expect(dE).toBeCloseTo(expected);
    });
  });
});

// ── buildRunLen ────────────────────────────────────────────────────────────

describe("buildRunLen", () => {
  it("all monomers → all zeros", () => {
    const rl = buildRunLen([M, M, M, M]);
    expect(Array.from(rl)).toEqual([0, 0, 0, 0]);
  });

  it("single fibril → [1]", () => {
    const rl = buildRunLen([M, F, M]);
    expect(Array.from(rl)).toEqual([0, 1, 0]);
  });

  it("run of 3 → all cells report 3", () => {
    const rl = buildRunLen([M, F, F, F, M]);
    expect(Array.from(rl)).toEqual([0, 3, 3, 3, 0]);
  });

  it("two separate runs", () => {
    const rl = buildRunLen([F, F, M, F, F, F]);
    expect(Array.from(rl)).toEqual([2, 2, 0, 3, 3, 3]);
  });
});

// ── countActiveRuns ────────────────────────────────────────────────────────

describe("countActiveRuns", () => {
  it("no runs → 0", () => {
    const rl = buildRunLen([M, M, M]);
    expect(countActiveRuns(rl, 3)).toBe(0);
  });

  it("one run below threshold → 0", () => {
    const rl = buildRunLen([F, F, M, M]);
    expect(countActiveRuns(rl, 3)).toBe(0);
  });

  it("one active run → 1", () => {
    const rl = buildRunLen([F, F, F, M]);
    expect(countActiveRuns(rl, 3)).toBe(1);
  });

  it("two active runs → 2", () => {
    const rl = buildRunLen([F, F, F, M, F, F, F]);
    expect(countActiveRuns(rl, 3)).toBe(2);
  });

  it("one active + one sub-threshold → 1", () => {
    const rl = buildRunLen([F, F, F, M, F, F]);
    expect(countActiveRuns(rl, 3)).toBe(1);
  });
});

// ── updateRunLen ───────────────────────────────────────────────────────────

describe("updateRunLen — matches buildRunLen after every flip", () => {
  function checkUpdate(startChain, idx, newState) {
    const chain = [...startChain];
    const runLen = buildRunLen(chain);
    chain[idx] = newState;
    updateRunLen(chain, runLen, idx);
    const expected = Array.from(buildRunLen(chain));
    expect(Array.from(runLen)).toEqual(expected);
  }

  it("M→F merges two runs into one", () =>
    checkUpdate([F,F,F,M,F,F,F,F,F], 3, F));

  it("F→M splits one run into two", () =>
    checkUpdate([F,F,F,F,F,F,F], 3, M));

  it("M→F creates first run", () =>
    checkUpdate([M,M,M,M,M], 2, F));

  it("F→M removes only fibril site", () =>
    checkUpdate([M,F,M,M,M], 1, M));

  it("M→F extends run on left", () =>
    checkUpdate([F,F,F,M,M], 3, F));

  it("M→F extends run on right", () =>
    checkUpdate([M,M,F,F,F], 2, F));

  it("F→M shrinks run from left", () =>
    checkUpdate([F,F,F,F,F], 0, M));

  it("F→M shrinks run from right", () =>
    checkUpdate([F,F,F,F,F], 4, M));

  it("M→F merges a very long run", () => {
    const chain = new Array(20).fill(F);
    chain[10] = M;
    checkUpdate(chain, 10, F);
  });
});

// ── deltaEnergyFast — consistency with computeEnergy ──────────────────────

function checkDeltaFast(chain, idx, newState, params) {
  const oldState = chain[idx];
  const { nF, nActiveRuns, runLen } = buildSimState(chain, params.minRun);
  const oldE = computeEnergy(chain, params);
  chain[idx] = newState;
  updateRunLen(chain, runLen, idx);
  const dE = deltaEnergyFast(chain, idx, oldState, newState, params, nF, nActiveRuns, runLen);
  const newE = computeEnergy(chain, params);
  chain[idx] = oldState; // restore
  return { dE, expected: newE - oldE };
}

describe("deltaEnergyFast — exact consistency with computeEnergy", () => {
  const cases = [
    { desc: "M→F merging two active runs",       chain: [F,F,F,M,F,F,F], idx: 3, to: F },
    { desc: "F→M breaking an active run",        chain: [F,F,F,F,F],     idx: 2, to: M },
    { desc: "M→D with no interactions",          chain: [M,M,M,M,M],     idx: 2, to: D },
    { desc: "D→M inside D block",                chain: [D,D,D,M,F,F,F], idx: 1, to: M },
    { desc: "M→F creating first active run",     chain: [F,F,M,M,M],     idx: 2, to: F },
    { desc: "F→M removing only active run",      chain: [F,F,F,M,M],     idx: 1, to: M },
    { desc: "M→D next to D neighbour",           chain: [D,M,M,M,M],     idx: 1, to: D },
    { desc: "D→M breaking D-D pair",             chain: [D,D,M,M,M],     idx: 0, to: M },
    { desc: "F→D inside active run",             chain: [F,F,F,F,F],     idx: 2, to: D },
    { desc: "M→F extending run past minRun",     chain: [M,F,F,M,M],     idx: 0, to: F },
    { desc: "F→M with hFF active elsewhere",     chain: [F,F,F,M,F,F,F,M,F], idx: 7, to: M },
  ];

  cases.forEach(({ desc, chain: c, idx, to }) => {
    it(desc, () => {
      const chain = [...c];
      const { dE, expected } = checkDeltaFast(chain, idx, to, P);
      expect(dE).toBeCloseTo(expected);
    });
  });
});

// ── simulation helpers ─────────────────────────────────────────────────────

describe("initChain", () => {
  it("returns all-monomer chain of length n", () => {
    const c = initChain(10);
    expect(c).toHaveLength(10);
    expect(c.every(s => s === STATES.M)).toBe(true);
  });
});

describe("countStates", () => {
  it("counts correctly for a mixed chain", () => {
    const { M, D, F } = STATES;
    const chain = [M, M, D, F, F, F];
    const [nM, nD, nF] = countStates(chain);
    expect(nM).toBe(2);
    expect(nD).toBe(1);
    expect(nF).toBe(3);
  });

  it("sums to chain length", () => {
    const chain = initChain(20);
    const counts = countStates(chain);
    expect(counts.reduce((a, b) => a + b, 0)).toBe(20);
  });
});

describe("fibrilRunLengths", () => {
  const { M, D, F } = STATES;

  it("returns empty array for no fibrils", () => {
    expect(fibrilRunLengths([M, M, D, M])).toEqual([]);
  });

  it("returns correct run lengths for multiple blocks", () => {
    // F F M F F F M F → runs of 2, 3, 1
    const chain = [F, F, M, F, F, F, M, F];
    expect(fibrilRunLengths(chain)).toEqual([2, 3, 1]);
  });
});

describe("mcStep", () => {
  it("returns a chain of the same length", () => {
    const chain = initChain(20);
    const E = computeEnergy(chain, DEFAULT_PARAMS);
    const { chain: nc } = mcStep(chain, DEFAULT_PARAMS, 1.0, null, E);
    expect(nc).toHaveLength(20);
  });

  it("never modifies locked sites", () => {
    const { M, F } = STATES;
    // All sites locked as fibril
    const chain  = new Array(10).fill(F);
    const locked = new Array(10).fill(true);
    const E = computeEnergy(chain, DEFAULT_PARAMS);
    const { chain: nc } = mcStep(chain, DEFAULT_PARAMS, 10.0, locked, E);
    expect(nc.every(s => s === F)).toBe(true);
  });

  it("at T→0 never increases energy", () => {
    // Very low temperature — uphill moves should be rejected
    const chain = initChain(30);
    const T = 0.01;
    let c = chain;
    for (let i = 0; i < 50; i++) {
      const before = computeEnergy(c, DEFAULT_PARAMS);
      const { chain: nc } = mcStep(c, DEFAULT_PARAMS, T, null, computeEnergy(c, DEFAULT_PARAMS));
      const after = computeEnergy(nc, DEFAULT_PARAMS);
      expect(after).toBeLessThanOrEqual(before + 1e-9);
      c = nc;
    }
  });
});
