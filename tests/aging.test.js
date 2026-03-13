import { describe, it, expect } from "vitest";
import { STATES, fibrilRunLength, computeEnergy } from "../src/model.js";
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
