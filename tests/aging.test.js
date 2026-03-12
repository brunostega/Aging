import { describe, it, expect } from "vitest";
import { STATES, fibrilRunLength, computeEnergy, localEnergy } from "../src/model.js";
import { initChain, countStates, fibrilRunLengths, mcStep, DEFAULT_PARAMS } from "../src/simulation.js";

const { M, D, F } = STATES;

// ── shared parameter set ────────────────────────────────────────────────────
const P = { eM: 1, eD: 2, eF: 4, jF: 3, jFF: 1, jD: 2, minRun: 3 };

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
  it("fibril run < minRun: no coupling applied", () => {
    // run of 2, minRun=3 → only intrinsic energy
    const chain = [F, F];
    const expected = 2 * P.eF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("fibril run exactly minRun=3: 2 nn pairs + 1 long-range pair (0,2)", () => {
    // [F,F,F]: nn pairs=(0,1),(1,2)=2; long-range |d|>1: (0,2)=1
    const chain = [F, F, F];
    const expected = 3 * P.eF - 2 * P.jF - 1 * P.jFF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("fibril run of 4: 3 short-range pairs + 2 long-range pairs", () => {
    // [F,F,F,F]: nn pairs = (0,1),(1,2),(2,3) = 3
    // long-range |d|>1 pairs = (0,2),(0,3),(1,3) = 3
    const chain = [F, F, F, F];
    const expected = 4 * P.eF - 3 * P.jF - 3 * P.jFF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("fibril run of 5: 4 nn pairs + 6 long-range pairs", () => {
    // nn pairs: (0,1),(1,2),(2,3),(3,4) = 4
    // long-range pairs: (0,2),(0,3),(0,4),(1,3),(1,4),(2,4) = 6
    const chain = [F, F, F, F, F];
    const expected = 5 * P.eF - 4 * P.jF - 6 * P.jFF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("two separate active fibril runs do not interact via jFF", () => {
    // [F,F,F, M, F,F,F]: two runs of 3 separated by M
    // each run: 2 nn pairs, 0 long-range within run (run=3, only |d|>1 is the pair 0-2 within run)
    // across runs: sites are fibril but fibrilRunLength check happens per-site — both runs active
    // within run of 3: long-range pairs = (0,2) only = 1 per run
    // cross-run pairs: all |d|>1 between runs, both active → counted
    const chain = [F, F, F, M, F, F, F];
    // run1 sites: 0,1,2  run2 sites: 4,5,6
    // nn pairs within runs: (0,1),(1,2),(4,5),(5,6) = 4 × jF
    // long-range within run1: (0,2)=1; within run2: (4,6)=1
    // long-range cross-run: (0,4),(0,5),(0,6),(1,4),(1,5),(1,6),(2,4),(2,5),(2,6) = 9
    // total long-range = 1+1+9 = 11 × jFF
    const expected = 6 * P.eF + P.eM - 4 * P.jF - 11 * P.jFF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });
});

// ── computeEnergy: mixed states ────────────────────────────────────────────
describe("computeEnergy — mixed states", () => {
  it("[F,F,F,M,M] with minRun=3: only fibril block active", () => {
    // fibril block [0,1,2]: run=3=minRun, nn pairs=2, long-range (0,2)=1
    const chain = [F, F, F, M, M];
    const expected = 3 * P.eF + 2 * P.eM - 2 * P.jF - 1 * P.jFF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("[M,D,D,F,F,F] with minRun=3: D-D and F-F interactions both active", () => {
    // D sites: 1,2 → 1 D-D nn pair → −jD
    // F sites: 3,4,5 → run=3=minRun, nn pairs=2 → −2×jF, long-range (3,5)=1 → −jFF
    const chain = [M, D, D, F, F, F];
    const expected = P.eM + 2*P.eD + 3*P.eF - P.jD - 2*P.jF - P.jFF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });
});

// ── localEnergy consistency ────────────────────────────────────────────────
// localEnergy counts the full interaction weight for site idx (not halved),
// so it can be used directly as ΔE in MC without double-counting —
// what matters is that Δ(localEnergy) == Δ(computeEnergy) on a single-site flip.

describe("localEnergy — consistency with computeEnergy", () => {
  it("localEnergy of isolated M site equals eM", () => {
    expect(localEnergy([M], 0, P)).toBeCloseTo(P.eM);
  });

  it("localEnergy of isolated D site equals eD", () => {
    expect(localEnergy([D], 0, P)).toBeCloseTo(P.eD);
  });

  it("localEnergy of isolated F site (run=1 < minRun) equals eF", () => {
    expect(localEnergy([M, F, M], 1, P)).toBeCloseTo(P.eF);
  });

  it("ΔE from localEnergy matches global ΔE on a single-site flip (M→D)", () => {
    const chain = [F, F, F, M, D, D];
    const idx = 3; // M site, no neighbours that interact with it
    const oldE = localEnergy(chain, idx, P);
    const flipped = [...chain]; flipped[idx] = D;
    const newE = localEnergy(flipped, idx, P);
    expect(newE - oldE).toBeCloseTo(computeEnergy(flipped, P) - computeEnergy(chain, P));
  });

  it("ΔE from localEnergy matches global ΔE on a single-site flip (M→F inside active run)", () => {
    // Flipping the middle M to F extends a run and activates new couplings
    const chain = [F, F, M, F, F, F];
    const idx = 2;
    const oldE = localEnergy(chain, idx, P);
    const flipped = [...chain]; flipped[idx] = F;
    const newE = localEnergy(flipped, idx, P);
    expect(newE - oldE).toBeCloseTo(computeEnergy(flipped, P) - computeEnergy(chain, P));
  });

  it("ΔE from localEnergy matches global ΔE on a single-site flip (D→M inside D block)", () => {
    const chain = [D, D, D, M, F, F, F];
    const idx = 1;
    const oldE = localEnergy(chain, idx, P);
    const flipped = [...chain]; flipped[idx] = M;
    const newE = localEnergy(flipped, idx, P);
    expect(newE - oldE).toBeCloseTo(computeEnergy(flipped, P) - computeEnergy(chain, P));
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
    const { chain: nc } = mcStep(chain, DEFAULT_PARAMS, 1.0, null);
    expect(nc).toHaveLength(20);
  });

  it("never modifies locked sites", () => {
    const { M, F } = STATES;
    // All sites locked as fibril
    const chain  = new Array(10).fill(F);
    const locked = new Array(10).fill(true);
    const { chain: nc } = mcStep(chain, DEFAULT_PARAMS, 10.0, locked);
    expect(nc.every(s => s === F)).toBe(true);
  });

  it("at T→0 never increases energy", () => {
    // Very low temperature — uphill moves should be rejected
    const chain = initChain(30);
    const T = 0.01;
    let c = chain;
    for (let i = 0; i < 50; i++) {
      const before = computeEnergy(c, DEFAULT_PARAMS);
      const { chain: nc } = mcStep(c, DEFAULT_PARAMS, T, null);
      const after = computeEnergy(nc, DEFAULT_PARAMS);
      expect(after).toBeLessThanOrEqual(before + 1e-9);
      c = nc;
    }
  });
});
