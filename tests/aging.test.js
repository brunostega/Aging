import { describe, it, expect } from "vitest";
import { STATES, fibrilRunLength, computeEnergy, localEnergy } from "../src/model.js";

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

  it("fibril run exactly minRun=3: short-range only, no long-range", () => {
    // [F,F,F]: 2 nn pairs, 0 long-range pairs (all |d|≤1 or same site)
    const chain = [F, F, F];
    const expected = 3 * P.eF - 2 * P.jF;
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
    // fibril block [0,1,2]: run=3=minRun, nn pairs=2, long-range=(0,2)=1
    // no D-D interactions
    const chain = [F, F, F, M, M];
    const expected = 3 * P.eF + 2 * P.eM - 2 * P.jF - 1 * P.jFF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });

  it("[M,D,D,F,F,F] with minRun=3: D-D and F-F interactions both active", () => {
    // D sites: 1,2 → 1 D-D pair → −jD
    // F sites: 3,4,5 → run=3, nn pairs=2 → −2×jF, long-range (3,5)=1 → −jFF
    const chain = [M, D, D, F, F, F];
    const expected = P.eM + 2*P.eD + 3*P.eF - P.jD - 2*P.jF - P.jFF;
    expect(computeEnergy(chain, P)).toBeCloseTo(expected);
  });
});

// ── localEnergy consistency ────────────────────────────────────────────────
describe("localEnergy — consistency with computeEnergy", () => {
  it("sum of localEnergy over all sites equals 2 × computeEnergy", () => {
    // Each pairwise interaction is counted once in computeEnergy (with /2 symmetry factor)
    // but twice in sum of localEnergy (once per site). So sum = 2 × E.
    const chain = [F, F, F, D, D, M, F, F, F];
    const total = chain.reduce((acc, _, i) => acc + localEnergy(chain, i, P), 0);
    expect(total).toBeCloseTo(2 * computeEnergy(chain, P));
  });

  it("localEnergy of isolated M site equals eM", () => {
    expect(localEnergy([M], 0, P)).toBeCloseTo(P.eM);
  });

  it("localEnergy of isolated D site equals eD", () => {
    expect(localEnergy([D], 0, P)).toBeCloseTo(P.eD);
  });

  it("localEnergy of isolated F site (run=1 < minRun) equals eF", () => {
    expect(localEnergy([M, F, M], 1, P)).toBeCloseTo(P.eF);
  });

  it("ΔE from localEnergy matches actual energy difference on flip", () => {
    const chain = [F, F, F, M, D, D];
    const idx = 3; // M site
    const oldE = localEnergy(chain, idx, P);
    const flipped = [...chain];
    flipped[idx] = D;
    const newE = localEnergy(flipped, idx, P);
    const dE = newE - oldE;
    const globalDelta = computeEnergy(flipped, P) - computeEnergy(chain, P);
    expect(dE).toBeCloseTo(globalDelta);
  });
});
