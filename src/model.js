// Pure energy model — no React dependencies.
// Imported by simulation.js, App.jsx, and the test suite.
//
// Hamiltonian:
//   H = Σᵢ ε[sᵢ]
//     − J_D Σ⟨i,j⟩ δ(D,D)
//     − J_F Σ⟨i,j⟩ δ(F,F) · 𝟙[runLen ≥ minRun]
//     − h_FF · N_F · 𝟙[∃ active run]

export const STATES = { M: 0, D: 1, F: 2 };

// ── Run-length cache ─────────────────────────────────────────────────────────
//
// runLen[i] = length of the contiguous F-run containing site i, or 0 if i is
// not F. Built in O(N); updated incrementally in O(runLength) after each flip.

// O(runLength) — used for one-off queries outside the MC loop.
export function fibrilRunLength(chain, i) {
  if (chain[i] !== STATES.F) return 0;
  const N = chain.length;
  let lo = i, hi = i;
  while (lo > 0     && chain[lo - 1] === STATES.F) lo--;
  while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;
  return hi - lo + 1;
}

// O(N) — builds the full runLen cache from scratch.
export function buildRunLen(chain) {
  const N = chain.length;
  const runLen = new Int32Array(N);
  let i = 0;
  while (i < N) {
    if (chain[i] !== STATES.F) { i++; continue; }
    let j = i;
    while (j < N && chain[j] === STATES.F) j++;
    const len = j - i;
    for (let k = i; k < j; k++) runLen[k] = len;
    i = j;
  }
  return runLen;
}

// O(runLength) — updates runLen in-place after flipping chain[idx].
// Must be called AFTER the flip has been applied to chain.
export function updateRunLen(chain, runLen, idx) {
  const N = chain.length;
  let lo = idx, hi = idx;
  while (lo > 0     && chain[lo - 1] === STATES.F) lo--;
  while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;

  if (chain[idx] !== STATES.F) {
    // idx left a run — split into left [lo..idx-1] and right [idx+1..hi]
    runLen[idx] = 0;
    if (lo < idx) { const len = idx - lo;     for (let k = lo;      k < idx; k++) runLen[k] = len; }
    if (idx < hi) { const len = hi  - idx;    for (let k = idx + 1; k <= hi; k++) runLen[k] = len; }
  } else {
    // idx joined a run — merge everything in [lo..hi]
    const len = hi - lo + 1;
    for (let k = lo; k <= hi; k++) runLen[k] = len;
  }
}

// O(N/avgRunLen) — counts runs whose length ≥ minRun by jumping run-by-run.
export function countActiveRuns(runLen, minRun) {
  const N = runLen.length;
  let count = 0, i = 0;
  while (i < N) {
    if (runLen[i] === 0) { i++; continue; }
    const len = runLen[i];
    if (len >= minRun) count++;
    i += len;  // jump to next run
  }
  return count;
}

// O(N) — builds runLen, nF, and nActiveRuns in a single pass.
// Returns everything mcStep needs to maintain incremental state.
export function buildSimState(chain, minRun) {
  const N = chain.length;
  const runLen = new Int32Array(N);
  let nF = 0, i = 0;
  while (i < N) {
    if (chain[i] !== STATES.F) { i++; continue; }
    let j = i;
    while (j < N && chain[j] === STATES.F) { nF++; j++; }
    const len = j - i;
    for (let k = i; k < j; k++) runLen[k] = len;
    i = j;
  }
  const nActiveRuns = countActiveRuns(runLen, minRun);
  return { nF, nActiveRuns, runLen };
}

// ── Energy ───────────────────────────────────────────────────────────────────

// O(N) — computes full energy from scratch. Used for initialisation and tests;
// the MC loop uses deltaEnergyFast instead.
export function computeEnergy(chain, params) {
  const { eM, eD, eF, jF, jD, hFF, minRun } = params;
  const N = chain.length;
  const baseE = [eM, eD, eF];
  const runLen = buildRunLen(chain);
  let E = 0, nF = 0, hasActive = false;

  for (let i = 0; i < N; i++) {
    E += baseE[chain[i]];
    if (chain[i] === STATES.F) {
      nF++;
      if (runLen[i] >= minRun) hasActive = true;
    } else if (chain[i] === STATES.D) {
      if (i > 0     && chain[i - 1] === STATES.D) E -= jD / 2;
      if (i < N - 1 && chain[i + 1] === STATES.D) E -= jD / 2;
    }
  }

  // Second pass for jF: requires knowing which sites are in active runs
  for (let i = 0; i < N; i++) {
    if (chain[i] === STATES.F && runLen[i] >= minRun) {
      if (i > 0     && chain[i - 1] === STATES.F && runLen[i - 1] >= minRun) E -= jF / 2;
      if (i < N - 1 && chain[i + 1] === STATES.F && runLen[i + 1] >= minRun) E -= jF / 2;
    }
  }

  if (hasActive && nF > 0) E -= hFF * nF;
  return E;
}

// ── Local ΔE ─────────────────────────────────────────────────────────────────
//
// When chain[idx] flips from oldState to newState, only these terms change:
//   1. Intrinsic energy of idx
//   2. D-D coupling of idx with its neighbours
//   3. jF coupling within the affected fibril window
//   4. h_FF term (depends on nF and whether any active run exists)
//
// Strategy: compute patch energy of the affected window before and after,
// then compute the h_FF delta from the before/after nActiveRuns counters.

// Returns [lo, hi]: the window of sites whose patchEnergy may change.
// Expands to cover any adjacent fibril run (jF is short-range within runs).
// D-D coupling is handled via boundary terms, so D neighbours stay outside.
function affectedWindow(chain, idx) {
  const N = chain.length;
  let lo = idx, hi = idx;
  if (idx > 0     && chain[idx - 1] === STATES.F) { lo = idx - 1; while (lo > 0     && chain[lo - 1] === STATES.F) lo--; }
  if (idx < N - 1 && chain[idx + 1] === STATES.F) { hi = idx + 1; while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++; }
  return [lo, hi];
}

// Energy of the patch [lo..hi] plus its boundary couplings with [lo-1] and [hi+1].
function patchEnergy(chain, lo, hi, baseE, jD, jF, minRun, runLen) {
  const N = chain.length;
  let E = 0;
  for (let i = lo; i <= hi; i++) {
    E += baseE[chain[i]];
    if (i < hi) {
      if (chain[i] === STATES.D && chain[i + 1] === STATES.D) E -= jD;
      if (chain[i] === STATES.F && chain[i + 1] === STATES.F &&
          runLen[i] >= minRun && runLen[i + 1] >= minRun) E -= jF;
    }
  }
  // Left boundary coupling
  if (lo > 0) {
    if (chain[lo] === STATES.D && chain[lo - 1] === STATES.D) E -= jD;
    if (chain[lo] === STATES.F && chain[lo - 1] === STATES.F &&
        runLen[lo] >= minRun && runLen[lo - 1] >= minRun) E -= jF;
  }
  // Right boundary coupling
  if (hi < N - 1) {
    if (chain[hi] === STATES.D && chain[hi + 1] === STATES.D) E -= jD;
    if (chain[hi] === STATES.F && chain[hi + 1] === STATES.F &&
        runLen[hi] >= minRun && runLen[hi + 1] >= minRun) E -= jF;
  }
  return E;
}

// Exact ΔE for flipping chain[idx] from oldState → newState.
// Preconditions:
//   - chain[idx] already holds newState
//   - runLen already reflects newState
//   - nF and nActiveRuns are the BEFORE-flip values
export function deltaEnergyFast(chain, idx, oldState, newState, params, nF, nActiveRuns, runLen) {
  const { eM, eD, eF, jD, jF, hFF, minRun } = params;
  const baseE = [eM, eD, eF];
  const N = chain.length;

  // ── Fast path: M↔D flip — no fibril run affected, no hFF change ──────────
  if (oldState !== STATES.F && newState !== STATES.F) {
    let dE = baseE[newState] - baseE[oldState];
    // D-D coupling with left neighbour
    if (idx > 0 && chain[idx - 1] === STATES.D) {
      if (newState === STATES.D) dE -= jD;
      if (oldState === STATES.D) dE += jD;
    }
    // D-D coupling with right neighbour
    if (idx < N - 1 && chain[idx + 1] === STATES.D) {
      if (newState === STATES.D) dE -= jD;
      if (oldState === STATES.D) dE += jD;
    }
    return dE;
  }

  // ── General path: at least one state is F ────────────────────────────────
  const [lo, hi] = affectedWindow(chain, idx);
  const newPatch = patchEnergy(chain, lo, hi, baseE, jD, jF, minRun, runLen);

  // Temporarily revert to compute old patch
  chain[idx] = oldState;
  updateRunLen(chain, runLen, idx);
  const oldPatch = patchEnergy(chain, lo, hi, baseE, jD, jF, minRun, runLen);

  // Restore new state
  chain[idx] = newState;
  updateRunLen(chain, runLen, idx);

  const newNF          = nF + (newState === STATES.F ? 1 : 0) - (oldState === STATES.F ? 1 : 0);
  const newNActiveRuns = countActiveRuns(runLen, minRun);
  const oldHFF         = nActiveRuns    > 0 ? hFF * nF    : 0;
  const newHFF         = newNActiveRuns > 0 ? hFF * newNF : 0;

  return (newPatch - oldPatch) - (newHFF - oldHFF);
}
