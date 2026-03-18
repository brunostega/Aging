// Pure model functions — no React dependencies.
// Imported by both the UI and the test suite.

export const STATES = { M: 0, D: 1, F: 2 };

// Length of the contiguous fibril run containing site i — O(runLength).
// Used for one-off calls outside the MC loop. mcStep uses the runLen cache.
export function fibrilRunLength(chain, i) {
  if (chain[i] !== STATES.F) return 0;
  const N = chain.length;
  let lo = i, hi = i;
  while (lo > 0 && chain[lo - 1] === STATES.F) lo--;
  while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;
  return hi - lo + 1;
}

// ── Run-length cache ─────────────────────────────────────────────────────────

// Build runLen[i] in one O(N) pass.
// runLen[i] = length of the fibril run at site i, or 0 if site is not F.
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

// Update runLen in-place after flipping site idx — O(runLength).
export function updateRunLen(chain, runLen, idx) {
  const N = chain.length;
  let lo = idx, hi = idx;
  while (lo > 0 && chain[lo - 1] === STATES.F) lo--;
  while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;

  if (chain[idx] !== STATES.F) {
    runLen[idx] = 0;
    // Left flank
    if (lo < idx) { const len = idx - lo; for (let k = lo; k < idx; k++) runLen[k] = len; }
    // Right flank
    if (idx < hi) { const len = hi - idx; for (let k = idx + 1; k <= hi; k++) runLen[k] = len; }
  } else {
    // Merged run
    const len = hi - lo + 1;
    for (let k = lo; k <= hi; k++) runLen[k] = len;
  }
}

// Count active runs (length >= minRun) from a runLen array — O(N/avgRunLen).
export function countActiveRuns(runLen, minRun) {
  let count = 0, i = 0;
  const N = runLen.length;
  while (i < N) {
    if (runLen[i] === 0) { i++; continue; }
    const len = runLen[i];
    if (len >= minRun) count++;
    i += len;
  }
  return count;
}

// Build all incremental state needed by mcStep in one O(N) pass.
// Returns { nF, nActiveRuns, runLen }
export function buildSimState(chain, minRun) {
  const N = chain.length;
  const runLen = new Int32Array(N);
  let nF = 0;
  let i = 0;
  while (i < N) {
    if (chain[i] !== STATES.F) { i++; continue; }
    nF++;
    let j = i + 1;
    while (j < N && chain[j] === STATES.F) { nF++; j++; }
    const len = j - i;
    for (let k = i; k < j; k++) runLen[k] = len;
    i = j;
  }
  const nActiveRuns = countActiveRuns(runLen, minRun);
  return { nF, nActiveRuns, runLen };
}

// ── Energy ───────────────────────────────────────────────────────────────────

export function computeEnergy(chain, params) {
  const { eM, eD, eF, jF, jD, hFF, minRun } = params;
  const N = chain.length;
  const baseE = [eM, eD, eF];
  let E = 0, nF = 0, hasActive = false;

  const runLen = buildRunLen(chain);

  for (let i = 0; i < N; i++) {
    E += baseE[chain[i]];
    if (chain[i] === STATES.F) {
      nF++;
      if (runLen[i] >= minRun) hasActive = true;
    }
    if (chain[i] === STATES.D) {
      if (i > 0     && chain[i-1] === STATES.D) E -= jD / 2;
      if (i < N - 1 && chain[i+1] === STATES.D) E -= jD / 2;
    }
  }

  for (let i = 0; i < N; i++) {
    if (chain[i] === STATES.F && runLen[i] >= minRun) {
      if (i > 0     && chain[i-1] === STATES.F && runLen[i-1] >= minRun) E -= jF / 2;
      if (i < N - 1 && chain[i+1] === STATES.F && runLen[i+1] >= minRun) E -= jF / 2;
    }
  }

  if (hasActive && nF > 0) E -= hFF * nF;
  return E;
}

// ── Local ΔE ─────────────────────────────────────────────────────────────────
//
// When site idx changes state, only these contributions to E can change:
//   1. Intrinsic energy of idx
//   2. D-D coupling of idx with its two neighbours
//   3. jF coupling: the run containing idx plus immediately adjacent runs
//   4. h_FF term: depends on nF and nActiveRuns
//
// Strategy: evaluate patch energy of the affected window before and after.
// ΔE = patchE_new − patchE_old exactly.

function affectedWindow(chain, idx) {
  const N = chain.length;
  let lo = idx, hi = idx;
  if (idx > 0) {
    if      (chain[idx-1] === STATES.F) { lo = idx-1; while (lo > 0 && chain[lo-1] === STATES.F) lo--; }
    else if (chain[idx-1] === STATES.D) { lo = idx-1; }
  }
  if (idx < N-1) {
    if      (chain[idx+1] === STATES.F) { hi = idx+1; while (hi < N-1 && chain[hi+1] === STATES.F) hi++; }
    else if (chain[idx+1] === STATES.D) { hi = idx+1; }
  }
  return [lo, hi];
}

function patchEnergyFast(chain, lo, hi, params, runLen) {
  const { eM, eD, eF, jF, jD, minRun } = params;
  const baseE = [eM, eD, eF];
  const N = chain.length;
  let E = 0;
  for (let i = lo; i <= hi; i++) {
    E += baseE[chain[i]];
    if (i < hi) {
      if (chain[i] === STATES.D && chain[i+1] === STATES.D) E -= jD;
      if (chain[i] === STATES.F && chain[i+1] === STATES.F &&
          runLen[i] >= minRun && runLen[i+1] >= minRun) E -= jF;
    }
  }
  if (lo > 0) {
    if (chain[lo] === STATES.D && chain[lo-1] === STATES.D) E -= jD;
    if (chain[lo] === STATES.F && chain[lo-1] === STATES.F &&
        runLen[lo] >= minRun && runLen[lo-1] >= minRun) E -= jF;
  }
  if (hi < N - 1) {
    if (chain[hi] === STATES.D && chain[hi+1] === STATES.D) E -= jD;
    if (chain[hi] === STATES.F && chain[hi+1] === STATES.F &&
        runLen[hi] >= minRun && runLen[hi+1] >= minRun) E -= jF;
  }
  return E;
}

// Exact ΔE for flipping chain[idx] from oldState to newState.
// chain must contain newState at idx; runLen must reflect newState.
// nActiveRuns is the BEFORE-flip count.
export function deltaEnergyFast(chain, idx, oldState, newState, params, nF, nActiveRuns, runLen) {
  const { hFF, minRun } = params;
  const [lo, hi] = affectedWindow(chain, idx);

  const newPatch = patchEnergyFast(chain, lo, hi, params, runLen);

  // Temporarily revert
  chain[idx] = oldState;
  updateRunLen(chain, runLen, idx);
  const oldPatch = patchEnergyFast(chain, lo, hi, params, runLen);

  // Restore
  chain[idx] = newState;
  updateRunLen(chain, runLen, idx);

  const newNF = nF + (newState === STATES.F ? 1 : 0) - (oldState === STATES.F ? 1 : 0);
  const newNActiveRuns = countActiveRuns(runLen, minRun);

  const oldHFF = nActiveRuns    > 0 ? hFF * nF    : 0;
  const newHFF = newNActiveRuns > 0 ? hFF * newNF : 0;

  return (newPatch - oldPatch) - (newHFF - oldHFF);
}
