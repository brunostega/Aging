// Pure model functions — no React dependencies.
// Imported by both the UI and the test suite.

export const STATES = { M: 0, D: 1, F: 2 };

// Length of the contiguous fibril run containing site i
export function fibrilRunLength(chain, i) {
  if (chain[i] !== STATES.F) return 0;
  const N = chain.length;
  let lo = i, hi = i;
  while (lo > 0 && chain[lo - 1] === STATES.F) lo--;
  while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;
  return hi - lo + 1;
}

// Returns true if the chain contains at least one active fibril run
function hasActiveRun(chain, minRun) {
  const N = chain.length;
  let run = 0;
  for (let i = 0; i < N; i++) {
    run = chain[i] === STATES.F ? run + 1 : 0;
    if (run >= minRun) return true;
  }
  return false;
}

export function computeEnergy(chain, params) {
  const { eM, eD, eF, jF, jD, hFF, minRun } = params;
  const N = chain.length;
  const baseE = [eM, eD, eF];
  let E = 0, nF = 0, hasActive = false, run = 0;

  for (let i = 0; i < N; i++) {
    E += baseE[chain[i]];
    if (chain[i] === STATES.F) { nF++; run++; if (run >= minRun) hasActive = true; }
    else run = 0;
    if (chain[i] === STATES.D) {
      if (i > 0 && chain[i-1] === STATES.D) E -= jD / 2;
      if (i < N-1 && chain[i+1] === STATES.D) E -= jD / 2;
    }
  }

  // Second pass for jF: only now we know which sites are in active runs
  for (let i = 0; i < N; i++) {
    if (chain[i] === STATES.F && fibrilRunLength(chain, i) >= minRun) {
      if (i > 0     && chain[i-1] === STATES.F && fibrilRunLength(chain, i-1) >= minRun) E -= jF / 2;
      if (i < N - 1 && chain[i+1] === STATES.F && fibrilRunLength(chain, i+1) >= minRun) E -= jF / 2;
    }
  }

  if (hasActive && nF > 0) E -= hFF * nF;
  return E;
}

// ── Local ΔE for a single-site flip ────────────────────────────────────────
//
// When site idx changes state, only these contributions to E can change:
//   1. Intrinsic energy of idx
//   2. D-D coupling of idx with its two neighbours
//   3. jF coupling: the run containing idx, plus the immediately adjacent runs
//   4. h_FF term: depends on nF and hasActiveRun
//
// Strategy: evaluate the patch energy of the affected window before and after
// the flip. ΔE = patchE_new − patchE_old exactly.

function patchEnergy(chain, lo, hi, params) {
  const { eM, eD, eF, jF, jD, minRun } = params;
  const baseE = [eM, eD, eF];
  const N = chain.length;
  let E = 0;
  for (let i = lo; i <= hi; i++) {
    E += baseE[chain[i]];
    if (i < hi) {
      if (chain[i] === STATES.D && chain[i + 1] === STATES.D) E -= jD;
      if (chain[i] === STATES.F && chain[i + 1] === STATES.F &&
          fibrilRunLength(chain, i) >= minRun &&
          fibrilRunLength(chain, i + 1) >= minRun) E -= jF;
    }
  }
  // Left-boundary interaction between lo-1 and lo (outside window, unchanged)
  if (lo > 0) {
    if (chain[lo] === STATES.D && chain[lo - 1] === STATES.D) E -= jD;
    if (chain[lo] === STATES.F && chain[lo - 1] === STATES.F &&
        fibrilRunLength(chain, lo) >= minRun &&
        fibrilRunLength(chain, lo - 1) >= minRun) E -= jF;
  }
  // Right-boundary: interaction between hi and hi+1
  if (hi < N - 1) {
    if (chain[hi] === STATES.D && chain[hi + 1] === STATES.D) E -= jD;
    if (chain[hi] === STATES.F && chain[hi + 1] === STATES.F &&
        fibrilRunLength(chain, hi) >= minRun &&
        fibrilRunLength(chain, hi + 1) >= minRun) E -= jF;
  }
  return E;
}

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

// Exact ΔE for flipping chain[idx] from oldState to newState.
// chain must contain newState at idx when called.
// nF and hadActive are the BEFORE-flip values (passed in to avoid recount).
export function deltaEnergy(chain, idx, oldState, newState, params, nF, hadActive) {
  const { hFF, minRun } = params;
  const [lo, hi] = affectedWindow(chain, idx);

  const newPatch = patchEnergy(chain, lo, hi, params);

  chain[idx] = oldState;
  const oldPatch = patchEnergy(chain, lo, hi, params);
  chain[idx] = newState; // restore

  const newNF = nF + (newState === STATES.F ? 1 : 0) - (oldState === STATES.F ? 1 : 0);
  const newHasActive = hasActiveRun(chain, minRun); // chain has newState at idx
  const oldHFF = hadActive ? hFF * nF : 0;
  const newHFF = newHasActive ? hFF * newNF : 0;

  return (newPatch - oldPatch) - (newHFF - oldHFF);
}

// ── O(1) run-length cache ───────────────────────────────────────────────────

// Build runLen[i] in one O(N) pass.
// runLen[i] = length of the fibril run containing site i, or 0 if not F.
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

// Count active runs (runs with length >= minRun) from a runLen array — O(N).
export function countActiveRuns(runLen, minRun) {
  const N = runLen.length;
  let count = 0, i = 0;
  while (i < N) {
    if (runLen[i] === 0) { i++; continue; }
    const len = runLen[i];
    if (len >= minRun) count++;
    i += len; // jump to end of run
  }
  return count;
}

// Build all sim state needed by mcStep in one O(N) pass.
// Returns { nF, nActiveRuns, runLen }
export function buildSimState(chain, minRun) {
  const runLen = buildRunLen(chain);
  let nF = 0;
  for (let i = 0; i < chain.length; i++) if (chain[i] === STATES.F) nF++;
  const nActiveRuns = countActiveRuns(runLen, minRun);
  return { nF, nActiveRuns, runLen };
}

// Update runLen in-place after flipping site idx.
// Only recomputes runs that overlap the affected window.
// Returns the new run length at idx (useful for incremental nActiveRuns update).
export function updateRunLen(chain, runLen, idx) {
  const N = chain.length;
  // Find the full extent of the region that needs recomputing:
  // walk left/right from idx as far as F sites go (in the new chain state).
  let lo = idx, hi = idx;
  while (lo > 0 && chain[lo - 1] === STATES.F) lo--;
  while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;

  if (chain[idx] !== STATES.F) {
    // idx is no longer F — zero it, then recompute the two flanking runs
    runLen[idx] = 0;
    // Left flank: [lo .. idx-1]
    if (lo < idx) {
      const len = idx - lo;
      for (let k = lo; k < idx; k++) runLen[k] = len;
    }
    // Right flank: [idx+1 .. hi]
    if (idx < hi) {
      const len = hi - idx;
      for (let k = idx + 1; k <= hi; k++) runLen[k] = len;
    }
  } else {
    // idx is now F — the run spans [lo .. hi]
    const len = hi - lo + 1;
    for (let k = lo; k <= hi; k++) runLen[k] = len;
  }
}

// Incrementally update nActiveRuns after a flip.
// Called AFTER chain[idx] and runLen have been updated to the new state.
// oldRunLen: the runLen[idx] value BEFORE the flip (pass 0 if oldState !== F).
// Returns the updated nActiveRuns.
export function updateActiveRuns(runLen, nActiveRuns, idx, oldRunLen, minRun) {
  const newRunLen = runLen[idx];
  // Did the old run cross the threshold?
  const oldWasActive = oldRunLen >= minRun;
  const newIsActive  = newRunLen >= minRun;
  if (oldWasActive === newIsActive) return nActiveRuns; // no change
  // The run at idx changed active status — but we need to know if it's a new
  // run or the same run that was modified. Since runLen reflects the new merged/
  // split state, we check the flanking sites to count distinct runs affected.
  // Simpler: recount from runLen — O(N/runLen) amortised but exact.
  // For correctness we just recount — it's still O(N) worst case but only
  // called when active status changes (rare once mostly fibril).
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

// patchEnergy variant that uses a pre-built runLen for O(1) lookups.
export function patchEnergyFast(chain, lo, hi, params, runLen) {
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

// deltaEnergy variant using the runLen cache.
// chain must contain newState at idx. runLen must reflect newState.
// nActiveRuns: number of active runs BEFORE the flip.
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
  // nActiveRuns after: recount from updated runLen
  let newNActiveRuns = 0;
  { let i = 0; const N = runLen.length;
    while (i < N) { if (runLen[i] === 0) { i++; continue; }
      const len = runLen[i]; if (len >= minRun) newNActiveRuns++; i += len; } }

  const oldHFF = nActiveRuns > 0 ? hFF * nF : 0;
  const newHFF = newNActiveRuns > 0 ? hFF * newNF : 0;

  return (newPatch - oldPatch) - (newHFF - oldHFF);
}
