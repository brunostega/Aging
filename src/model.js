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

// Returns true if the chain contains at least one active fibril run (length >= minRun)
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
  let E = 0;
  let nF = 0;
  let hasActive = false;

  // Single pass: intrinsic energies, D-D and F-F nn couplings, count F sites and active runs
  let run = 0;
  for (let i = 0; i < N; i++) {
    E += baseE[chain[i]];

    if (chain[i] === STATES.F) {
      nF++;
      run++;
      if (run >= minRun) hasActive = true;
      // Short-range F-F nn coupling (active runs only, checked below)
    } else {
      run = 0;
    }

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

  // Background field: acts on all F sites when at least one active run exists
  if (hasActive && nF > 0) E -= hFF * nF;

  return E;
}

// ── Local ΔE for a single-site flip ────────────────────────────────────────
//
// When site idx changes state, only these contributions to E can change:
//   1. Intrinsic energy of idx
//   2. D-D coupling of idx with its two neighbours
//   3. jF coupling: the run containing idx, plus the immediately adjacent runs
//      (their length — and therefore active status — may change)
//   4. h_FF term: depends on nF (count of F sites) and hasActiveRun
//
// Strategy: evaluate the patch energy of the affected window before and after
// the flip. ΔE = patchE_new − patchE_old exactly.

// Sum of all energy contributions from sites in [lo..hi], counting each
// pairwise interaction once (left-to-right), plus the boundary interaction
// between lo-1 and lo.
function patchEnergy(chain, lo, hi, params) {
  const { eM, eD, eF, jF, jD, minRun } = params;
  const baseE = [eM, eD, eF];
  const N = chain.length;
  let E = 0;
  for (let i = lo; i <= hi; i++) {
    E += baseE[chain[i]];
    // Right-neighbour interactions (within window)
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

// Window that must be re-evaluated when idx is flipped.
// Expands to cover:
//   - any fibril run immediately adjacent (run length and active status may change)
//   - any immediate D neighbour (D-D coupling with idx may appear or disappear)
// Note: we only need one step for D since D-D is nearest-neighbour only.
function affectedWindow(chain, idx) {
  const N = chain.length;
  let lo = idx, hi = idx;

  // Expand left
  if (idx > 0) {
    if (chain[idx - 1] === STATES.F) {
      lo = idx - 1;
      while (lo > 0 && chain[lo - 1] === STATES.F) lo--;
    } else if (chain[idx - 1] === STATES.D) {
      lo = idx - 1;
    }
  }

  // Expand right
  if (idx < N - 1) {
    if (chain[idx + 1] === STATES.F) {
      hi = idx + 1;
      while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;
    } else if (chain[idx + 1] === STATES.D) {
      hi = idx + 1;
    }
  }

  return [lo, hi];
}

// Exact ΔE for flipping chain[idx] from oldState to newState.
// chain must contain newState at idx when called.
// nF and hadActive are the BEFORE-flip values (passed in to avoid recount).
export function deltaEnergy(chain, idx, oldState, newState, params, nF, hadActive) {
  const { hFF, minRun } = params;

  // Compute affected window using the new chain state
  const [lo, hi] = affectedWindow(chain, idx);

  // Patch energy with newState (chain already has newState at idx)
  const newPatch = patchEnergy(chain, lo, hi, params);

  // Temporarily revert to compute old patch
  chain[idx] = oldState;
  const oldPatch = patchEnergy(chain, lo, hi, params);
  chain[idx] = newState; // restore

  // h_FF delta: field acts on all F sites when any active run exists
  const newNF = nF + (newState === STATES.F ? 1 : 0) - (oldState === STATES.F ? 1 : 0);
  const newHasActive = hasActiveRun(chain, minRun); // chain has newState at idx
  const oldHFF = hadActive ? hFF * nF : 0;
  const newHFF = newHasActive ? hFF * newNF : 0;

  return (newPatch - oldPatch) + (newHFF - oldHFF) * -1;
  //                              ^ sign: h_FF term in H is −hFF*nF
}

