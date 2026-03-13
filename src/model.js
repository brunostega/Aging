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
