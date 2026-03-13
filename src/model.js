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

// Start index of the contiguous fibril run containing site i.
// Used to determine whether two fibril sites belong to the same run.
function fibrilRunStart(chain, i) {
  let lo = i;
  while (lo > 0 && chain[lo - 1] === STATES.F) lo--;
  return lo;
}

export function computeEnergy(chain, params) {
  const { eM, eD, eF, jF, jFF, jD, minRun } = params;
  const N = chain.length;
  const baseE = [eM, eD, eF];
  let E = 0;
  for (let i = 0; i < N; i++) {
    E += baseE[chain[i]];
    if (chain[i] === STATES.F && fibrilRunLength(chain, i) >= minRun) {
      // Short-range: nearest neighbours in the same active run
      if (i > 0     && chain[i-1] === STATES.F && fibrilRunLength(chain, i-1) >= minRun) E -= jF / 2;
      if (i < N - 1 && chain[i+1] === STATES.F && fibrilRunLength(chain, i+1) >= minRun) E -= jF / 2;
      // Long-range: any active fibril site belonging to a DIFFERENT run
      const runI = fibrilRunStart(chain, i);
      for (let j = 0; j < N; j++) {
        if (j !== i && chain[j] === STATES.F && fibrilRunLength(chain, j) >= minRun
            && fibrilRunStart(chain, j) !== runI) {
          E -= jFF / 2;
        }
      }
    }
    if (chain[i] === STATES.D) {
      if (i > 0     && chain[i-1] === STATES.D) E -= jD / 2;
      if (i < N - 1 && chain[i+1] === STATES.D) E -= jD / 2;
    }
  }
  return E;
}

