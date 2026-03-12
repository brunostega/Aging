// Pure model functions — no React dependencies.
// Imported by both the UI (fray.jsx) and the test suite.

export const STATES = { M: 0, D: 1, F: 2 };

export function fibrilRunLength(chain, i) {
  if (chain[i] !== STATES.F) return 0;
  const N = chain.length;
  let lo = i, hi = i;
  while (lo > 0 && chain[lo - 1] === STATES.F) lo--;
  while (hi < N - 1 && chain[hi + 1] === STATES.F) hi++;
  return hi - lo + 1;
}

export function computeEnergy(chain, params) {
  const { eM, eD, eF, jF, jFF, jD, minRun } = params;
  const N = chain.length;
  const baseE = [eM, eD, eF];
  let E = 0;
  for (let i = 0; i < N; i++) {
    E += baseE[chain[i]];
    if (chain[i] === STATES.F && fibrilRunLength(chain, i) >= minRun) {
      if (i > 0     && chain[i-1] === STATES.F && fibrilRunLength(chain, i-1) >= minRun) E -= jF / 2;
      if (i < N - 1 && chain[i+1] === STATES.F && fibrilRunLength(chain, i+1) >= minRun) E -= jF / 2;
      for (let j = 0; j < N; j++) {
        if (Math.abs(j - i) > 1 && chain[j] === STATES.F && fibrilRunLength(chain, j) >= minRun) E -= jFF / 2;
      }
    }
    if (chain[i] === STATES.D) {
      if (i > 0     && chain[i-1] === STATES.D) E -= jD / 2;
      if (i < N - 1 && chain[i+1] === STATES.D) E -= jD / 2;
    }
  }
  return E;
}

export function localEnergy(chain, idx, params) {
  const { eM, eD, eF, jF, jFF, jD, minRun } = params;
  const N = chain.length;
  const baseE = [eM, eD, eF];
  let E = baseE[chain[idx]];
  if (chain[idx] === STATES.F && fibrilRunLength(chain, idx) >= minRun) {
    if (idx > 0     && chain[idx-1] === STATES.F && fibrilRunLength(chain, idx-1) >= minRun) E -= jF;
    if (idx < N - 1 && chain[idx+1] === STATES.F && fibrilRunLength(chain, idx+1) >= minRun) E -= jF;
    for (let j = 0; j < N; j++) {
      if (Math.abs(j - idx) > 1 && chain[j] === STATES.F && fibrilRunLength(chain, j) >= minRun) E -= jFF;
    }
  }
  if (chain[idx] === STATES.D) {
    if (idx > 0     && chain[idx-1] === STATES.D) E -= jD;
    if (idx < N - 1 && chain[idx+1] === STATES.D) E -= jD;
  }
  return E;
}
