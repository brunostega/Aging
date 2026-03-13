// Simulation logic — MC step, state helpers, chain initialisation.
// No React dependencies.

import { STATES, fibrilRunLength, computeEnergy } from "./model.js";

export { STATES };

export const DEFAULT_N  = 80;
export const MIN_N      = 10;
export const MAX_N      = 300;

export const STATE_NAMES  = ["Monomer", "Disordered", "Fibril"];
export const STATE_COLORS = ["#56B4E9", "#E69F00", "#CC79A7"];

export const DEFAULT_PARAMS = {
  eM: 0, eD: 1.5, eF: 3.0, jF: 2.5, jFF: 0.5, jD: 1.2, minRun: 3,
};

export function initChain(n) {
  return new Array(n).fill(STATES.M);
}

export function countStates(chain) {
  const counts = [0, 0, 0];
  chain.forEach(s => counts[s]++);
  return counts;
}

export function fibrilRunLengths(chain) {
  const n = chain.length;
  const runs = [];
  let i = 0;
  while (i < n) {
    if (chain[i] === STATES.F) {
      let j = i;
      while (j < n && chain[j] === STATES.F) j++;
      runs.push(j - i);
      i = j;
    } else {
      i++;
    }
  }
  return runs;
}

export function mcStep(chain, params, T, locked) {
  const N = chain.length;
  const c = [...chain];

  for (let _ = 0; _ < N; _++) {
    const idx = Math.floor(Math.random() * N);
    if (locked && locked[idx]) continue;
    const oldState = c[idx];
    const newState = Math.floor(Math.random() * 3);
    if (newState === oldState) continue;
    //const oldE = localEnergy(c, idx, params);
    const oldE = computeEnergy(c, params);
    c[idx] = newState;
    //const newE = localEnergy(c, idx, params);
    const newE = computeEnergy(c, params);
    const dE = newE - oldE;
    if (dE > 0 && Math.random() >= Math.exp(-dE / T)) c[idx] = oldState;
  }

  // Lock any newly active fibril sites when irreversible mode is on
  const newLocked = locked ? [...locked] : new Array(N).fill(false);
  if (locked !== null) {
    for (let i = 0; i < N; i++) {
      if (!newLocked[i] && c[i] === STATES.F && fibrilRunLength(c, i) >= params.minRun) {
        newLocked[i] = true;
      }
    }
  }

  return { chain: c, locked: newLocked };
}
