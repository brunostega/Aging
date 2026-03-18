// Simulation logic — MC step, state helpers, chain initialisation.
// No React dependencies.

import {
  STATES, computeEnergy,
  deltaEnergyFast, buildSimState, updateRunLen, countActiveRuns,
} from "./model.js";

export { STATES };

export const DEFAULT_N  = 100;
export const MIN_N      = 10;
export const MAX_N      = 10000;

export const STATE_NAMES  = ["Monomer", "Disordered", "Fibril"];
export const STATE_COLORS = ["#56B4E9", "#E69F00", "#CC79A7"];

export const DEFAULT_PARAMS = {
  eM: 0, eD: 1.4, eF: 6.0, jF: 3.5, hFF: 3.5, jD: 1.0, minRun: 3,
};

// Parse a string like "MMMFFFDDD" into a chain array.
// Returns { chain, error } — chain is null if the string is invalid.
export function parseChain(str) {
  const map = { M: STATES.M, D: STATES.D, F: STATES.F };
  const upper = str.toUpperCase().trim();
  if (upper.length === 0) return { chain: null, error: "String is empty" };
  if (upper.length < MIN_N) return { chain: null, error: `Too short — minimum ${MIN_N} sites` };
  if (upper.length > MAX_N) return { chain: null, error: `Too long — maximum ${MAX_N} sites` };
  const chain = [];
  for (const ch of upper) {
    if (!(ch in map)) return { chain: null, error: `Invalid character '${ch}' — use M, D, F only` };
    chain.push(map[ch]);
  }
  return { chain, error: null };
}

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

// mcStep accepts the current global energy and returns the updated energy.
// Builds runLen cache once per sweep — O(N) — then uses O(window) deltaEnergyFast
// per flip instead of O(N) fibrilRunLength calls or full recomputes.
export function mcStep(chain, params, T, locked, currentE) {
  const N = chain.length;
  const c = [...chain];
  let E = currentE;

  let { nF, nActiveRuns, runLen } = buildSimState(c, params.minRun);

  for (let _ = 0; _ < N; _++) {
    const idx = Math.floor(Math.random() * N);
    if (locked && locked[idx]) continue;
    const oldState = c[idx];
    const newState = Math.floor(Math.random() * 3);
    if (newState === oldState) continue;

    c[idx] = newState;
    // Only update runLen when F is involved — otherwise runLen is unchanged
    const fInvolved = oldState === STATES.F || newState === STATES.F;
    if (fInvolved) updateRunLen(c, runLen, idx);

    const dE = deltaEnergyFast(c, idx, oldState, newState, params, nF, nActiveRuns, runLen);

    if (dE > 0 && Math.random() >= Math.exp(-dE / T)) {
      // Reject — revert chain and runLen
      c[idx] = oldState;
      if (fInvolved) updateRunLen(c, runLen, idx);
    } else {
      // Accept — update incremental counters
      E += dE;
      if (oldState === STATES.F) nF--;
      if (newState === STATES.F) nF++;
      if (fInvolved) nActiveRuns = countActiveRuns(runLen, params.minRun);
    }
  }

  // Lock any newly active fibril sites when irreversible mode is on.
  // runLen is up to date so no need to call fibrilRunLength.
  const newLocked = locked ? [...locked] : new Array(N).fill(false);
  if (locked !== null) {
    for (let i = 0; i < N; i++) {
      if (!newLocked[i] && c[i] === STATES.F && runLen[i] >= params.minRun) {
        newLocked[i] = true;
      }
    }
  }

  return { chain: c, locked: newLocked, E };
}
