// MC engine, chain helpers, and constants.
// No React dependencies.

import {
  STATES,
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

// Returns { chain, error } — chain is null if the string is invalid.
// Accepts any mix of upper/lower case M, D, F.
export function parseChain(str) {
  const map = { M: STATES.M, D: STATES.D, F: STATES.F };
  const upper = str.toUpperCase().trim();
  if (upper.length === 0)     return { chain: null, error: "String is empty" };
  if (upper.length < MIN_N)   return { chain: null, error: `Too short — minimum ${MIN_N} sites` };
  if (upper.length > MAX_N)   return { chain: null, error: `Too long — maximum ${MAX_N} sites` };
  const chain = [];
  for (const ch of upper) {
    if (!(ch in map)) return { chain: null, error: `Invalid character '${ch}' — use M, D, F only` };
    chain.push(map[ch]);
  }
  return { chain, error: null };
}

// Returns an all-monomer chain of length n.
export function initChain(n) {
  return new Array(n).fill(STATES.M);
}

// Returns [nM, nD, nF].
export function countStates(chain) {
  const counts = [0, 0, 0];
  for (let i = 0; i < chain.length; i++) counts[chain[i]]++;
  return counts;
}

// Returns an array of run lengths for each contiguous F block.
export function fibrilRunLengths(chain) {
  const N = chain.length;
  const runs = [];
  let i = 0;
  while (i < N) {
    if (chain[i] === STATES.F) {
      let j = i;
      while (j < N && chain[j] === STATES.F) j++;
      runs.push(j - i);
      i = j;
    } else {
      i++;
    }
  }
  return runs;
}

// One MC sweep of N attempted single-site flips using the Metropolis criterion.
//
// Performance notes:
//   - buildSimState builds runLen, nF, nActiveRuns in one O(N) pass per sweep.
//   - M↔D flips use a fast O(1) path in deltaEnergyFast (no runLen update).
//   - F-involved flips update runLen O(runLength) and recount active runs.
//   - locked sites are skipped without any energy computation.
//
// Parameters:
//   chain    — input chain (not mutated)
//   params   — model parameters
//   T        — temperature k_B T
//   locked   — boolean array of locked sites, or null to disable locking
//   currentE — current energy (avoids recomputing it from scratch)
//
// Returns { chain, locked, E } with the updated state.
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

    const fInvolved = oldState === STATES.F || newState === STATES.F;

    c[idx] = newState;
    if (fInvolved) updateRunLen(c, runLen, idx);

    const dE = deltaEnergyFast(c, idx, oldState, newState, params, nF, nActiveRuns, runLen);

    if (dE > 0 && Math.random() >= Math.exp(-dE / T)) {
      // Reject — restore chain and runLen
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

  // In irreversible mode, lock any site that just became part of an active run.
  // runLen is already up to date so no extra traversal needed.
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
