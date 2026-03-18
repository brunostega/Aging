// src/worker.js — MC simulation worker.
//
// The worker owns the simulation loop so the main thread stays responsive.
// Each batch runs SWEEPS_PER_BATCH sweeps, then posts results and yields via
// setTimeout(loop, 0) so the browser can process UI events between batches.
//
// Messages IN  (from UI):
//   { type: "start",  chain, locked, params, T, irreversible, stopOnFibril,
//                     E, step, sweepsPerBatch, generation }
//   { type: "stop" }
//   { type: "step",   chain, locked, params, T, irreversible, E, step, generation }
//   { type: "params", params, T, E }
//
// Messages OUT (to UI):
//   { type: "batch", chain, locked, E, step, snapshots: [{step, chain, E}], generation }
//   { type: "done",  generation }

import { STATES, mcStep } from "./simulation.js";

let running       = false;
let currentState  = null;

self.onmessage = (e) => {
  const msg = e.data;

  if (msg.type === "stop") {
    running      = false;
    currentState = null;
    return;
  }

  if (msg.type === "params") {
    // Hot-swap params/T without stopping — E is recomputed on the main thread
    if (currentState) {
      currentState.params = msg.params;
      currentState.T      = msg.T;
      currentState.E      = msg.E;
    }
    return;
  }

  if (msg.type === "step") {
    // Single-step from paused state
    const activeLocked = msg.irreversible ? msg.locked : null;
    const result = mcStep(msg.chain, msg.params, msg.T, activeLocked, msg.E);
    const step   = msg.step + 1;
    postMessage({
      type:      "batch",
      chain:     result.chain,
      locked:    result.locked,
      E:         result.E,
      step,
      snapshots: [{ step, chain: result.chain.slice(), E: result.E }],
      generation: msg.generation ?? 0,
    });
    return;
  }

  if (msg.type === "start") {
    running = true;
    const generation     = msg.generation ?? 0;
    const sweepsPerBatch = msg.sweepsPerBatch ?? 10;

    let state = {
      chain:        msg.chain,
      locked:       msg.locked,
      params:       msg.params,
      T:            msg.T,
      irreversible: msg.irreversible,
      stopOnFibril: msg.stopOnFibril,
      E:            msg.E,
      step:         msg.step,
    };
    currentState = state;

    function loop() {
      if (!running) return;

      const { params, T, irreversible, stopOnFibril } = state;
      let   { chain, locked, E, step }                = state;
      const snapshots = [];

      for (let i = 0; i < sweepsPerBatch; i++) {
        const activeLocked = irreversible ? locked : null;
        const result = mcStep(chain, params, T, activeLocked, E);
        chain  = result.chain;
        locked = result.locked;
        E      = result.E;
        step  += 1;
        snapshots.push({ step, chain: chain.slice(), E });

        // Stop early if the requested fibril fraction has been reached
        if (stopOnFibril) {
          let nF = 0;
          for (let j = 0; j < chain.length; j++) if (chain[j] === STATES.F) nF++;
          if (nF / chain.length >= 1.0) {
            running = false;
            postMessage({ type: "batch", chain, locked, E, step, snapshots, generation });
            postMessage({ type: "done",  generation });
            return;
          }
        }
      }

      state        = { ...state, chain, locked, E, step };
      currentState = state;
      postMessage({ type: "batch", chain, locked, E, step, snapshots, generation });
      setTimeout(loop, 0);
    }

    loop();
  }
};
