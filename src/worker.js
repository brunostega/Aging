// src/worker.js — MC simulation worker.
// Runs in a separate thread. Communicates with App.jsx via postMessage.
//
// Messages IN (from UI):
//   { type: "start",  chain, locked, params, T, irreversible, stopOnFibril,
//                     E, step, energySum, sweepsPerBatch }
//   { type: "stop" }
//   { type: "step",   chain, locked, params, T, irreversible, E, step, energySum }
//   { type: "config", sweepsPerBatch }
//   { type: "params", params, T, E }
//
// Messages OUT (to UI):
//   { type: "batch",  chain, locked, E, step, energySum,
//                     snapshots: [{ step, chain, E }, ...] }
//   { type: "done" }   — emitted after final batch when stopOnFibril triggers

import { STATES, mcStep } from "./simulation.js";

let running = false;
let sweepsPerBatch = 10;
let currentState = null;  // shared with params message handler

self.onmessage = (e) => {
  const msg = e.data;

  if (msg.type === "stop") {
    running = false;
    currentState = null;
    return;
  }

  if (msg.type === "config") {
    sweepsPerBatch = msg.sweepsPerBatch;
    return;
  }

  if (msg.type === "params") {
    // Update params and T in the running state so the next batch picks them up
    if (currentState) {
      currentState.params = msg.params;
      currentState.T      = msg.T;
      currentState.E      = msg.E;  // energy was recomputed on the JS side
    }
    return;
  }

  if (msg.type === "step") {
    const { params, T, irreversible } = msg;
    const activeLocked = irreversible ? msg.locked : null;
    const result = mcStep(msg.chain, params, T, activeLocked, msg.E);
    const step = msg.step + 1;
    const energySum = msg.energySum + result.E;
    postMessage({
      type: "batch",
      chain: result.chain,
      locked: result.locked,
      E: result.E,
      step,
      energySum,
      snapshots: [{ step, chain: result.chain.slice(), E: result.E }],
    });
    return;
  }

  if (msg.type === "start") {
    running = true;
    if (msg.sweepsPerBatch) sweepsPerBatch = msg.sweepsPerBatch;

    let state = {
      chain:        msg.chain,
      locked:       msg.locked,
      params:       msg.params,
      T:            msg.T,
      irreversible: msg.irreversible,
      stopOnFibril: msg.stopOnFibril,
      E:            msg.E,
      step:         msg.step,
      energySum:    msg.energySum,
    };
    currentState = state;

    function loop() {
      if (!running) return;

      const { params, T, irreversible, stopOnFibril } = state;
      let { chain, locked, E, step, energySum } = state;
      const snapshots = [];

      for (let i = 0; i < sweepsPerBatch; i++) {
        const activeLocked = irreversible ? locked : null;
        const result = mcStep(chain, params, T, activeLocked, E);
        chain  = result.chain;
        locked = result.locked;
        E      = result.E;
        step  += 1;
        energySum += E;
        snapshots.push({ step, chain: chain.slice(), E });

        if (stopOnFibril && chain.every(s => s === STATES.F)) {
          running = false;
          postMessage({ type: "batch", chain, locked, E, step, energySum, snapshots });
          postMessage({ type: "done" });
          return;
        }
      }

      state = { ...state, chain, locked, E, step, energySum };
      currentState = state;
      postMessage({ type: "batch", chain, locked, E, step, energySum, snapshots });
      setTimeout(loop, 0); // yield so incoming "stop" messages can be processed
    }

    loop();
  }
};
