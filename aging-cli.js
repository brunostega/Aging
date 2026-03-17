#!/usr/bin/env node
// aging-cli.js — command-line interface for the Aging simulation.
//
// Usage:
//   node aging-cli.js [options]
//
// Options:
//   --steps N          Number of MC steps to run (default: 1000)
//   --n N              Chain length (default: 80)
//   --seq MMMFFF...    Starting conformation string (overrides --n)
//   --eM N             Intrinsic energy — Monomer       (default: 0)
//   --eD N             Intrinsic energy — Disordered    (default: 1.5)
//   --eF N             Intrinsic energy — Fibril        (default: 3.0)
//   --jD N             J_Disordered coupling            (default: 1.2)
//   --jF N             J_Fibril coupling                (default: 2.5)
//   --hFF N            h_FF background field            (default: 0.5)
//   --minRun N         Min fibril run length            (default: 3)
//   --T N              Temperature k_BT                 (default: 1.0)
//   --irreversible     Enable irreversible fibril locking
//   --stopOnFibril     Stop when chain is 100% fibril
//   --trace FILE       Output file for time trace       (default: trace.tsv)
//   --traj FILE        Output file for trajectory       (default: trajectory.tsv)
//   --help             Show this help message

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

import {
  STATES, STATE_NAMES, DEFAULT_N, DEFAULT_PARAMS, MIN_N, MAX_N,
  initChain, parseChain, countStates, mcStep,
} from "./src/simulation.js";
import { computeEnergy, fibrilRunLength } from "./src/model.js";

// ── Argument parsing ────────────────────────────────────────────────────────

function parseArgs(argv) {
  const args = argv.slice(2);
  const opts = {
    steps:         1000,
    n:             DEFAULT_N,
    seq:           null,
    T:             1.0,
    irreversible:  false,
    stopOnFibril:  false,
    trace:         "trace.tsv",
    traj:          "trajectory.tsv",
    params:        { ...DEFAULT_PARAMS },
  };

  const paramKeys = ["eM", "eD", "eF", "jD", "jF", "hFF", "minRun"];

  for (let i = 0; i < args.length; i++) {
    const a = args[i];
    if (a === "--help" || a === "-h") {
      const src = fs.readFileSync(fileURLToPath(import.meta.url), "utf8");
      const helpLines = src.split("\n")
        .filter(l => l.startsWith("//"))
        .map(l => l.slice(3));
      console.log(helpLines.join("\n"));
      process.exit(0);
    }
    else if (a === "--irreversible")  opts.irreversible = true;
    else if (a === "--stopOnFibril")  opts.stopOnFibril = true;
    else if (a === "--steps")   opts.steps  = parseInt(args[++i]);
    else if (a === "--n")       opts.n      = parseInt(args[++i]);
    else if (a === "--T")       opts.T      = parseFloat(args[++i]);
    else if (a === "--seq")     opts.seq    = args[++i];
    else if (a === "--trace")   opts.trace  = args[++i];
    else if (a === "--traj")    opts.traj   = args[++i];
    else if (paramKeys.some(k => a === `--${k}`)) {
      const key = a.slice(2);
      opts.params[key] = key === "minRun" ? parseInt(args[++i]) : parseFloat(args[++i]);
    }
    else {
      console.error(`Unknown option: ${a}. Use --help for usage.`);
      process.exit(1);
    }
  }
  return opts;
}

// ── Active fibril fraction helper ───────────────────────────────────────────

function activeFibrilFrac(chain, minRun) {
  const nF = chain.filter(s => s === STATES.F).length;
  if (nF === 0) return 0;
  let active = 0;
  const N = chain.length;
  let i = 0;
  while (i < N) {
    if (chain[i] === STATES.F) {
      let j = i;
      while (j < N && chain[j] === STATES.F) j++;
      const len = j - i;
      if (len >= minRun) active += len;
      i = j;
    } else i++;
  }
  return active / nF;
}

// ── Main ────────────────────────────────────────────────────────────────────

const opts = parseArgs(process.argv);
const { steps, T, params, irreversible, stopOnFibril } = opts;

// Initialise chain
let chain;
if (opts.seq) {
  const { chain: parsed, error } = parseChain(opts.seq);
  if (error) { console.error(`--seq error: ${error}`); process.exit(1); }
  chain = parsed;
} else {
  const n = Math.max(MIN_N, Math.min(MAX_N, opts.n));
  chain = initChain(n);
}

let locked = new Array(chain.length).fill(false);
let E = computeEnergy(chain, params);
let energySum = E;

// Open output files — write headers
const traceFd  = fs.openSync(opts.trace, "w");
const trajFd   = fs.openSync(opts.traj,  "w");

fs.writeSync(traceFd, ["step", "fM", "fD", "fF", "fActiveF", "E", "avgE"].join("\t") + "\n");
fs.writeSync(trajFd,  ["step", "chain", "E"].join("\t") + "\n");

// Write step 0
function writeStep(step) {
  const counts = countStates(chain);
  const N = chain.length;
  const fM = counts[0] / N, fD = counts[1] / N, fF = counts[2] / N;
  const fAF = activeFibrilFrac(chain, params.minRun);
  const avgE = energySum / (step + 1);
  fs.writeSync(traceFd,
    [step, fM.toFixed(4), fD.toFixed(4), fF.toFixed(4), fAF.toFixed(4),
     E.toFixed(4), avgE.toFixed(4)].join("\t") + "\n"
  );
  fs.writeSync(trajFd,
    [step, chain.map(s => STATE_NAMES[s][0]).join(""), E.toFixed(4)].join("\t") + "\n"
  );
}

writeStep(0);

// Print config summary
console.log("Aging CLI");
console.log(`  Chain length : ${chain.length}`);
console.log(`  Steps        : ${steps}`);
console.log(`  T            : ${T}`);
console.log(`  Params       : ${JSON.stringify(params)}`);
console.log(`  Irreversible : ${irreversible}`);
console.log(`  Stop@100%F   : ${stopOnFibril}`);
console.log(`  Trace file   : ${opts.trace}`);
console.log(`  Traj file    : ${opts.traj}`);
console.log("");

// Progress bar helper
const BAR_WIDTH = 40;
function progress(step) {
  const pct = step / steps;
  const filled = Math.round(pct * BAR_WIDTH);
  const bar = "█".repeat(filled) + "░".repeat(BAR_WIDTH - filled);
  process.stdout.write(`\r  [${bar}] ${step}/${steps}`);
}

// ── MC loop ─────────────────────────────────────────────────────────────────

for (let s = 1; s <= steps; s++) {
  const activeLocked = irreversible ? locked : null;
  const result = mcStep(chain, params, T, activeLocked, E);
  chain  = result.chain;
  locked = result.locked;
  E      = result.E;
  energySum += E;

  writeStep(s);

  if (s % 100 === 0 || s === steps) progress(s);

  if (stopOnFibril && chain.every(s => s === STATES.F)) {
    progress(s);
    console.log(`\n  Stopped early: 100% fibril at step ${s}`);
    break;
  }
}

fs.closeSync(traceFd);
fs.closeSync(trajFd);

console.log(`\n  Done. Trace → ${opts.trace}  |  Trajectory → ${opts.traj}`);
