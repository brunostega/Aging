# Aging

A 1D Ising-like Monte Carlo simulation of protein aggregation, running live in the browser.

Each site on a linear chain represents a protein monomer that can exist in one of three states: **Monomer**, **Disordered Oligomer**, or **Fibril**. The system evolves via the Metropolis algorithm, and you can watch nucleation, aggregation, and phase-like transitions unfold in real time.

## Live Demo

👉 [carlocamilloni.github.io/Aging](https://carlocamilloni.github.io/Aging)

---

## The Model

### States

| State | Color | Description |
|---|---|---|
| Monomer (M) | 🔵 Sky blue | Free, soluble protein |
| Disordered Oligomer (D) | 🟠 Orange | Loosely associated aggregate |
| Fibril (F) | 🩷 Pink-violet | Ordered amyloid-like aggregate — stabilised by cooperative interactions |

### Hamiltonian

```
H = Σᵢ εₛᵢ
  − J_D  Σ⟨i,j⟩  δ(D,D)
  − J_F  Σ⟨i,j⟩  δ(F,F) · 1[run ≥ minRun]
  − h_FF · N_F   · 1[exists active run]
```

- **εₛ** — intrinsic energy of each state
- **J_D** — disordered oligomer coupling: nearest neighbours only
- **J_F** — short-range fibril coupling: nearest neighbours, active only when both sites belong to a contiguous run of length ≥ `minRun`
- **h_FF** — background field acting on all fibril sites once at least one active run exists anywhere on the chain; models the catalytic effect of an established fibril nucleus on the rest of the chain

### Key physics

**Cooperative fibril nucleation** — fibril coupling is gated by a minimum run length (`minRun`). Isolated fibril sites gain no stabilisation energy, so nucleation requires building a critical seed before the fibril state becomes thermodynamically favourable. This mimics the nucleation barrier seen in amyloid formation.

**Catalytic background field** — once any active fibril run forms, `h_FF` lowers the energy of every fibril site on the chain regardless of its local context. This captures secondary nucleation and seeding: an established fibril template changes the free-energy landscape for the whole chain. The field is O(N) and avoids the O(N²) cost of explicit pairwise long-range couplings.

**Irreversible fibril mode** — an optional toggle locks any fibril site once it joins an active run, preventing reversion. This captures kinetic trapping of mature amyloid cores.

---

## Parameters

All energy parameters share the range **0–8**.

| Parameter | Description |
|---|---|
| Chain length N | Number of monomers (10–300) |
| k_B T | Temperature — controls the scale of thermal fluctuations |
| E_Monomer | Intrinsic energy of the monomer state |
| E_Disordered | Intrinsic energy of the disordered oligomer state |
| E_Fibril | Intrinsic energy of the fibril state |
| J_Disordered | Disordered–disordered nearest-neighbour coupling |
| J_Fibril | Short-range fibril–fibril coupling (nearest neighbours, run ≥ minRun) |
| h_FF background | Background field on all F sites once an active run exists |
| Min Fibril Run | Minimum run length to activate fibril coupling and the background field |

---

## Trajectory export

Clicking **Save (n)** downloads a JSON file containing every recorded snapshot:

```json
{
  "metadata": {
    "created": "2026-03-13T10:23:00.000Z",
    "software": "Aging",
    "n": 80,
    "totalSteps": 342,
    "snapshots": 342,
    "params": { "eM": 0, "eD": 1.5, "eF": 3.0, "jD": 1.2, "jF": 2.5, "hFF": 0.5, "minRun": 3, "T": 1.0 },
    "states": { "0": "Monomer", "1": "Disordered", "2": "Fibril" }
  },
  "trajectory": [
    { "step": 1, "chain": [0, 0, 2, 0, 1], "E": 12.4 },
    { "step": 2, "chain": [0, 1, 2, 0, 1], "E": 11.1 }
  ]
}
```

Each entry contains the step index, the full chain as an array of integers (0 = Monomer, 1 = Disordered, 2 = Fibril), and the total energy. The buffer is cleared on Reset.

---

## Running locally

You need [Node.js](https://nodejs.org) v18 or later.

```bash
git clone https://github.com/carlocamilloni/Aging.git
cd Aging
npm install
npm run dev
```

Then open [http://localhost:5173](http://localhost:5173).

### On macOS with MacPorts

```bash
sudo port install nodejs22 npm10
```

---

## Command-line interface

The simulation can be run headlessly from the terminal, writing the time trace and full trajectory to TSV files on the fly.

```bash
node aging-cli.js [options]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 1000 | Number of MC steps |
| `--n N` | 80 | Chain length |
| `--seq MMMFFF...` | — | Starting conformation string (overrides `--n`) |
| `--T N` | 1.0 | Temperature k_BT |
| `--eM N` | 0 | Intrinsic energy — Monomer |
| `--eD N` | 1.5 | Intrinsic energy — Disordered |
| `--eF N` | 3.0 | Intrinsic energy — Fibril |
| `--jD N` | 1.2 | J_Disordered coupling |
| `--jF N` | 2.5 | J_Fibril coupling |
| `--hFF N` | 0.5 | h_FF background field |
| `--minRun N` | 3 | Min fibril run length |
| `--irreversible` | off | Enable irreversible fibril locking |
| `--stopOnFibril` | off | Stop when chain reaches 100% fibril |
| `--trace FILE` | trace.tsv | Time trace output file |
| `--traj FILE` | trajectory.tsv | Trajectory output file |
| `--help` | — | Show usage |

Both output files are written incrementally — safe to inspect mid-run with `tail -f`.

**trace.tsv** — one row per step:
```
step  fM      fD      fF      fActiveF  E        avgE
0     1.0000  0.0000  0.0000  0.0000    80.0000  80.0000
1     0.9875  0.0125  0.0000  0.0000    79.5000  79.7500
...
```

**trajectory.tsv** — one row per step, chain encoded as a string of M/D/F characters:
```
step  chain           E
0     MMMMMMMM...     80.0000
1     MMDMMMMD...     79.5000
...
```

Example — run 5000 steps at low temperature starting from all monomers:
```bash
node aging-cli.js --steps 5000 --n 100 --T 0.5 --jF 3.0 --hFF 1.0 --stopOnFibril
```

Or via npm:
```bash
npm run cli -- --steps 5000 --T 0.5
```

---

## Tests

The energy functions and MC engine are covered by a suite of regression tests using [Vitest](https://vitest.dev). All tests are single-point energy calculations with hand-verified expected values — fully deterministic, no MC randomness.

```bash
npm test
```

Tests live in `tests/aging.test.js` and cover:

- `fibrilRunLength` — boundary cases and run detection
- `computeEnergy` — intrinsic energies, D-D coupling, F-F short-range coupling, `h_FF` background field (firing conditions, sub-threshold sites, multiple runs)
- `mcStep` — chain length preservation, irreversible locking, energy non-increase at low temperature

The CI workflow (`.github/workflows/test.yml`) runs the suite on every push and pull request.

---

## Building for production

```bash
npm run build
```

The static output lands in `dist/` and can be served from any static host.

---

## Deployment

Pushes to `main` deploy automatically to GitHub Pages via GitHub Actions. See [`.github/workflows/deploy.yml`](.github/workflows/deploy.yml).

After cloning, enable Pages in your repo under **Settings → Pages → Source → GitHub Actions**.

---

## Project structure

```
Aging/
├── .github/workflows/
│   ├── deploy.yml       # deploy to GitHub Pages on push to main
│   └── test.yml         # run test suite on push and pull request
├── src/
│   ├── model.js         # pure energy functions (computeEnergy, fibrilRunLength)
│   ├── simulation.js    # MC engine, helpers, constants
│   ├── aging.jsx        # React UI
│   └── main.jsx         # entry point
├── tests/
│   └── aging.test.js    # regression tests
├── aging-cli.js         # headless CLI runner
├── index.html
├── package.json
├── vite.config.js       # Vite + Vitest config
└── README.md
```

---

## Stack

- [React 18](https://react.dev) — UI and state management
- [Vite](https://vitejs.dev) + [Vitest](https://vitest.dev) — build tooling and testing
- Vanilla JS — Monte Carlo engine, no simulation libraries

---

## License

MIT © The Authors
