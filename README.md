# Aging

A 1D Ising-like Monte Carlo simulation of protein aggregation, running live in the browser.

Each site on a linear chain represents a protein monomer that can exist in one of three states: **Monomer**, **Disordered Oligomer**, or **Fibril**. The system evolves via the Metropolis algorithm, and you can watch nucleation, aggregation, and phase-like transitions unfold in real time.

## Live Demo

👉 [carlocamilloni.github.io/aging](https://carlocamilloni.github.io/aging)

---

## The Model

### States

| State | Color | Description |
|---|---|---|
| Monomer (M) | 🟢 Green | Free, soluble protein — lowest intrinsic energy |
| Disordered Oligomer (D) | 🟠 Orange | Loosely associated aggregate |
| Fibril (F) | 🔴 Red | Ordered amyloid-like aggregate — highest intrinsic cost, stabilised by interactions |

### Hamiltonian

```
H = Σᵢ εₛᵢ
  − J_F  Σ⟨i,j⟩    δ(F,F)  · 𝟙[run ≥ minRun]
  − J_FF Σ|i−j|>1  δ(F*,F*)
  − Σ|i−j|≤r  J̃ᵢⱼ δ(D,D)
```

- **εₛ** — intrinsic energy of each state (ε_M < ε_D < ε_F by default)
- **J_F** — short-range fibril coupling: direct neighbors only (|d| = 1), active only when both sites belong to a run of length ≥ `minRun`
- **J_FF** — long-range fibril coupling: between any two active fibril sites at |d| > 1
- **J̃ᵢⱼ ~ U[0, J_D max]** — quenched random disordered-oligomer coupling, drawn once per pair at initialisation and fixed for the lifetime of the simulation

### Key features

**Cooperative fibril nucleation** — fibril coupling is gated by a minimum run length (`minRun`). Isolated fibril sites gain no stabilisation energy, so nucleation requires building a critical seed before the fibril state becomes thermodynamically favourable. This mimics the nucleation barrier seen in amyloid formation.

**Quenched disorder on oligomer contacts** — each pair of disordered sites has its own fixed random coupling strength, reflecting heterogeneous sequence-dependent interactions in real disordered aggregates. The disorder matrix is resampled on every Reset.

**Irreversible fibril mode** — an optional toggle locks any fibril site once it joins an active run, preventing reversion. This captures the kinetic trapping of mature amyloid cores observed experimentally.

**Long-range fibril templating** — `J_FF` allows spatially separated fibril clusters to attract each other across the full chain, modelling secondary nucleation and seeding effects.

---

## Parameters

| Parameter | Description |
|---|---|
| Chain length N | Number of monomers (10–300) |
| k_B T | Temperature — controls the scale of thermal fluctuations |
| E_Monomer | Intrinsic energy of the monomer state |
| E_Disordered | Intrinsic energy of the disordered oligomer state |
| E_Fibril | Intrinsic energy of the fibril state |
| J_Fibril | Short-range fibril–fibril coupling (nearest neighbors) |
| J_Fibril long-range | Long-range fibril–fibril coupling (any distance) |
| J_D max | Maximum disordered coupling strength (quenched random per pair) |
| r_Disordered | Interaction range for disordered oligomers |
| Min Fibril Run | Minimum contiguous run length to activate fibril coupling |

---

## Running locally

You need [Node.js](https://nodejs.org) v18 or later.

```bash
git clone https://github.com/yourname/Aging.git
cd Aging
npm install
npm run dev
```

Then open [http://localhost:5173](http://localhost:5173).

### On macOS with MacPorts

```bash
sudo port install nodejs22 npm10
```

Then follow the steps above.

---

## Building for production

```bash
npm run build
```

The static output lands in `dist/` and can be served from any static host.

---

## Deployment

The repository is configured to deploy automatically to GitHub Pages via GitHub Actions on every push to `main`. See [`.github/workflows/deploy.yml`](.github/workflows/deploy.yml).

After cloning, enable Pages in your repo settings under **Settings → Pages → Source → GitHub Actions**.

---

## Project structure

```
Aging/
├── .github/
│   └── workflows/
│       └── deploy.yml   # CI/CD to GitHub Pages
├── src/
│   ├── aging.jsx         # simulation + UI
│   └── main.jsx         # React entry point
├── index.html
├── package.json
├── vite.config.js
└── README.md
```

---

## Stack

- [React 18](https://react.dev) — UI and state management
- [Vite](https://vitejs.dev) — build tooling
- Vanilla JS — Monte Carlo engine, no simulation libraries

---

## License

MIT © The Authors
