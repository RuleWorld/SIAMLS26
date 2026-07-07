# Structure-based Modeling of Biochemical Reaction Networks using BioNetGen

**A minitutorial at SIAM Life Sciences 2026**

This page has everything you need before, during, and after the tutorial: a schedule, software setup instructions, and links to all example models we'll walk through.

> **Note:** This is not a hands-on coding session — you don't need to install anything to follow along. But if you'd like to try the examples yourself afterward (or follow along live), see [Getting Started](#getting-started) below.

---

## Schedule

| Topic | Presenter |
|-------|-----------|
| Introduction to rule-based modeling, BioNetGen VSCode extension, SBML | James Faeder |
| Demonstration in VSCode, Interfacing with MATLAB | Laura Strube (LFStrube@gmail.com) |
| PyBioNetGen, Stochastic Simulation, Sensitivity Analysis, Parameter Estimation | James Faeder |
| Network-based vs. network-free simulation, NFsim | Alex DiBiasi |
| Model integration and simulation pipelines with BNGPlayground (INDRA, EGFR case study) | Achyudhan Kutuva |
| Applications: WESTPA/WEBNG, RuleHub, and the broader ecosystem | Various |

---

## What we'll cover

- **Rule-based modeling concepts** — why some biochemical systems (like multivalent ligand-receptor binding) can't be described with a fixed set of ODEs
- **BioNetGen** — building rule-based models with `.bngl` files, generating reaction networks, and simulating with ODE/SSA solvers
- **Python ecosystem** — PyBioNetGen, Jupyter-based workflows, sensitivity analysis, parameter estimation, identifiability
- **NFsim** — simulating systems whose reaction networks are too large (or infinite) to enumerate, using network-free simulation
- **Applications** — weighted ensemble simulation (WESTPA/WEBNG), model integration with INDRA, and the broader RuleHub model repository

---

## Getting started

To follow along with the example models (optional, not required for the talk itself), consider one of avenues for accessing BioNetGen:

1. **BNGPlayground** Simply click this [link](https://ruleworld.github.io/bngplayground/)

2. **Install BioNetGen.** Follow the instructions from the [BioNetGen website](http://bionetgen.readthedocs.io/en/latest/install.html) or install the [BioNetGen VSCode extension](https://marketplace.visualstudio.com/items?itemName=BioNetGen-ext.biontegen-vscode-extension) for an integrated editing/simulation experience.
3. **Install PyBioNetGen** for Python-based workflows:
   ```bash
   pip install bionetgen
   ```

---

## Example models

All example `.bngl` files discussed in the tutorial are in the [`models/`](./models) folder of this repo:

- `Example_SimpleSTAT.bngl`— Simple model of IL10-induced STAT phosphorylation
- `LR.bngl` — basic L + R ⇌ LR ligand-receptor model
- `TLBR.bngl` — trivalent ligand / bivalent receptor model (NFsim)
- `LAT.bngl` — LAT-Grb2-SOS1 aggregation model (NFsim)
- `translation.bngl` — ribosome translation model (NFsim)

---

## Resources

- [BioNetGen documentation](https://bionetgen.org)
- [BNGPlayground](https://ruleworld.github.io/bngplayground/)
- [RuleHub](https://github.com/RuleWorld/RuleHub) — repository of curated rule-based models
- [NFsim](https://github.com/RuleWorld/nfsim)
- Slides: *(link once finalized)*

---

## Questions?

We'll have time for Q&A throughout the session — feel free to interrupt with questions. After the tutorial, you can reach us at *[bionetgen.main@gmail.com]* or open an issue on this repo.

---

*Faeder Lab, University of Pittsburgh*
