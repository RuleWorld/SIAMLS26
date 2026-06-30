---
title: Rule-Based Modeling Minitutorial — SIAM LS26
---

# Rule-Based Modeling with BioNetGen and NFsim

**A minitutorial at SIAM Life Sciences 2026 (SIAM LS26)**

This page has everything you need before, during, and after the tutorial: a schedule, software setup instructions, and links to all example models we'll walk through.

> **Note:** This is not a hands-on coding session — you don't need to install anything to follow along. But if you'd like to try the examples yourself afterward (or follow along live), see [Getting Started](#getting-started) below.

---

## Schedule

| Topic | Presenter |
|-------|-----------|
| Introduction to rule-based modeling, BioNetGen VSCode extension | James Faeder |
| Network-based vs. network-free simulation, NFsim | Alex DiBiasi |
| Python tools for rule-based modeling (PyBioNetGen, sensitivity analysis, parameter estimation) | Laura Strube |
| Model integration and simulation pipelines (INDRA, EGFR case study) | Achyudhan Kutuva |
| Applications: WESTPA/WEBNG, RuleHub, and the broader ecosystem | Various |

*(We'll fill in exact times and last names before this goes live — for now this just needs to look complete enough to discuss with Jim.)*

---

## What we'll cover

- **Rule-based modeling concepts** — why some biochemical systems (like multivalent ligand-receptor binding) can't be described with a fixed set of ODEs
- **BioNetGen** — building rule-based models with `.bngl` files, generating reaction networks, and simulating with ODE/SSA solvers
- **NFsim** — simulating systems whose reaction networks are too large (or infinite) to enumerate, using network-free simulation
- **Python ecosystem** — PyBioNetGen, Jupyter-based workflows, sensitivity analysis, parameter estimation, identifiability
- **Applications** — weighted ensemble simulation (WESTPA/WEBNG), model integration with INDRA, and the broader RuleHub model repository

---

## Getting started

To follow along with the example models (optional, not required for the talk itself):

1. **Install BioNetGen.** Download from the [BioNetGen website](https://bionetgen.org) or install the [BioNetGen VSCode extension](#) for an integrated editing/simulation experience.
2. **(Optional) Install PyBioNetGen** for Python-based workflows:
   ```bash
   pip install bionetgen
   ```
3. **(Optional) Install NFsim** if you want to run the network-free examples directly — see the [NFsim GitHub repo](https://github.com/RuleWorld/nfsim) for build instructions.

---

## Example models

All example `.bngl` files discussed in the tutorial are in the [`models/`](./models) folder of this repo:

- `one_arm_model.bngl` — minimal model for learning the BNGL syntax
- `lr_simple.bngl` — basic L + R ⇌ LR ligand-receptor model
- `tlbr.bngl` — trivalent ligand / bivalent receptor model (network-free example)
- `lat_aggregation.bngl` — LAT-Grb2-SOS1 aggregation model
- `translation_model.bngl` — ribosome translation model (network-free example)
- `egfr.bngl` — EGFR pathway model (INDRA case study)

---

## Resources

- [BioNetGen documentation](https://bionetgen.org)
- [RuleHub](https://github.com/RuleWorld/RuleHub) — repository of curated rule-based models
- [NFsim](https://github.com/RuleWorld/nfsim)
- [WESTPA](https://westpa.github.io/westpa/)
- Slides: *(link once finalized)*

---

## Questions?

We'll have time for Q&A throughout the session — feel free to interrupt with questions. After the tutorial, you can reach us at *[lab contact email]* or open an issue on this repo.

---

*Faeder Lab, University of Pittsburgh*