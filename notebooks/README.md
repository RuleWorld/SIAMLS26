# BioNetGen Parameter Estimation Tutorials

This directory contains Jupyter notebooks demonstrating parameter estimation and sensitivity analysis for rule-based models using BioNetGen.

---

## Notebooks

### 1. `parameter_estimation_demo.ipynb` ⭐ **Main Tutorial**

**Purpose**: Complete end-to-end workflow for parameter estimation in BioNetGen models.

**What it demonstrates**:
- Synthetic data generation with realistic noise models
- PETAb standardized format for parameter estimation problems
- Multi-trial global optimization using differential evolution
- Parameter identifiability analysis using Fisher Information Matrix
- Automated diagnostics for common parameter estimation issues

**Dependencies**:
- **Model file**: `fceri_starter.bngl` (FcεRI receptor signaling model)
- **Python modules** (in `python/` folder):
  - `bngl_simulator.py` - BioNetGen simulation wrapper
  - `parameter_estimator.py` - Parameter estimation framework
  - `debug_fixed_params.py` - Diagnostics for fixed parameter handling
  - `identifiability_diagnosis.py` - Sensitivity and correlation analysis
  - `fim_analysis_fixed.py` - Fisher Information Matrix computation
- **Documentation**:
  - `PARAMETER_FIXING_GUIDE.md` - How to handle non-identifiable parameters
  - `IDENTIFIABILITY_HELP_README.md` - Quick start guide for identifiability

**Key sections**:
1. **Setup and Imports** - Configure environment
2. **Model Simulation** - Load and test BNGL model
3. **Synthetic Data Generation** - Create noisy measurements
4. **PETAb File Creation** - Standardized problem specification
5. **Parameter Estimation** - Multi-trial optimization
6. **Visualize Results** - Assess fit quality
7. **Identifiability Analysis** - Determine which parameters are estimable

**Configuration cells** (modify these for your model):
- Section 3.1: Experimental design (observables, conditions, timepoints)
- Section 4.1: Parameters to estimate (bounds, initial guesses)
- Section 4.2: Observable formulas (scaling, transformations)
- Section 5.3: Optimization settings (trials, iterations)

**Diagnostic cells** (run when troubleshooting):
- Section 5.1.1: Debug fixed parameter issues
- Section 7.1: Diagnose identifiability problems
- Section 7.3: Fix ill-conditioned FIM

---

### 2. `bngl_sensitivity_analysis.ipynb`

**Purpose**: Parameter sensitivity analysis for BioNetGen models.

**What it demonstrates**:
- Local sensitivity analysis (one-at-a-time parameter perturbations)
- Global sensitivity using Monte Carlo sampling
- Visualization of parameter effects on observables
- Time-dependent sensitivity profiles

**Dependencies**:
- **Model file**: Any BNGL model
- **Python modules**:
  - `bngl_simulator.py` - BioNetGen simulation wrapper
- **External packages**: numpy, pandas, matplotlib, scipy

**Key analyses**:
- Normalized sensitivity coefficients
- Parameter ranking by influence
- Observable-specific sensitivity
- Time-course sensitivity profiles

---

### 3. `run_bngl.ipynb`

**Purpose**: Basic BioNetGen model simulation and exploration.

**What it demonstrates**:
- Loading and running BNGL models
- Parameter manipulation
- Time-course simulations
- Basic plotting of observables

**Dependencies**:
- **Model file**: `fceri_starter.bngl` or any BNGL model
- **Python modules**:
  - `bngl_simulator.py` - BioNetGen simulation wrapper
- **External packages**: numpy, pandas, matplotlib

**Use cases**:
- Testing new BNGL models
- Quick parameter exploration
- Generating reference trajectories
- Model validation

---

### 4. `RasModelAnalyzer.ipynb`

**Purpose**: Analysis of the Ras signaling pathway model.

**What it demonstrates**:
- Model-specific analysis workflows
- Application of simulation and analysis to Ras signaling
- Custom analysis for specific biological questions

**Dependencies**:
- **Model file**: `RasModel.bngl`
- **Python modules**: Various (model-specific)

---

## Python Modules (`python/` folder)

### Core Simulation and Estimation

1. **`bngl_simulator.py`** ⭐ **Core**
   - Wrapper for BioNetGen command-line interface
   - Handles BNGL model loading, parameter setting, simulation
   - Provides convenient interface for time-course simulations
   - Used by: All notebooks

2. **`parameter_estimator.py`** ⭐ **Core**
   - PETAb-based parameter estimation framework
   - Negative log-likelihood computation
   - Parameter vector/dict conversions
   - Supports fixed and estimated parameters
   - Used by: `parameter_estimation_demo.ipynb`

### Diagnostic Tools

3. **`debug_fixed_params.py`**
   - Diagnoses mismatches between fixed parameter nominal values and true values
   - Provides quick-fix for estimator configuration issues
   - Prevents NLLH calculation errors
   - Used by: `parameter_estimation_demo.ipynb` (Section 5.1.1)

4. **`identifiability_diagnosis.py`**
   - Parameter sensitivity analysis
   - Parameter correlation estimation
   - Automated recommendations for which parameters to fix
   - Comprehensive visualization
   - Used by: `parameter_estimation_demo.ipynb` (Section 7.1)

5. **`fim_analysis_fixed.py`**
   - Corrected Fisher Information Matrix computation
   - Excludes noise parameters for numerical stability
   - Proper parameter space transformations (log10 ↔ natural)
   - Uncertainty quantification via Cramér-Rao bound
   - Used by: `parameter_estimation_demo.ipynb` (Section 7.2)

### Experimental Tools (Not Currently Used in Notebooks)

6. **`diagnose_nllh.py`**
   - Decomposes NLLH into residual, sigma, and constant terms
   - Helps understand why NLLH behaves unexpectedly
   - Identifies when optimizer is "gaming" likelihood via sigma

7. **`diagnose_sigma_issue.py`**
   - Explains why sigma is overestimated
   - Decomposes residual variance into noise vs model inadequacy
   - Diagnostic for identifiability issues

8. **`identifiability_helper_cell.py`**
   - Helper functions for identifiability analysis
   - May contain utility functions used by other modules

---

## Model Files

### `fceri_starter.bngl`
FcεRI receptor signaling model - used as the primary example in `parameter_estimation_demo.ipynb` and `run_bngl.ipynb`.

**Parameters** (11 total):
- Initial concentrations: `Lig_tot`, `Rec_tot`, `Lyn_tot`
- Rate constants: `kp1`, `km1`, `kp2`, `km2`, `kpL`, `kmL`, `pLb`, `dm`

**Observables** (6 total):
- `RecDimer1`, `RecDimer2` - Receptor dimerization
- `LynFree`, `LynBoundU` - Lyn kinase states
- `BetaP1`, `BetaP2` - Receptor phosphorylation

### `fceri_starter_sbml.xml`
SBML export of the FcεRI model (generated from BNGL).

### `RasModel.bngl`
Ras signaling pathway model - used in `RasModelAnalyzer.ipynb`.

---

## Documentation Files

### `PARAMETER_FIXING_GUIDE.md` ⭐ **Essential**
Comprehensive guide on handling non-identifiable parameters.

**Contents**:
- Why parameters become non-identifiable
- Strategies for fixing parameters
- How to choose which parameters to fix
- Step-by-step workflow
- Model-specific recommendations for FcεRI

**When to read**: Before running parameter estimation with limited data.

### `IDENTIFIABILITY_HELP_README.md`
Quick-start guide for identifiability analysis.

**Contents**:
- How to use Section 7.1 diagnostics
- Interpreting sensitivity and correlation results
- Applying automated recommendations
- Troubleshooting guide

**When to read**: When encountering FIM issues or poor parameter recovery.

---

## Output Directories

### `petab_files/`
PETAb standardized files generated by `parameter_estimation_demo.ipynb`:
- `conditions.tsv` - Experimental conditions
- `measurements.tsv` - Measurement data
- `observables.tsv` - Observable definitions
- `parameters.tsv` - Parameter specifications

### `estimation_results/`
Results from parameter estimation runs:
- `trial_results_TIMESTAMP.json` - All optimization trials
- `best_parameters_TIMESTAMP.csv` - Best parameter estimates
- `best_fit_TIMESTAMP.png` - Fit visualization
- `trial_comparison_TIMESTAMP.png` - Multi-trial comparison
- `fim_analysis_fixed_TIMESTAMP.png` - FIM analysis plots
- `param_correlation_TIMESTAMP.png` - Parameter correlations
- `identifiability_diagnosis_TIMESTAMP.png` - Identifiability analysis

### `results_fceri_starter/`
Legacy results directory from `run_bngl.ipynb`.

### `results_RasModel/`
Results directory for `RasModelAnalyzer.ipynb`.

---

## Quick Start Guide

### For First-Time Users

1. **Start with basic simulation**: Run `run_bngl.ipynb` to understand how BNGL models work
2. **Learn parameter estimation**: Open `parameter_estimation_demo.ipynb` and follow sections 1-7
3. **Handle identifiability issues**: Use Section 7.1 diagnostics to fix problematic parameters
4. **Explore sensitivity**: Run `bngl_sensitivity_analysis.ipynb` to understand parameter effects

### For Parameter Estimation

1. Open `parameter_estimation_demo.ipynb`
2. **Modify Section 3.1**: Define your experimental design
   - Which observables to measure
   - Experimental conditions
   - Measurement timepoints
3. **Modify Section 4.1**: Configure parameters
   - Set `estimate=1` for parameters to fit
   - Set `estimate=0` for fixed parameters
   - Define bounds and initial guesses
4. **Run Sections 1-5**: Generate data and optimize
5. **Check Section 7.1**: Diagnose identifiability
6. **Iterate**: Fix parameters based on recommendations and re-run

### Common Workflows

**Workflow 1: Estimate all parameters (first attempt)**
- Keep all parameters with `estimate=1` in Section 4.1
- Run full workflow
- Section 7.1 will tell you which parameters to fix

**Workflow 2: Fixed parameter workflow**
- Set some parameters to `estimate=0` in Section 4.1
- **CRITICAL**: Set `nominalValue` to true value for fixed params
- Run Section 5.1.1 to verify fixed parameter configuration
- Proceed with optimization

**Workflow 3: Limited data (1 observable)**
- Estimate ≤ 2-3 model parameters
- Fix noise parameter (sigma) if known
- Use Section 7.3 recommendations

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| **NLLH huge at true parameters** | Run Section 5.1.1 - fixed parameters have wrong nominal values |
| **Sigma massively overestimated** | Run Section 7.1 - identifiability issue, need to fix parameters |
| **FIM has negative eigenvalues** | Run Section 7.1 and 7.3 - estimating too many parameters |
| **Poor parameter recovery** | Increase trials in Section 5.3 OR fix more parameters |
| **Import errors** | Check that all files in `python/` folder exist and are accessible |
| **BioNetGen not found** | Install BioNetGen and ensure it's in your PATH |

---

## External Dependencies

All notebooks require the following Python packages:
- `numpy` - Numerical computations
- `pandas` - Data manipulation
- `matplotlib` - Plotting
- `scipy` - Optimization and linear algebra
- `jupyter` - Notebook environment

**BioNetGen**: Must be installed and accessible via command line
- Download: [bionetgen.org](https://bionetgen.org)
- After installation, ensure `BNG2.pl` is in your system PATH

---

## File Structure

```
notebooks/
├── README.md                              # This file
├── python/                                # Python modules
│   ├── bngl_simulator.py                  # Core simulation
│   ├── parameter_estimator.py             # Core estimation
│   ├── debug_fixed_params.py              # Diagnostic
│   ├── identifiability_diagnosis.py       # Diagnostic
│   ├── fim_analysis_fixed.py              # Diagnostic
│   ├── diagnose_nllh.py                   # Experimental
│   ├── diagnose_sigma_issue.py            # Experimental
│   └── identifiability_helper_cell.py     # Helper
├── parameter_estimation_demo.ipynb        # Main tutorial
├── bngl_sensitivity_analysis.ipynb        # Sensitivity analysis
├── run_bngl.ipynb                         # Basic simulation
├── RasModelAnalyzer.ipynb                 # Ras model analysis
├── fceri_starter.bngl                     # FcεRI model (BNGL)
├── fceri_starter_sbml.xml                 # FcεRI model (SBML)
├── RasModel.bngl                          # Ras model
├── PARAMETER_FIXING_GUIDE.md              # Documentation
├── IDENTIFIABILITY_HELP_README.md         # Documentation
├── petab_files/                           # PETAb outputs
├── estimation_results/                    # Estimation outputs
├── results_fceri_starter/                 # Legacy outputs
└── results_RasModel/                      # Ras analysis outputs
```

---

## Citation

If you use these notebooks in your research, please cite:

- **BioNetGen**: Harris et al., Bioinformatics 32(21):3366-3368, 2016
- **PETAb**: Schmiester et al., PLoS Comput Biol 17(1):e1008646, 2021

---

## License

These notebooks are provided for educational and research purposes.

---

## Contact

For issues or questions about these tutorials, please contact the repository maintainer or open an issue on GitHub.

---

**Last updated**: 2026-07-08
