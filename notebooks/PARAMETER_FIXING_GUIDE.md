# Guide: Choosing Which Parameters to Fix

## Why You Have Identifiability Issues

Your current setup:
- **1 observable** (BetaP2 only)
- **6 model parameters** (kp1, kp2, kpL, kmL, pLb, dm)
- **3 conditions** (different ligand concentrations)

This is a **parameter identifiability** problem: you're trying to estimate 6 parameters from limited data.

## Common Causes of Non-Identifiability

### 1. **Parameter Compensation (Correlation)**
Two parameters affect the output in similar ways, so changing one can be compensated by changing the other.

**Example from your model:**
- `kpL` and `kmL` define the Lyn binding equilibrium: `K_eq = kpL/kmL`
- Your data might only constrain the *ratio* `kpL/kmL`, not the individual rates
- Fixing one allows you to estimate the other

### 2. **Low Sensitivity**
Some parameters barely affect the observable you're measuring.

**Example:**
- If BetaP2 is insensitive to `kp1`, you can't estimate `kp1` from BetaP2 data alone
- These parameters should be fixed or measured using different observables

### 3. **Insufficient Data**
With only 1 observable, you can typically estimate ~3-4 parameters reliably.

## Decision Framework

### Step 1: Run the Diagnostic

Add this to your notebook (I've created the code in `identifiability_diagnosis.py`):

```python
from identifiability_diagnosis import (
    analyze_parameter_sensitivity,
    compute_parameter_correlations,
    recommend_parameters_to_fix,
    visualize_identifiability_analysis
)

# Analyze which parameters affect BetaP2
sensitivity_df = analyze_parameter_sensitivity(...)

# Find highly correlated parameter pairs
correlation_df = compute_parameter_correlations(...)

# Get automated recommendation
recommendation = recommend_parameters_to_fix(
    sensitivity_df, correlation_df,
    param_names=model_params,
    n_to_estimate=3  # Adjust this!
)
```

### Step 2: Interpret the Results

The diagnostic will tell you:

1. **Sensitivity ranking**: Which parameters BetaP2 is most sensitive to
   - High sensitivity (>0.5): Strongly affects observable, good to estimate
   - Low sensitivity (<0.1): Barely affects observable, good to fix

2. **Correlation pairs**: Which parameters compensate for each other
   - High correlation (|r| > 0.8): Fix one of the pair
   - Low correlation (|r| < 0.5): Can estimate both

3. **Recommendation**: Automated suggestion of which to estimate/fix

### Step 3: Make Your Decision

**Option A: Follow the automated recommendation**
- Safest for beginners
- Based on sensitivity and correlation analysis
- Usually picks the most identifiable subset

**Option B: Use domain knowledge**

Consider:
- **Measurable independently?** Fix it to the measured value
  - Example: Total concentrations (Lig_tot, Rec_tot) can be measured

- **Known from literature?** Fix it to the literature value
  - Example: If km1 and km2 are known to be ~0 for this system

- **Most biologically important?** Keep it estimable
  - Example: If pLb (phosphorylation rate) is your key question, estimate it

### Step 4: How Many Parameters to Estimate?

**Rule of thumb:** `n_params ≤ n_observables × n_conditions / 2`

For you:
- 1 observable × 3 conditions / 2 = **1-2 parameters** (conservative)
- But with good coverage, you can push to **3-4 parameters**

**Start conservative (3 parameters), then add more if FIM looks good.**

## Practical Strategies for Your Model

### Strategy 1: Fix Equilibrium Constants (Recommended)

The parameters come in pairs that define equilibria:

| Pair | Defines | Fix Which? |
|------|---------|------------|
| kpL/kmL | Lyn binding equilibrium | Fix kmL, estimate kpL |
| pLb/dm | Phosphorylation equilibrium | Estimate both (they're central!) |
| kp1/kp2 | Receptor binding/aggregation | Fix kp1, estimate kp2 |

**Rationale:**
- Fixing one of each pair resolves the compensation
- Keep the *forward* rates (kpL, pLb, kp2) as they're often more interpretable
- Fix the *reverse* rates (kmL, km1/km2=0) as they may be known

**Result: Estimate 3 parameters** (kpL, kp2, pLb) or (kp2, pLb, dm)

### Strategy 2: Fix Low-Sensitivity Parameters

Run the sensitivity analysis first, then fix the bottom 3 parameters.

**Likely candidates to fix (need to verify with analysis):**
- `kp1`: Initial binding might not affect steady-state BetaP2 much
- `kmL`: Reverse rate for Lyn binding
- `kpL`: If Lyn is saturated, its binding rate doesn't matter

**Result: Estimate 3 parameters** (the most sensitive ones)

### Strategy 3: Add More Observables (Best Long-Term)

Instead of fixing parameters, measure more observables!

**Easy wins:**
- Add `RecDimer2` (receptor dimerization) → helps identify kp1, kp2
- Add `LynBoundU` (Lyn binding) → helps identify kpL, kmL

With 3 observables × 3 conditions, you could estimate all 6 parameters.

## Implementation

### To Fix a Parameter:

In Section 4.1 of the notebook, change:

```python
{
    'parameterId': 'kmL',
    'parameterName': 'Lyn unbinding rate',
    'nominalValue': 20,              # ← Set to true/known value
    'lowerBound': 1.0,
    'upperBound': 100,
    'parameterScale': 'log10',
    'estimate': 0,                    # ← Change from 1 to 0
},
```

### Workflow:

1. **Run diagnostic** (using `identifiability_diagnosis.py`)
2. **Choose strategy** (fix 3 parameters based on sensitivity/correlation)
3. **Update configuration** (Section 4.1, set `estimate: 0` for fixed params)
4. **Re-run optimization** (Section 5)
5. **Check FIM** (Section 7 - should now be well-conditioned)
6. **Validate**:
   - All FIM eigenvalues > 0
   - Condition number < 1000
   - Parameter uncertainties < 50%

## Expected Outcomes

### Before Fixing (Current State)
- FIM: Negative eigenvalues, singular matrix
- Uncertainties: NaN or inf
- Condition number: inf
- Parameter recovery: Poor (17-163% error)

### After Fixing 3 Parameters
- FIM: All positive eigenvalues ✓
- Uncertainties: 10-30% for estimated parameters ✓
- Condition number: 10-100 ✓
- Parameter recovery: Good (5-20% error) ✓

## Quick Start

Just want to get started? Try this:

**Fix these 3 (likely low sensitivity or correlated):**
- `kmL = 20` (reverse rate)
- `kp1 = 1.328e-7` (initial binding)
- `kpL = 0.05` (or estimate this one if you care about Lyn)

**Estimate these 3 (likely high impact on BetaP2):**
- `kp2` (receptor aggregation)
- `pLb` (phosphorylation - the star of the show!)
- `dm` (dephosphorylation - balances pLb)

Then run the diagnostic to see if this was a good choice!
