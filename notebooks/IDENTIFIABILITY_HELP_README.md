# Identifiability Analysis - Quick Start

## What I've Added to Help You

### 1. New Notebook Section: 7.1 Diagnose Identifiability

Run this **before** the FIM analysis (now in Section 7.2) to get:

**Automated analysis showing:**
- Which parameters BetaP2 is sensitive to
- Which parameters are highly correlated (compensate for each other)
- **Specific recommendations** on which 3 parameters to estimate and which 3 to fix

**Output visualization:**
- Sensitivity heatmap
- Correlation matrix
- Sensitivity ranking
- Recommendation summary

### 2. Helper Modules Created

**`identifiability_diagnosis.py`** - Main diagnostic tool with 4 functions:
- `analyze_parameter_sensitivity()` - Computes how much each parameter affects observables
- `compute_parameter_correlations()` - Estimates parameter correlations via Monte Carlo
- `recommend_parameters_to_fix()` - Automated decision algorithm
- `visualize_identifiability_analysis()` - Creates comprehensive visualization

**`PARAMETER_FIXING_GUIDE.md`** - Detailed guide explaining:
- Why you have identifiability issues
- How to interpret the diagnostic results
- Three practical strategies for your specific model
- Step-by-step implementation instructions

### 3. Fixed FIM Analysis

The FIM computation (now Section 7.2) has been corrected to:
- Exclude noise parameters (prevents singularity)
- Properly handle log10 → natural space transformations
- Use adaptive step sizes for numerical stability

## How to Use

### Quick Workflow

1. **Run the notebook up to Section 7.1** (Diagnose Identifiability)

2. **Review the output** - Look for:
   ```
   FINAL RECOMMENDATION

   Estimate these 3 parameters:
     ✓ pLb    (sensitivity: 2.5)
     ✓ dm     (sensitivity: 1.8)
     ✓ kp2    (sensitivity: 0.7)

   Fix these 3 parameters:
     ✗ kmL    - correlated with kpL (r=0.92)
     ✗ kp1    - low sensitivity (0.05)
     ✗ kpL    - correlated with kmL (r=0.92)
   ```

3. **Apply the recommendation**:
   - Go to **Section 4.1** (Parameter Estimation Configuration)
   - For each parameter to FIX, change:
     ```python
     'estimate': 1   →   'estimate': 0
     ```
   - Set `nominalValue` to the true value

4. **Re-run from Section 4.3** onwards:
   - PETAb file creation
   - Optimization
   - FIM analysis

5. **Verify success** - Section 7.2 should now show:
   - ✓ All FIM eigenvalues > 0
   - ✓ Condition number < 1000
   - ✓ Reasonable parameter uncertainties (10-50%)

### Example: Fixing Parameters

**Before** (Section 4.1):
```python
{
    'parameterId': 'kmL',
    'parameterName': 'Lyn unbinding rate',
    'nominalValue': 10,
    'lowerBound': 1.0,
    'upperBound': 100,
    'parameterScale': 'log10',
    'estimate': 1,  # ← Estimating
},
```

**After** (following recommendation to fix kmL):
```python
{
    'parameterId': 'kmL',
    'parameterName': 'Lyn unbinding rate',
    'nominalValue': 20,  # ← Set to true value
    'lowerBound': 1.0,
    'upperBound': 100,
    'parameterScale': 'log10',
    'estimate': 0,  # ← Fixed (not estimating)
},
```

## Understanding the Diagnostic

### Sensitivity Analysis

**High sensitivity (>0.5):**
- Parameter strongly affects observable
- **Good candidate to estimate**
- Example: `pLb` (phosphorylation rate) → BetaP2 (phosphorylated receptors)

**Low sensitivity (<0.1):**
- Parameter barely affects observable
- **Good candidate to fix**
- Example: `kp1` might have low sensitivity if dimerization saturates quickly

### Correlation Analysis

**High correlation (|r| > 0.8):**
- Parameters compensate for each other
- Can only estimate their ratio or sum, not individually
- **Fix one of the pair**
- Example: `kpL`/`kmL` likely correlated (define equilibrium constant)

**Low correlation (|r| < 0.3):**
- Parameters affect output independently
- **Can estimate both**

### The Recommendation Algorithm

1. **Identify low-sensitivity parameters** → candidates to fix
2. **Identify correlated pairs** → fix one of each pair (the less sensitive one)
3. **Rank remaining by sensitivity** → estimate top N (you choose N)

**You can override** the automated recommendation if you have:
- Literature values for specific parameters
- Biological reasons to estimate certain parameters
- Additional measurements for some parameters

## Troubleshooting

### Problem: Diagnostic says to estimate parameter X, but I want to estimate Y instead

**Solution:** Manually override the recommendation
- The diagnostic is a suggestion, not a requirement
- If you have domain knowledge, use it!
- Just make sure total estimated ≤ 3-4 for 1 observable

### Problem: Even after fixing parameters, FIM is still poorly conditioned

**Possible causes:**
1. Not enough parameters fixed (try fixing 4, estimate only 2)
2. Wrong parameters fixed (try different combination)
3. Insufficient data (need more timepoints or observables)

**Solutions:**
- Try the "Quick Start" suggestion in `PARAMETER_FIXING_GUIDE.md`
- Add more observables (RecDimer2, LynBoundU)
- Use profile likelihood instead of FIM

### Problem: Diagnostic shows all parameters have low correlation

**This is good!** But if FIM is still singular:
- Parameters might have low *individual* correlation
- But high *collective* non-identifiability
- Solution: Still need to fix some based on sensitivity alone

## Expected Results

### With All 6 Parameters (Current)
- Fit errors: 17-163%
- FIM: Negative eigenvalues ❌
- Condition number: inf ❌
- Uncertainties: NaN ❌

### After Fixing 3 Parameters (Recommended)
- Fit errors: 5-20%
- FIM: All positive eigenvalues ✓
- Condition number: 10-100 ✓
- Uncertainties: 10-30% ✓

## Files Reference

| File | Purpose |
|------|---------|
| `identifiability_diagnosis.py` | Diagnostic functions |
| `PARAMETER_FIXING_GUIDE.md` | Detailed explanation and strategies |
| `FIM_ANALYSIS_FIX_README.md` | Technical details on FIM fixes |
| `fim_analysis_fixed.py` | Fixed FIM computation |

## Next Steps

1. ✅ Reload the notebook (to see new Section 7.1)
2. ✅ Run through Section 7.1 to get recommendations
3. ✅ Apply recommendations in Section 4.1
4. ✅ Re-run optimization and FIM
5. ✅ Enjoy well-conditioned parameter estimation!

Questions? Check `PARAMETER_FIXING_GUIDE.md` for detailed explanations.
