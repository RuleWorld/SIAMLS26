"""
Add this cell to your notebook to diagnose identifiability issues and
get recommendations on which parameters to fix.
"""

# Cell to add to notebook:
NOTEBOOK_CELL = """
# =============================================================================
# DIAGNOSE IDENTIFIABILITY ISSUES
# Run this to get recommendations on which parameters to fix
# =============================================================================

from identifiability_diagnosis import (
    analyze_parameter_sensitivity,
    compute_parameter_correlations,
    recommend_parameters_to_fix,
    visualize_identifiability_analysis
)

# 1. Analyze parameter sensitivity
# This tells you which parameters actually affect your measured observables
sensitivity_df = analyze_parameter_sensitivity(
    simulator=simulator,
    param_names=[p for p in estimator.estimable_params if not p.startswith('sigma_')],
    observables=measured_observables,
    conditions=experimental_conditions,
    true_params=true_params,
    t_end=max(measurement_times),
    n_steps=200
)

print("\\nSensitivity Analysis Complete!")
print("\\nAverage sensitivity by parameter:")
print(sensitivity_df.groupby('parameter')['avg_sensitivity'].mean().sort_values(ascending=False))

# 2. Estimate parameter correlations
# This tells you which parameters compensate for each other
model_params = [p for p in estimator.estimable_params if not p.startswith('sigma_')]
correlation_df = compute_parameter_correlations(
    estimator=estimator,
    param_vector=best_result['param_vector'],
    param_names=model_params,
    n_samples=200,
    perturbation=0.05
)

print("\\nParameter Correlations (log scale):")
print(correlation_df.round(3))

# 3. Get recommendation on which parameters to estimate vs fix
# You can adjust n_to_estimate based on your preference
recommendation = recommend_parameters_to_fix(
    sensitivity_df=sensitivity_df,
    correlation_df=correlation_df,
    param_names=model_params,
    n_to_estimate=3  # Change this number based on how many you want to estimate
)

# 4. Visualize everything
visualize_identifiability_analysis(
    sensitivity_df=sensitivity_df,
    correlation_df=correlation_df,
    recommendation=recommendation,
    save_path=f'{results_dir}/identifiability_diagnosis_{timestamp}.png'
)

# =============================================================================
# APPLY RECOMMENDATION
# Update parameters_to_estimate based on recommendation
# =============================================================================

print("\\n" + "="*80)
print("TO APPLY THIS RECOMMENDATION:")
print("="*80)
print("\\n1. Go back to Section 4.1 (PARAMETER ESTIMATION CONFIGURATION)")
print("\\n2. For parameters to FIX, change 'estimate': 1 to 'estimate': 0")
print("   and set 'nominalValue' to the true value")
print("\\nExample:")
print("    {")
print("        'parameterId': 'kmL',")
print("        'parameterName': 'Lyn unbinding rate',")
print("        'nominalValue': 20,          # Set to true/known value")
print("        'lowerBound': 1.0,")
print("        'upperBound': 100,")
print("        'parameterScale': 'log10',")
print("        'estimate': 0,               # Changed from 1 to 0")
print("    },")
print("\\n3. Re-run the optimization (Section 5)")
print("\\n4. Re-run the FIM analysis (Section 7)")
print("\\nExpected improvement:")
print("  - FIM should be positive definite (no negative eigenvalues)")
print("  - Smaller condition number")
print("  - Better parameter recovery")
"""

if __name__ == "__main__":
    print("Copy this cell into your Jupyter notebook after Section 6 (before FIM analysis):")
    print("\n" + "="*80)
    print(NOTEBOOK_CELL)
    print("="*80)
