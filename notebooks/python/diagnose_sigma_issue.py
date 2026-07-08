"""
Diagnose why sigma is being estimated so poorly.

The noise parameter should be well-constrained by the residuals,
but we're seeing sigma_BetaP2 = 43.6 when true value is 0.323.

This usually indicates:
1. Poor model fit (large residuals)
2. Model structure error (wrong parameters → wrong predictions)
3. Optimizer is using sigma to hide model inadequacy
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict


def diagnose_sigma_estimation(
    estimator,
    best_result: Dict,
    true_param_dict: Dict,
    measurements_df: pd.DataFrame,
    results_dir: str,
    timestamp: str
):
    """
    Diagnose why sigma is poorly estimated.

    Computes:
    1. Actual residuals at optimum
    2. What sigma SHOULD be (RMS of residuals)
    3. What sigma IS estimated as
    4. Residuals at TRUE parameters (to separate noise from model error)
    """

    print("="*80)
    print("DIAGNOSING SIGMA ESTIMATION")
    print("="*80)

    # 1. Get parameter vectors
    estimated_vector = best_result['param_vector']
    estimated_params = estimator.vector_to_params(estimated_vector)

    true_vector = estimator.params_to_vector(true_param_dict)

    # 2. Simulate with ESTIMATED parameters
    print("\n1. Simulating with ESTIMATED parameters...")
    results_estimated = estimator.simulate_conditions(estimated_vector, n_steps=200)

    # 3. Simulate with TRUE parameters
    print("2. Simulating with TRUE parameters...")
    results_true = estimator.simulate_conditions(true_vector, n_steps=200)

    # 4. Compute residuals for each condition
    residuals_estimated = []
    residuals_true = []
    measurements = []
    times = []
    conditions = []

    for _, cond_row in estimator.conditions_df.iterrows():
        cond_id = cond_row['conditionId']

        result_est = results_estimated[cond_id]
        result_true = results_true[cond_id]

        # Get measurements for this condition
        cond_data = measurements_df[
            measurements_df['simulationConditionId'] == cond_id
        ]

        for _, meas_row in cond_data.iterrows():
            obs_id = meas_row['observableId']
            time = meas_row['time']
            measurement = meas_row['measurement']

            # Find time index
            time_idx = np.argmin(np.abs(result_est['time'] - time))

            # Get observable formula and evaluate
            obs_formula = estimator.observables_df[
                estimator.observables_df['observableId'] == obs_id
            ]['observableFormula'].iloc[0]

            # For simple case where formula = observable name
            if obs_formula == obs_id:
                pred_est = result_est[obs_id][time_idx]
                pred_true = result_true[obs_id][time_idx]
            else:
                # Would need to evaluate formula, but for now assume direct
                pred_est = result_est[obs_id][time_idx]
                pred_true = result_true[obs_id][time_idx]

            residuals_estimated.append(measurement - pred_est)
            residuals_true.append(measurement - pred_true)
            measurements.append(measurement)
            times.append(time)
            conditions.append(cond_id)

    residuals_estimated = np.array(residuals_estimated)
    residuals_true = np.array(residuals_true)
    measurements = np.array(measurements)

    # 5. Compute statistics
    print("\n" + "="*80)
    print("RESIDUAL STATISTICS")
    print("="*80)

    # RMS of residuals
    rms_estimated = np.sqrt(np.mean(residuals_estimated**2))
    rms_true = np.sqrt(np.mean(residuals_true**2))

    # Mean absolute residual
    mae_estimated = np.mean(np.abs(residuals_estimated))
    mae_true = np.mean(np.abs(residuals_true))

    # What sigma should be (RMS of residuals)
    sigma_from_residuals_est = rms_estimated
    sigma_from_residuals_true = rms_true

    # What sigma was estimated as
    sigma_estimated = estimated_params.get('sigma_BetaP2', np.nan)
    sigma_true = true_param_dict.get('sigma_BetaP2', np.nan)

    print(f"\nWith ESTIMATED parameters:")
    print(f"  RMS residual:        {rms_estimated:.4f}")
    print(f"  Mean |residual|:     {mae_estimated:.4f}")
    print(f"  Max |residual|:      {np.max(np.abs(residuals_estimated)):.4f}")
    print(f"  → Sigma should be:   {sigma_from_residuals_est:.4f}")
    print(f"  → Sigma estimated:   {sigma_estimated:.4f}")
    print(f"  → Ratio:             {sigma_estimated / sigma_from_residuals_est:.2f}x")

    print(f"\nWith TRUE parameters:")
    print(f"  RMS residual:        {rms_true:.4f}")
    print(f"  Mean |residual|:     {mae_true:.4f}")
    print(f"  Max |residual|:      {np.max(np.abs(residuals_true)):.4f}")
    print(f"  → Sigma should be:   {sigma_from_residuals_true:.4f}")
    print(f"  → True sigma used:   {sigma_true:.4f}")
    print(f"  → Ratio:             {sigma_true / sigma_from_residuals_true:.2f}x")

    # 6. Decompose variance
    print("\n" + "="*80)
    print("VARIANCE DECOMPOSITION")
    print("="*80)

    # Total variance in residuals (estimated params)
    total_var = np.var(residuals_estimated)

    # Variance due to model error (true params give us "pure" noise)
    noise_var = np.var(residuals_true)

    # Variance due to model inadequacy
    model_error_var = total_var - noise_var

    print(f"\nVariance in residuals (estimated params): {total_var:.4f}")
    print(f"  Due to measurement noise:               {noise_var:.4f} ({100*noise_var/total_var:.1f}%)")
    print(f"  Due to model inadequacy:                {model_error_var:.4f} ({100*model_error_var/total_var:.1f}%)")

    print(f"\nConclusion:")
    if model_error_var > noise_var:
        print(f"  ⚠ Model error dominates! The poor parameter estimates cause")
        print(f"    large residuals, and sigma is inflated to compensate.")
        print(f"    Sigma is absorbing MODEL ERROR, not just measurement noise.")
    else:
        print(f"  ✓ Noise dominates, sigma should be well-estimated.")

    # 7. Visualize
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Residuals over time (estimated params)
    ax = axes[0, 0]
    for cond in np.unique(conditions):
        mask = np.array(conditions) == cond
        ax.scatter(np.array(times)[mask], residuals_estimated[mask],
                  label=cond, alpha=0.7, s=50)
    ax.axhline(0, color='k', linestyle='--', linewidth=1)
    ax.axhline(sigma_estimated, color='r', linestyle='--', linewidth=2,
              label=f'Estimated σ = {sigma_estimated:.2f}')
    ax.axhline(-sigma_estimated, color='r', linestyle='--', linewidth=2)
    ax.set_xlabel('Time (s)', fontweight='bold')
    ax.set_ylabel('Residual (data - model)', fontweight='bold')
    ax.set_title('Residuals with ESTIMATED Parameters', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Residuals over time (true params)
    ax = axes[0, 1]
    for cond in np.unique(conditions):
        mask = np.array(conditions) == cond
        ax.scatter(np.array(times)[mask], residuals_true[mask],
                  label=cond, alpha=0.7, s=50)
    ax.axhline(0, color='k', linestyle='--', linewidth=1)
    ax.axhline(sigma_true, color='g', linestyle='--', linewidth=2,
              label=f'True σ = {sigma_true:.2f}')
    ax.axhline(-sigma_true, color='g', linestyle='--', linewidth=2)
    ax.set_xlabel('Time (s)', fontweight='bold')
    ax.set_ylabel('Residual (data - model)', fontweight='bold')
    ax.set_title('Residuals with TRUE Parameters', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Residual distribution (estimated)
    ax = axes[1, 0]
    ax.hist(residuals_estimated, bins=15, alpha=0.7, color='blue', edgecolor='black')
    ax.axvline(0, color='k', linestyle='--', linewidth=1)
    ax.axvline(rms_estimated, color='r', linestyle='--', linewidth=2,
              label=f'RMS = {rms_estimated:.2f}')
    ax.axvline(-rms_estimated, color='r', linestyle='--', linewidth=2)
    ax.set_xlabel('Residual', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title(f'Residual Distribution (Estimated)\nσ_est={sigma_estimated:.2f}, RMS={rms_estimated:.2f}',
                fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: Residual distribution (true)
    ax = axes[1, 1]
    ax.hist(residuals_true, bins=15, alpha=0.7, color='green', edgecolor='black')
    ax.axvline(0, color='k', linestyle='--', linewidth=1)
    ax.axvline(rms_true, color='g', linestyle='--', linewidth=2,
              label=f'RMS = {rms_true:.2f}')
    ax.axvline(-rms_true, color='g', linestyle='--', linewidth=2)
    ax.set_xlabel('Residual', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title(f'Residual Distribution (True)\nσ_true={sigma_true:.2f}, RMS={rms_true:.2f}',
                fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_path = f'{results_dir}/sigma_diagnosis_{timestamp}.png'
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved figure: {save_path}")
    plt.show()

    # 8. Recommendations
    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)

    if sigma_estimated > 2 * sigma_from_residuals_est:
        print("\n⚠ Sigma is significantly overestimated!")
        print("\nPossible causes:")
        print("1. Poor parameter fit → large residuals → sigma absorbs model error")
        print("2. Optimizer got stuck in local minimum")
        print("3. Some parameters are non-identifiable (causing poor fit)")

        print("\nSolutions:")
        print("1. Fix some parameters (use Section 7.1 recommendations)")
        print("2. Run more optimization iterations or trials")
        print("3. Try fixing sigma to true value and re-optimize model parameters")
        print("4. Add more observables to better constrain the model")

    return {
        'rms_estimated': rms_estimated,
        'rms_true': rms_true,
        'sigma_estimated': sigma_estimated,
        'sigma_true': sigma_true,
        'residuals_estimated': residuals_estimated,
        'residuals_true': residuals_true,
        'noise_fraction': noise_var / total_var if total_var > 0 else 0,
        'model_error_fraction': model_error_var / total_var if total_var > 0 else 0,
    }
