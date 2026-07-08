"""
Diagnose why NLLH behaves strangely.

Key questions:
1. Why is NLLH negative? (Answer: log(sigma) term dominates when sigma is small)
2. Why is NLLH at optimum better than at true parameters?
3. Is the optimizer gaming the likelihood by manipulating sigma?
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Tuple


def decompose_nllh(
    estimator,
    param_vector: np.ndarray,
    verbose: bool = True
) -> Dict:
    """
    Decompose NLLH into its components to understand what's driving it.

    NLLH = Σ[ 0.5*(y-ŷ)²/σ² + log(σ) + 0.5*log(2π) ]
           \_____________/   \______/   \__________/
           residual term   sigma term   constant

    Returns breakdown of each component.
    """
    param_dict = estimator.vector_to_params(param_vector)

    # Track contributions
    total_nllh = 0.0
    residual_term_sum = 0.0
    sigma_term_sum = 0.0
    constant_term_sum = 0.0
    n_points = 0

    residuals = []
    predictions = []
    measurements = []
    sigma_values = []

    # Loop through conditions
    for _, cond_row in estimator.conditions_df.iterrows():
        cond_id = cond_row['conditionId']

        # Get condition parameters
        cond_params = {}
        for col in cond_row.index:
            if col != 'conditionId' and col != 'conditionName' and pd.notna(cond_row[col]):
                cond_params[col] = cond_row[col]

        # Get measurements
        cond_measurements = estimator.measurements_df[
            estimator.measurements_df['simulationConditionId'] == cond_id
        ]

        if len(cond_measurements) == 0:
            continue

        # Simulate
        meas_times = np.sort(cond_measurements['time'].unique())
        t_end = float(meas_times[-1])
        n_steps = int(max(int(t_end * 2), 100))

        simulator = estimator._get_simulator(param_dict)
        result = simulator.simulate(t_end=t_end, n_steps=n_steps,
                                   reference_values=cond_params, reset=True)

        # Compute likelihood for each measurement
        for _, meas_row in cond_measurements.iterrows():
            obs_id = meas_row['observableId']
            time = meas_row['time']
            measurement = meas_row['measurement']

            time_idx = np.argmin(np.abs(result['time'] - time))

            # Get observable formula
            obs_row = estimator.observables_df[
                estimator.observables_df['observableId'] == obs_id
            ]
            obs_formula = obs_row.iloc[0]['observableFormula']
            noise_formula = obs_row.iloc[0]['noiseFormula']

            # Evaluate prediction
            prediction = estimator._eval_observable_formula(
                obs_formula, result, param_dict, time_idx
            )

            # Get noise
            obs_params_str = meas_row.get('observableParameters', '')
            obs_params = estimator._parse_observable_parameters(obs_params_str)

            if 'SEM' in obs_params:
                noise_std = obs_params['SEM']
            else:
                noise_std = estimator._eval_observable_formula(
                    noise_formula, result, param_dict, time_idx
                )

            # Compute NLLH components
            residual = measurement - prediction

            residual_term = 0.5 * (residual / noise_std) ** 2
            sigma_term = np.log(noise_std)
            constant_term = 0.5 * np.log(2 * np.pi)

            residual_term_sum += residual_term
            sigma_term_sum += sigma_term
            constant_term_sum += constant_term
            total_nllh += residual_term + sigma_term + constant_term

            residuals.append(residual)
            predictions.append(prediction)
            measurements.append(measurement)
            sigma_values.append(noise_std)
            n_points += 1

    result = {
        'total_nllh': total_nllh,
        'residual_term': residual_term_sum,
        'sigma_term': sigma_term_sum,
        'constant_term': constant_term_sum,
        'n_points': n_points,
        'residuals': np.array(residuals),
        'predictions': np.array(predictions),
        'measurements': np.array(measurements),
        'sigma_values': np.array(sigma_values),
        'rms_residual': np.sqrt(np.mean(np.array(residuals)**2)),
        'mean_sigma': np.mean(sigma_values),
    }

    if verbose:
        print("NLLH Decomposition:")
        print(f"  Total NLLH:        {result['total_nllh']:10.4f}")
        print(f"  Residual term:     {result['residual_term']:10.4f}  (how well model fits)")
        print(f"  Sigma term:        {result['sigma_term']:10.4f}  (penalty for large σ)")
        print(f"  Constant term:     {result['constant_term']:10.4f}  (π normalization)")
        print(f"\n  N datapoints:      {result['n_points']}")
        print(f"  RMS residual:      {result['rms_residual']:10.4f}")
        print(f"  Mean sigma:        {result['mean_sigma']:10.4f}")
        print(f"  Ratio (RMS/σ):     {result['rms_residual']/result['mean_sigma']:10.4f}")

    return result


def compare_parameter_sets(
    estimator,
    param_vectors: Dict[str, np.ndarray],
    param_dicts: Dict[str, Dict],
    results_dir: str,
    timestamp: str
):
    """
    Compare NLLH and fit quality for different parameter sets.

    Args:
        param_vectors: Dict of {name: param_vector}
        param_dicts: Dict of {name: param_dict} (for display)
    """
    print("="*80)
    print("COMPARING PARAMETER SETS")
    print("="*80)

    results = {}

    for name, vector in param_vectors.items():
        print(f"\n{name}:")
        print("-"*40)

        # Show parameter values
        pdict = param_dicts[name]
        model_params = {k: v for k, v in pdict.items() if not k.startswith('sigma_')}
        print("Model parameters:")
        for k, v in model_params.items():
            print(f"  {k:15s} = {v:.6e}")
        print(f"Noise parameter:")
        for k, v in pdict.items():
            if k.startswith('sigma_'):
                print(f"  {k:15s} = {v:.6e}")

        # Decompose NLLH
        result = decompose_nllh(estimator, vector, verbose=True)
        results[name] = result

    # Compare
    print("\n" + "="*80)
    print("COMPARISON SUMMARY")
    print("="*80)

    print(f"\n{'Parameter Set':<20} {'NLLH':<12} {'Residual':<12} {'Sigma':<12} {'RMS/σ':<12}")
    print("-"*70)

    for name, res in results.items():
        print(f"{name:<20} {res['total_nllh']:<12.4f} {res['residual_term']:<12.4f} "
              f"{res['sigma_term']:<12.4f} {res['rms_residual']/res['mean_sigma']:<12.4f}")

    # Explanation
    print("\n" + "="*80)
    print("INTERPRETATION")
    print("="*80)

    print("\nNLLH = Residual term + Sigma term + Constant")
    print("\nResidual term: (y-ŷ)²/σ²  [want low = good fit]")
    print("Sigma term:    log(σ)     [negative for small σ, positive for large σ]")
    print("\nKey insights:")
    print("1. Negative NLLH is OK! Just means log(σ) dominates (small σ)")
    print("2. Lower NLLH = better (we minimize)")
    print("3. But beware: optimizer can game NLLH by making σ very small")
    print("   - This inflates residual term but sigma term becomes very negative")
    print("   - Check RMS/σ ratio: should be ~1 for good fit")
    print("   - If RMS/σ >> 1: residuals too large (poor fit)")
    print("   - If RMS/σ << 1: sigma too small (overfitting)")

    # Find best and worst
    best_name = min(results.keys(), key=lambda k: results[k]['total_nllh'])
    best = results[best_name]

    print(f"\nBest NLLH: {best_name} = {best['total_nllh']:.4f}")

    # Check if it's gaming sigma
    rms_sigma_ratio = best['rms_residual'] / best['mean_sigma']
    if rms_sigma_ratio < 0.5:
        print(f"\n⚠ WARNING: RMS/σ = {rms_sigma_ratio:.2f} << 1")
        print("  Sigma is too small! Optimizer may be gaming the likelihood.")
        print("  The fit might look good but sigma is unrealistic.")
    elif rms_sigma_ratio > 2:
        print(f"\n⚠ WARNING: RMS/σ = {rms_sigma_ratio:.2f} >> 1")
        print("  Residuals are large compared to sigma.")
        print("  Either sigma is too small OR model fit is poor.")
    else:
        print(f"\n✓ RMS/σ = {rms_sigma_ratio:.2f} is reasonable (0.5-2 range)")
        print("  Sigma appears consistent with actual residuals.")

    # Visualize
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: NLLH comparison
    ax = axes[0, 0]
    names = list(results.keys())
    nllhs = [results[n]['total_nllh'] for n in names]
    colors = ['green' if n == best_name else 'gray' for n in names]
    ax.bar(names, nllhs, color=colors, alpha=0.7)
    ax.set_ylabel('NLLH', fontweight='bold')
    ax.set_title('NLLH Comparison (Lower = Better)', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 2: NLLH decomposition
    ax = axes[0, 1]
    residual_terms = [results[n]['residual_term'] for n in names]
    sigma_terms = [results[n]['sigma_term'] for n in names]

    x = np.arange(len(names))
    width = 0.35
    ax.bar(x - width/2, residual_terms, width, label='Residual term', alpha=0.7)
    ax.bar(x + width/2, sigma_terms, width, label='Sigma term', alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('Contribution to NLLH', fontweight='bold')
    ax.set_title('NLLH Decomposition', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 3: Residuals distribution
    ax = axes[1, 0]
    for name, res in results.items():
        ax.hist(res['residuals'], bins=10, alpha=0.5, label=name)
    ax.set_xlabel('Residual', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Residual Distributions', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: RMS/sigma ratio
    ax = axes[1, 1]
    ratios = [results[n]['rms_residual']/results[n]['mean_sigma'] for n in names]
    colors = ['green' if 0.5 < r < 2 else 'orange' if 0.3 < r < 3 else 'red'
              for r in ratios]
    ax.bar(names, ratios, color=colors, alpha=0.7)
    ax.axhline(1, color='k', linestyle='--', linewidth=2, label='Ideal = 1')
    ax.axhline(0.5, color='orange', linestyle='--', linewidth=1, alpha=0.5)
    ax.axhline(2, color='orange', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_ylabel('RMS residual / σ', fontweight='bold')
    ax.set_title('Fit Quality (should be ~1)', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    save_path = f'{results_dir}/nllh_diagnosis_{timestamp}.png'
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved figure: {save_path}")
    plt.show()

    return results
