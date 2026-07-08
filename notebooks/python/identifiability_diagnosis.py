"""
Diagnose identifiability issues and recommend which parameters to fix.

This module helps identify:
1. Which parameters are most/least sensitive for each observable
2. Which parameters are highly correlated
3. Which parameters should be fixed to resolve identifiability issues
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
from .bngl_simulator import BNGLSimulator
from .parameter_estimator import ParameterEstimator


def analyze_parameter_sensitivity(
    simulator: BNGLSimulator,
    param_names: List[str],
    observables: List[str],
    conditions: List[dict],
    true_params: dict,
    t_end: float = 1200,
    n_steps: int = 200
) -> pd.DataFrame:
    """
    Compute sensitivity of each observable to each parameter.

    Returns DataFrame with sensitivity metrics for each parameter-observable pair.
    """
    results = []

    print("Computing parameter sensitivities...")
    print(f"Parameters: {param_names}")
    print(f"Observables: {observables}")

    for obs in observables:
        print(f"\nAnalyzing sensitivity for {obs}:")

        for param in param_names:
            # Set all parameters to true values
            for p, v in true_params.items():
                if p in simulator.parameters:
                    simulator.set_parameter(p, v)

            # Get nominal value
            nominal_value = simulator.get_parameter(param)

            # Run simulations across all conditions
            sensitivities = []

            for condition in conditions:
                cond_id = condition['conditionId']
                cond_params = {k: v for k, v in condition.items()
                              if k not in ['conditionId', 'conditionName']}

                # Nominal simulation
                result_nom = simulator.simulate(
                    t_end=t_end, n_steps=n_steps,
                    reference_values=cond_params, reset=True
                )

                # Perturbed simulation (+10%)
                simulator.set_parameter(param, nominal_value * 1.1)
                result_pert = simulator.simulate(
                    t_end=t_end, n_steps=n_steps,
                    reference_values=cond_params, reset=True
                )

                # Compute normalized sensitivity at final time
                obs_nom = result_nom[obs][-1]
                obs_pert = result_pert[obs][-1]

                # Normalized sensitivity: (dO/O) / (dp/p)
                if obs_nom != 0:
                    norm_sens = ((obs_pert - obs_nom) / obs_nom) / 0.1
                else:
                    norm_sens = 0

                sensitivities.append(abs(norm_sens))

                # Reset parameter
                simulator.set_parameter(param, nominal_value)

            # Average sensitivity across conditions
            avg_sensitivity = np.mean(sensitivities)
            max_sensitivity = np.max(sensitivities)

            results.append({
                'parameter': param,
                'observable': obs,
                'avg_sensitivity': avg_sensitivity,
                'max_sensitivity': max_sensitivity,
            })

            print(f"  {param:15s}: avg={avg_sensitivity:.4f}, max={max_sensitivity:.4f}")

    return pd.DataFrame(results)


def compute_parameter_correlations(
    estimator: ParameterEstimator,
    param_vector: np.ndarray,
    param_names: List[str],
    n_samples: int = 100,
    perturbation: float = 0.05
) -> pd.DataFrame:
    """
    Estimate parameter correlations by sampling around the optimum.

    This is a Monte Carlo approach to understand parameter compensation.
    """
    print(f"\nComputing parameter correlations via Monte Carlo sampling...")
    print(f"Sampling {n_samples} points around optimum")

    # Get indices for parameters of interest (exclude noise params)
    param_indices = [i for i, name in enumerate(estimator.estimable_params)
                     if name in param_names]

    # Generate samples around the optimum
    samples = []
    nllh_values = []

    for i in range(n_samples):
        # Perturb parameters
        perturbations = np.random.normal(0, perturbation, size=len(param_indices))
        perturbed_vector = param_vector.copy()

        for idx, param_idx in enumerate(param_indices):
            perturbed_vector[param_idx] += perturbations[idx]

        # Evaluate likelihood
        nllh = estimator.compute_nllh(perturbed_vector)

        # Only keep if not too far from optimum (within 10 NLLH units)
        if nllh < estimator.compute_nllh(param_vector) + 10:
            # Convert to natural space
            param_dict = estimator.vector_to_params(perturbed_vector)
            sample = [param_dict[name] for name in param_names]
            samples.append(sample)
            nllh_values.append(nllh)

    print(f"Accepted {len(samples)}/{n_samples} samples")

    # Compute correlation matrix
    samples_array = np.array(samples)

    # Log-transform for better correlation estimation
    log_samples = np.log10(samples_array + 1e-10)

    corr_matrix = np.corrcoef(log_samples.T)

    # Create DataFrame
    corr_df = pd.DataFrame(corr_matrix, index=param_names, columns=param_names)

    return corr_df


def recommend_parameters_to_fix(
    sensitivity_df: pd.DataFrame,
    correlation_df: pd.DataFrame,
    param_names: List[str],
    n_to_estimate: int = 3
) -> Dict:
    """
    Recommend which parameters to fix based on sensitivity and correlation.

    Strategy:
    1. Identify low-sensitivity parameters (good candidates to fix)
    2. Identify highly correlated parameter pairs (fix one of each pair)
    3. Rank remaining parameters by sensitivity
    4. Recommend top n_to_estimate to keep estimating

    Returns:
        Dict with 'to_estimate', 'to_fix', and 'reasoning'
    """
    print("\n" + "="*80)
    print("PARAMETER SELECTION RECOMMENDATION")
    print("="*80)

    # 1. Identify low-sensitivity parameters
    avg_sens_by_param = sensitivity_df.groupby('parameter')['avg_sensitivity'].mean()
    sensitivity_threshold = 0.1

    low_sensitivity = avg_sens_by_param[avg_sens_by_param < sensitivity_threshold].index.tolist()

    print(f"\n1. Low sensitivity parameters (avg < {sensitivity_threshold}):")
    for param in low_sensitivity:
        print(f"   {param}: {avg_sens_by_param[param]:.4f}")

    # 2. Identify highly correlated pairs
    correlation_threshold = 0.85
    correlated_pairs = []

    for i, param1 in enumerate(param_names):
        for j, param2 in enumerate(param_names[i+1:], start=i+1):
            corr = abs(correlation_df.loc[param1, param2])
            if corr > correlation_threshold:
                correlated_pairs.append((param1, param2, corr))

    print(f"\n2. Highly correlated parameter pairs (|r| > {correlation_threshold}):")
    for param1, param2, corr in correlated_pairs:
        print(f"   {param1} <-> {param2}: r = {corr:.3f}")

    # 3. Build recommendation
    candidates_to_fix = set(low_sensitivity)

    # For each correlated pair, fix the one with lower sensitivity
    for param1, param2, corr in correlated_pairs:
        sens1 = avg_sens_by_param[param1]
        sens2 = avg_sens_by_param[param2]

        if sens1 < sens2:
            candidates_to_fix.add(param1)
            print(f"\n   Recommend fixing {param1} (sens={sens1:.4f}) over {param2} (sens={sens2:.4f})")
        else:
            candidates_to_fix.add(param2)
            print(f"\n   Recommend fixing {param2} (sens={sens2:.4f}) over {param1} (sens={sens1:.4f})")

    # 4. Select top n_to_estimate by sensitivity
    remaining_params = [p for p in param_names if p not in candidates_to_fix]

    if len(remaining_params) < n_to_estimate:
        # Need to reclaim some from candidates_to_fix
        # Sort all by sensitivity
        all_sorted = avg_sens_by_param.sort_values(ascending=False)
        to_estimate = all_sorted.index[:n_to_estimate].tolist()
        to_fix = [p for p in param_names if p not in to_estimate]
    else:
        # More remaining than needed, pick top n_to_estimate by sensitivity
        remaining_sens = avg_sens_by_param[remaining_params].sort_values(ascending=False)
        to_estimate = remaining_sens.index[:n_to_estimate].tolist()
        to_fix = [p for p in param_names if p not in to_estimate]

    # Generate reasoning
    reasoning = []
    for param in to_fix:
        reasons = []
        if param in low_sensitivity:
            reasons.append(f"low sensitivity ({avg_sens_by_param[param]:.4f})")

        for p1, p2, corr in correlated_pairs:
            if param == p1:
                reasons.append(f"correlated with {p2} (r={corr:.3f})")
            elif param == p2:
                reasons.append(f"correlated with {p1} (r={corr:.3f})")

        reasoning.append((param, ", ".join(reasons) if reasons else "lower sensitivity than alternatives"))

    print("\n" + "="*80)
    print("FINAL RECOMMENDATION")
    print("="*80)
    print(f"\nEstimate these {len(to_estimate)} parameters:")
    for param in to_estimate:
        print(f"   ✓ {param:15s} (sensitivity: {avg_sens_by_param[param]:.4f})")

    print(f"\nFix these {len(to_fix)} parameters:")
    for param, reason in reasoning:
        print(f"   ✗ {param:15s} - {reason}")

    return {
        'to_estimate': to_estimate,
        'to_fix': to_fix,
        'reasoning': reasoning,
        'sensitivity': avg_sens_by_param.to_dict(),
        'correlated_pairs': correlated_pairs
    }


def visualize_identifiability_analysis(
    sensitivity_df: pd.DataFrame,
    correlation_df: pd.DataFrame,
    recommendation: Dict,
    save_path: str = None
):
    """Create comprehensive visualization of identifiability analysis."""

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    param_names = list(correlation_df.index)

    # 1. Sensitivity heatmap
    ax1 = fig.add_subplot(gs[0, 0])

    # Pivot sensitivity data
    sens_pivot = sensitivity_df.pivot(
        index='parameter',
        columns='observable',
        values='avg_sensitivity'
    )

    # Reorder by total sensitivity
    total_sens = sens_pivot.sum(axis=1).sort_values(ascending=False)
    sens_pivot = sens_pivot.loc[total_sens.index]

    im1 = ax1.imshow(sens_pivot.values, cmap='YlOrRd', aspect='auto')
    ax1.set_xticks(range(len(sens_pivot.columns)))
    ax1.set_yticks(range(len(sens_pivot.index)))
    ax1.set_xticklabels(sens_pivot.columns, fontsize=10)
    ax1.set_yticklabels(sens_pivot.index, fontsize=10)
    ax1.set_title('Parameter Sensitivity by Observable', fontweight='bold')

    # Add values
    for i in range(len(sens_pivot.index)):
        for j in range(len(sens_pivot.columns)):
            val = sens_pivot.values[i, j]
            color = 'white' if val > 0.5 else 'black'
            ax1.text(j, i, f'{val:.2f}', ha='center', va='center',
                    color=color, fontsize=9)

    plt.colorbar(im1, ax=ax1, label='Normalized Sensitivity')

    # 2. Correlation matrix
    ax2 = fig.add_subplot(gs[0, 1])

    im2 = ax2.imshow(correlation_df.values, cmap='RdBu_r',
                     vmin=-1, vmax=1, aspect='auto')
    ax2.set_xticks(range(len(param_names)))
    ax2.set_yticks(range(len(param_names)))
    ax2.set_xticklabels(param_names, rotation=45, ha='right', fontsize=10)
    ax2.set_yticklabels(param_names, fontsize=10)
    ax2.set_title('Parameter Correlation Matrix (log scale)', fontweight='bold')

    # Add values
    for i in range(len(param_names)):
        for j in range(len(param_names)):
            val = correlation_df.values[i, j]
            color = 'white' if abs(val) > 0.5 else 'black'
            ax2.text(j, i, f'{val:.2f}', ha='center', va='center',
                    color=color, fontsize=8)

    plt.colorbar(im2, ax=ax2, label='Correlation')

    # 3. Sensitivity ranking
    ax3 = fig.add_subplot(gs[1, 0])

    sens_values = [recommendation['sensitivity'][p] for p in param_names]
    colors = ['green' if p in recommendation['to_estimate'] else 'lightgray'
              for p in param_names]

    sorted_indices = np.argsort(sens_values)[::-1]
    sorted_params = [param_names[i] for i in sorted_indices]
    sorted_sens = [sens_values[i] for i in sorted_indices]
    sorted_colors = [colors[i] for i in sorted_indices]

    ax3.barh(range(len(sorted_params)), sorted_sens, color=sorted_colors, alpha=0.7)
    ax3.set_yticks(range(len(sorted_params)))
    ax3.set_yticklabels(sorted_params, fontsize=10)
    ax3.set_xlabel('Average Sensitivity', fontweight='bold')
    ax3.set_title('Parameter Sensitivity Ranking', fontweight='bold')
    ax3.axvline(0.1, color='red', linestyle='--', linewidth=1, alpha=0.5,
               label='Low sensitivity threshold')
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='x')

    # 4. Recommendation summary
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')

    # Create text summary
    summary_text = "RECOMMENDATION SUMMARY\n" + "="*40 + "\n\n"
    summary_text += f"Parameters to ESTIMATE ({len(recommendation['to_estimate'])}):\n"
    for param in recommendation['to_estimate']:
        sens = recommendation['sensitivity'][param]
        summary_text += f"  ✓ {param:12s} (sens: {sens:.4f})\n"

    summary_text += f"\nParameters to FIX ({len(recommendation['to_fix'])}):\n"
    for param, reason in recommendation['reasoning']:
        summary_text += f"  ✗ {param:12s}\n"
        summary_text += f"    → {reason}\n"

    if recommendation['correlated_pairs']:
        summary_text += f"\nHighly Correlated Pairs:\n"
        for p1, p2, corr in recommendation['correlated_pairs']:
            summary_text += f"  {p1} ↔ {p2} (r={corr:.3f})\n"

    ax4.text(0.1, 0.95, summary_text, transform=ax4.transAxes,
            fontsize=9, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\nSaved figure: {save_path}")

    plt.show()
