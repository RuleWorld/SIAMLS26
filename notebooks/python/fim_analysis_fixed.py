"""
Fixed FIM Analysis for Parameter Estimation

This module provides corrected FIM computation that:
1. Excludes noise parameters (sigma_*)
2. Properly handles parameter transformations (log10 space)
3. Uses adaptive step sizes for numerical stability
4. Transforms uncertainties correctly to natural space
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
from .parameter_estimator import ParameterEstimator


def compute_fim_robust(estimator: ParameterEstimator,
                       param_vector: np.ndarray,
                       exclude_noise: bool = True,
                       delta_rel: float = 1e-3,
                       delta_abs: float = 1e-6) -> Tuple[np.ndarray, List[str], List[int]]:
    """
    Compute Fisher Information Matrix with robust numerics.

    Key improvements:
    - Excludes noise parameters by default (they're not mechanistic)
    - Uses adaptive step sizes
    - FIM is in estimation space (matches optimization)

    Args:
        estimator: ParameterEstimator instance
        param_vector: Parameter values in ESTIMATION SPACE (log10 for most params)
        exclude_noise: If True, exclude sigma_* parameters from FIM
        delta_rel: Relative step size (fraction of parameter value)
        delta_abs: Absolute step size (for near-zero values)

    Returns:
        Tuple of:
        - FIM matrix (n_model_params × n_model_params)
        - List of parameter names included
        - List of indices in original param_vector
    """
    # Determine which parameters to include
    if exclude_noise:
        # Exclude noise parameters (sigma_*)
        param_indices = [i for i, name in enumerate(estimator.estimable_params)
                        if not name.startswith('sigma_')]
        param_names = [name for name in estimator.estimable_params
                      if not name.startswith('sigma_')]
    else:
        param_indices = list(range(len(estimator.estimable_params)))
        param_names = estimator.estimable_params.copy()

    n_params = len(param_indices)
    fim = np.zeros((n_params, n_params))

    print(f"Computing FIM for {n_params} parameters (excluding noise: {exclude_noise})")
    print(f"Parameters: {param_names}")
    print(f"Total function evaluations: ~{1 + 2*n_params + n_params*(n_params+1)//2}")

    # Helper function to evaluate objective with subset of parameters perturbed
    def eval_nllh(perturbations: Dict[int, float]) -> float:
        """Evaluate NLLH with specified parameter perturbations."""
        perturbed_vector = param_vector.copy()
        for idx, delta_val in perturbations.items():
            perturbed_vector[idx] += delta_val
        return estimator.compute_nllh(perturbed_vector)

    # Compute center value
    f_center = estimator.compute_nllh(param_vector)
    print(f"Center NLLH: {f_center:.4f}")

    # Compute Hessian using central finite differences
    print("\nComputing Hessian...")
    for idx_i, i in enumerate(param_indices):
        # Adaptive step size
        # In log space, a fixed relative step works better than scaling by value
        step_i = max(delta_rel, delta_abs)

        # Diagonal element (second derivative)
        f_plus = eval_nllh({i: step_i})
        f_minus = eval_nllh({i: -step_i})

        fim[idx_i, idx_i] = (f_plus - 2*f_center + f_minus) / (step_i**2)

        # Off-diagonal elements (mixed derivatives)
        for idx_j, j in enumerate(param_indices[idx_i+1:], start=idx_i+1):
            step_j = max(delta_rel, delta_abs)

            # Four-point stencil for mixed derivative
            f_pp = eval_nllh({i: step_i, j: step_j})
            f_pm = eval_nllh({i: step_i, j: -step_j})
            f_mp = eval_nllh({i: -step_i, j: step_j})
            f_mm = eval_nllh({i: -step_i, j: -step_j})

            mixed_deriv = (f_pp - f_pm - f_mp + f_mm) / (4 * step_i * step_j)
            fim[idx_i, idx_j] = mixed_deriv
            fim[idx_j, idx_i] = mixed_deriv  # Symmetric

        if (idx_i + 1) % max(1, n_params // 5) == 0:
            print(f"  Completed {idx_i + 1}/{n_params} rows")

    print("FIM computation complete!")

    # Check for numerical issues
    eigenvals = np.linalg.eigvalsh(fim)
    n_negative = np.sum(eigenvals < -1e-10)
    n_near_zero = np.sum(np.abs(eigenvals) < 1e-10)

    if n_negative > 0:
        print(f"\nWarning: FIM has {n_negative} negative eigenvalues (should be positive semidefinite)")
        print("This suggests numerical issues in Hessian computation.")
    if n_near_zero > 0:
        print(f"Note: FIM has {n_near_zero} near-zero eigenvalues")
        print("This indicates non-identifiability or strong parameter correlation.")

    return fim, param_names, param_indices


def transform_uncertainties_to_natural_space(
    param_vector: np.ndarray,
    param_names: List[str],
    param_cov: np.ndarray,
    estimator: ParameterEstimator
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Transform parameter uncertainties from estimation space to natural space.

    For log10-scaled parameters:
    - Estimation space: θ_est = log10(θ_nat)
    - Natural space: θ_nat = 10^θ_est
    - By delta method: σ_nat = θ_nat * ln(10) * σ_est

    Args:
        param_vector: Parameter values in estimation space
        param_names: Parameter names (corresponding to param_vector subset)
        param_cov: Covariance matrix in estimation space
        estimator: ParameterEstimator instance

    Returns:
        Tuple of (natural_values, natural_std_errors)
    """
    n_params = len(param_names)
    natural_values = np.zeros(n_params)
    natural_std = np.zeros(n_params)

    # Get original indices in full parameter vector
    param_indices = [estimator.estimable_params.index(name) for name in param_names]

    for i, (name, idx) in enumerate(zip(param_names, param_indices)):
        param_info = estimator.param_info[name]
        scale = param_info['scale']

        # Get value in estimation space
        est_value = param_vector[idx]

        if scale == 'log10':
            # Transform to natural space
            nat_value = 10 ** est_value

            # Transform std error using delta method
            # d(10^x)/dx = 10^x * ln(10)
            # So: σ_nat = |d(10^x)/dx| * σ_est = 10^x * ln(10) * σ_est
            est_std = np.sqrt(param_cov[i, i])
            nat_std_err = nat_value * np.log(10) * est_std

            natural_values[i] = nat_value
            natural_std[i] = nat_std_err
        else:
            # Linear scale - no transformation needed
            natural_values[i] = est_value
            natural_std[i] = np.sqrt(param_cov[i, i])

    return natural_values, natural_std


def analyze_fim(estimator: ParameterEstimator,
                best_result: dict,
                true_param_dict: dict,
                results_dir: str,
                timestamp: str):
    """
    Complete FIM analysis with proper handling of parameter spaces.

    Args:
        estimator: ParameterEstimator instance
        best_result: Dictionary with 'param_vector' and 'param_dict'
        true_param_dict: Dictionary of true parameter values
        results_dir: Directory to save results
        timestamp: Timestamp string for file naming
    """
    print("="*80)
    print("FISHER INFORMATION MATRIX ANALYSIS")
    print("="*80)

    # 1. Compute FIM (excluding noise parameters)
    fim, param_names, param_indices = compute_fim_robust(
        estimator,
        best_result['param_vector'],
        exclude_noise=True,
        delta_rel=1e-3
    )

    print(f"\nFIM shape: {fim.shape}")
    print(f"FIM symmetry: max|FIM - FIM.T| = {np.max(np.abs(fim - fim.T)):.2e}")

    # 2. Eigenvalue analysis
    print("\n" + "="*80)
    print("EIGENVALUE ANALYSIS")
    print("="*80)

    eigenvalues, eigenvectors = np.linalg.eigh(fim)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Condition number
    if eigenvalues[-1] > 1e-10:
        condition_number = eigenvalues[0] / eigenvalues[-1]
    else:
        condition_number = np.inf

    print(f"\n{'Index':<8} {'Eigenvalue':<15} {'Relative':<15}")
    print("-"*45)
    for i, ev in enumerate(eigenvalues):
        rel = ev / eigenvalues[0] if eigenvalues[0] > 0 else 0
        print(f"{i+1:<8} {ev:<15.4e} {rel:<15.4e}")

    print(f"\nCondition number: {condition_number:.2e}")

    # 3. Compute parameter uncertainties
    print("\n" + "="*80)
    print("PARAMETER UNCERTAINTIES")
    print("="*80)

    # Try to invert FIM
    try:
        # For near-singular matrices, use regularization
        if condition_number > 1e6:
            print("Warning: FIM is poorly conditioned. Using regularization.")
            # Add small regularization to diagonal
            reg_fim = fim + 1e-6 * np.mean(np.diag(fim)) * np.eye(len(fim))
            param_cov = np.linalg.inv(reg_fim)
        else:
            param_cov = np.linalg.inv(fim)
    except np.linalg.LinAlgError:
        print("Error: FIM is singular. Using pseudo-inverse.")
        param_cov = np.linalg.pinv(fim)

    # Transform to natural space
    natural_values, natural_std = transform_uncertainties_to_natural_space(
        best_result['param_vector'],
        param_names,
        param_cov,
        estimator
    )

    # Print results
    print(f"\n{'Parameter':<15} {'Estimate':<15} {'Std Error':<15} {'Rel Error %':<15} {'True Value':<15} {'Fit Error %':<15}")
    print("-"*95)

    results_list = []
    for i, name in enumerate(param_names):
        est_val = natural_values[i]
        std_err = natural_std[i]
        rel_err = (std_err / est_val * 100) if est_val != 0 else np.inf
        true_val = true_param_dict.get(name, np.nan)

        if not np.isnan(true_val):
            fit_err = abs(est_val - true_val) / true_val * 100
        else:
            fit_err = np.nan

        print(f"{name:<15} {est_val:<15.6e} {std_err:<15.6e} {rel_err:<15.2f} {true_val:<15.6e} {fit_err:<15.2f}")

        results_list.append({
            'parameterId': name,
            'estimate': est_val,
            'std_error': std_err,
            'rel_error_pct': rel_err,
            'true_value': true_val,
            'fit_error_pct': fit_err
        })

    # 4. Create visualizations
    print("\n" + "="*80)
    print("CREATING VISUALIZATIONS")
    print("="*80)

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # 1. FIM heatmap
    ax1 = fig.add_subplot(gs[0, :])
    # Use symmetric log scale for better visualization
    fim_plot = np.sign(fim) * np.log10(np.abs(fim) + 1e-10)
    im = ax1.imshow(fim_plot, cmap='RdBu_r', aspect='auto')
    ax1.set_xticks(range(len(param_names)))
    ax1.set_yticks(range(len(param_names)))
    ax1.set_xticklabels(param_names, rotation=45, ha='right', fontsize=10)
    ax1.set_yticklabels(param_names, fontsize=10)
    ax1.set_title('Fisher Information Matrix (sign×log₁₀|FIM|)', fontweight='bold', fontsize=12)
    plt.colorbar(im, ax=ax1, label='sign×log₁₀|FIM|')

    # 2. Eigenvalue spectrum
    ax2 = fig.add_subplot(gs[1, 0])
    # Use absolute value on log scale, mark negative eigenvalues
    colors = ['red' if ev < 0 else 'blue' for ev in eigenvalues]
    ax2.semilogy(range(1, len(eigenvalues)+1), np.abs(eigenvalues), 'o-',
                linewidth=2, markersize=8, color='blue')
    for i, (ev, c) in enumerate(zip(eigenvalues, colors)):
        if ev < 0:
            ax2.plot(i+1, abs(ev), 'ro', markersize=8, label='Negative' if i == 0 else '')
    ax2.set_xlabel('Eigenvalue Index', fontweight='bold')
    ax2.set_ylabel('|Eigenvalue| (log scale)', fontweight='bold')
    ax2.set_title('FIM Eigenvalue Spectrum', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    if any(ev < 0 for ev in eigenvalues):
        ax2.legend()

    # 3. Parameter uncertainties
    ax3 = fig.add_subplot(gs[1, 1])
    rel_errors = [(natural_std[i] / natural_values[i] * 100)
                  for i in range(len(param_names))]
    colors = ['green' if re < 30 else 'orange' if re < 100 else 'red'
              for re in rel_errors]
    ax3.barh(range(len(param_names)), rel_errors, color=colors, alpha=0.7)
    ax3.set_yticks(range(len(param_names)))
    ax3.set_yticklabels(param_names, fontsize=10)
    ax3.set_xlabel('Relative Uncertainty (%)', fontweight='bold')
    ax3.set_title('Parameter Identifiability (from FIM)', fontweight='bold')
    ax3.axvline(30, color='orange', linestyle='--', linewidth=1, alpha=0.5)
    ax3.axvline(100, color='red', linestyle='--', linewidth=1, alpha=0.5)
    ax3.grid(True, alpha=0.3, axis='x')

    # 4. Correlation matrix
    ax4 = fig.add_subplot(gs[2, :])
    param_std_est = np.sqrt(np.diag(param_cov))
    param_corr = param_cov / np.outer(param_std_est, param_std_est)

    im = ax4.imshow(param_corr, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax4.set_xticks(range(len(param_names)))
    ax4.set_yticks(range(len(param_names)))
    ax4.set_xticklabels(param_names, rotation=45, ha='right', fontsize=10)
    ax4.set_yticklabels(param_names, fontsize=10)
    ax4.set_title('Parameter Correlation Matrix', fontweight='bold', fontsize=12)

    # Add correlation values
    for i in range(len(param_names)):
        for j in range(len(param_names)):
            val = param_corr[i, j]
            color = 'white' if abs(val) > 0.5 else 'black'
            ax4.text(j, i, f'{val:.2f}', ha="center", va="center",
                    color=color, fontsize=8)

    plt.colorbar(im, ax=ax4, label='Correlation')

    plt.savefig(f'{results_dir}/fim_analysis_fixed_{timestamp}.png',
                dpi=300, bbox_inches='tight')
    print(f"Saved: {results_dir}/fim_analysis_fixed_{timestamp}.png")
    plt.show()

    # 5. Save results
    results_df = pd.DataFrame(results_list)
    results_file = f'{results_dir}/fim_identifiability_{timestamp}.csv'
    results_df.to_csv(results_file, index=False)
    print(f"Saved: {results_file}")

    return {
        'fim': fim,
        'param_names': param_names,
        'eigenvalues': eigenvalues,
        'eigenvectors': eigenvectors,
        'param_cov': param_cov,
        'natural_values': natural_values,
        'natural_std': natural_std,
        'results_df': results_df
    }
