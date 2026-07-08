"""
Debug why fixed parameters are causing issues with NLLH calculation.
"""

import numpy as np
import pandas as pd
from typing import Dict


def diagnose_fixed_parameter_issue(
    estimator,
    parameters_to_estimate: list,
    true_params: dict,
    verbose: bool = True
):
    """
    Diagnose issues with fixed parameter handling.

    Shows:
    1. Which parameters are fixed vs estimated
    2. What values the estimator is actually using for fixed params
    3. Whether those match the true values
    4. How to construct proper true_param_dict for testing
    """

    if verbose:
        print("="*80)
        print("FIXED PARAMETER DIAGNOSTIC")
        print("="*80)

    # 1. Identify fixed vs estimated parameters
    fixed_params = []
    estimated_params = []

    for param in parameters_to_estimate:
        param_id = param['parameterId']
        if param_id.startswith('sigma_'):
            continue  # Skip noise params for now

        if param['estimate'] == 0:
            fixed_params.append({
                'id': param_id,
                'nominal': param['nominalValue'],
                'true': true_params.get(param_id, None)
            })
        else:
            estimated_params.append({
                'id': param_id,
                'nominal': param['nominalValue'],
                'true': true_params.get(param_id, None)
            })

    if verbose:
        print(f"\nParameter counts:")
        print(f"  Fixed:     {len(fixed_params)}")
        print(f"  Estimated: {len(estimated_params)}")
        print(f"  Total:     {len(fixed_params) + len(estimated_params)}")

    # 2. Check what estimator knows about
    if verbose:
        print(f"\nEstimator's estimable_params:")
        print(f"  {estimator.estimable_params}")

    # 3. Check estimator's param_info for fixed params
    if verbose:
        print(f"\nFixed parameters in detail:")
        print(f"  {'Parameter':<15} {'In Config':<15} {'In Estimator':<15} {'True Value':<15} {'Match?':<10}")
        print(f"  {'-'*75}")

    issues = []
    for fp in fixed_params:
        param_id = fp['id']
        config_nominal = fp['nominal']
        estimator_nominal = estimator.param_info.get(param_id, {}).get('nominal', None)
        true_val = fp['true']

        if estimator_nominal is None:
            match = "N/A"
            issues.append(f"{param_id}: Not in estimator!")
        elif true_val is None:
            match = "No true"
            issues.append(f"{param_id}: No true value to compare")
        else:
            rel_diff = abs(estimator_nominal - true_val) / abs(true_val)
            match = "✓" if rel_diff < 1e-6 else f"✗ ({rel_diff*100:.1f}%)"
            if rel_diff >= 1e-6:
                issues.append(f"{param_id}: Estimator nominal ({estimator_nominal:.6e}) != true ({true_val:.6e})")

        if verbose:
            print(f"  {param_id:<15} {config_nominal:<15.6e} {str(estimator_nominal):<15} {str(true_val):<15} {match:<10}")

    # 4. Show how to properly construct true_param_dict
    if verbose:
        print(f"\nHow to construct true_param_dict for testing:")
        print(f"  Option 1: Include ALL parameters (both fixed and estimated)")
        print(f"  Option 2: Include only ESTIMATED parameters")
        print(f"  → Estimator will add fixed params at nominal values automatically")

    # Construct proper true_param_dict (only estimated params)
    true_param_dict_estimated_only = {}
    for ep in estimated_params:
        if ep['true'] is not None:
            true_param_dict_estimated_only[ep['id']] = ep['true']

    # Add noise parameters
    for param in parameters_to_estimate:
        if param['parameterId'].startswith('sigma_') and param['estimate'] == 1:
            param_id = param['parameterId']
            if param_id in true_params:
                true_param_dict_estimated_only[param_id] = true_params[param_id]

    if verbose:
        print(f"\nCorrect true_param_dict (estimated params only):")
        for k, v in true_param_dict_estimated_only.items():
            print(f"  {k}: {v:.6e}")

    # 5. Test both approaches
    print(f"\n{'='*80}")
    print("TESTING NLLH WITH DIFFERENT CONFIGURATIONS")
    print(f"{'='*80}")

    # Test 1: Only estimated params in dict
    print(f"\nTest 1: true_param_dict with only ESTIMATED parameters")
    try:
        vec1 = estimator.params_to_vector(true_param_dict_estimated_only)
        print(f"  Vector shape: {vec1.shape}")
        print(f"  Vector values: {vec1}")

        # Reconstruct to see what's actually used
        reconstructed1 = estimator.vector_to_params(vec1)
        print(f"  Reconstructed params:")
        for k, v in reconstructed1.items():
            if not k.startswith('sigma_'):
                is_fixed = k in [fp['id'] for fp in fixed_params]
                marker = "[FIXED]" if is_fixed else "[ESTIMATED]"
                true_val = true_params.get(k, None)
                if true_val:
                    diff = abs(v - true_val) / abs(true_val) * 100
                    print(f"    {k:<15} = {v:<15.6e} {marker:<12} (diff from true: {diff:5.1f}%)")
                else:
                    print(f"    {k:<15} = {v:<15.6e} {marker}")

        nllh1 = estimator.compute_nllh(vec1)
        print(f"  NLLH: {nllh1:.4f}")
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        nllh1 = None

    # Test 2: All params in dict (how user was doing it)
    print(f"\nTest 2: true_param_dict with ALL parameters (fixed + estimated)")
    true_param_dict_all = true_param_dict_estimated_only.copy()
    for fp in fixed_params:
        if fp['true'] is not None:
            true_param_dict_all[fp['id']] = fp['true']

    try:
        vec2 = estimator.params_to_vector(true_param_dict_all)
        print(f"  Vector shape: {vec2.shape}")
        print(f"  Vector values: {vec2}")

        # Should be same as vec1 (fixed params ignored in vector)
        if nllh1 is not None and np.allclose(vec1, vec2):
            print(f"  ✓ Vector same as Test 1 (fixed params ignored as expected)")

        nllh2 = estimator.compute_nllh(vec2)
        print(f"  NLLH: {nllh2:.4f}")

        if nllh1 is not None and abs(nllh2 - nllh1) < 0.01:
            print(f"  ✓ NLLH same as Test 1")
        else:
            print(f"  ✗ NLLH different from Test 1!")
    except Exception as e:
        print(f"  ✗ Failed: {e}")

    # 6. Summary and recommendations
    print(f"\n{'='*80}")
    print("DIAGNOSIS SUMMARY")
    print(f"{'='*80}")

    if len(issues) > 0:
        print(f"\n⚠ {len(issues)} issue(s) found:")
        for issue in issues:
            print(f"  • {issue}")

        print(f"\nRECOMMENDATION:")
        print(f"  The estimator has cached the OLD nominal values from when it was initialized.")
        print(f"  To fix this, you need to:")
        print(f"  1. Go back to Section 4.1")
        print(f"  2. Set nominalValue = true value for all FIXED parameters")
        print(f"  3. Re-run Section 4.3 (Create PETAb DataFrames)")
        print(f"  4. Re-run Section 5.1 (Initialize Parameter Estimator)")
        print(f"  5. Then test NLLH in Section 5.2")

        print(f"\nOR use this quick fix:")
        print(f"  estimator.param_info['{issues[0].split(':')[0]}']['nominal'] = {fixed_params[0]['true']}")
        print(f"  (Repeat for each fixed parameter)")

    else:
        print(f"\n✓ No issues found!")
        print(f"  All fixed parameters have correct nominal values in the estimator.")

    return {
        'fixed_params': fixed_params,
        'estimated_params': estimated_params,
        'issues': issues,
        'true_param_dict_correct': true_param_dict_estimated_only,
    }


def quick_fix_estimator_nominals(estimator, parameters_to_estimate, true_params):
    """
    Directly update the estimator's param_info with correct nominal values.

    This is a quick fix that doesn't require re-initializing the estimator.
    """
    print("Applying quick fix to estimator.param_info...")

    fixed = 0
    for param in parameters_to_estimate:
        param_id = param['parameterId']

        if param['estimate'] == 0 and param_id in true_params:
            # This is a fixed parameter
            true_val = true_params[param_id]
            old_val = estimator.param_info.get(param_id, {}).get('nominal', None)

            if old_val is not None and abs(old_val - true_val) / abs(true_val) > 1e-6:
                estimator.param_info[param_id]['nominal'] = true_val
                print(f"  ✓ {param_id}: {old_val:.6e} → {true_val:.6e}")
                fixed += 1

    if fixed == 0:
        print("  No fixes needed (all nominal values already correct)")
    else:
        print(f"\n✓ Fixed {fixed} parameter(s) in estimator.param_info")
        print("  You can now test NLLH with correct values (no need to re-run initialization)")

    return fixed
