"""
Parameter estimation for BNGL models using PETAb-formatted data.

This module provides the ParameterEstimator class which computes goodness-of-fit
(negative log-likelihood) for BNGL models given experimental data in PETAb format.

REFACTORED VERSION (2025):
- No automatic normalization/scaling - uses parametric transformations from observables table
- Single simulator mode with condition parameters set via reference_values
- Parses and evaluates observable formulas (e.g., "scale * observable + baseline")
- Supports SEM-based noise from observableParameters
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from typing import Optional, Dict, List, Tuple, Union
from .bngl_simulator import BNGLSimulator


class ParameterEstimator:
    """
    Parameter estimation for BNGL models using PETAb-formatted data.

    This class computes goodness-of-fit (negative log-likelihood) for a given
    set of parameter values by:
    1. Running simulations for each experimental condition
    2. Applying observable transformations from PETAb observables table
    3. Computing likelihood using noise model (SEM or parametric sigma)

    Uses a single simulator for all conditions, with condition-specific parameters
    set dynamically via reference_values at simulation time.
    """

    def __init__(self,
                 bngl_file: str,
                 conditions_df: pd.DataFrame,
                 measurements_df: pd.DataFrame,
                 observables_df: pd.DataFrame,
                 parameters_df: pd.DataFrame,
                 estimable_params: Optional[List[str]] = None):
        """
        Initialize parameter estimator.

        Args:
            bngl_file: Path to BNGL model file
            conditions_df: PETAb conditions table
            measurements_df: PETAb measurements table
            observables_df: PETAb observables table
            parameters_df: PETAb parameters table
            estimable_params: List of parameter names to estimate (subset)
                             If None, uses all parameters with estimate=1
        """
        self.bngl_file = bngl_file
        self.conditions_df = conditions_df
        self.measurements_df = measurements_df
        self.observables_df = observables_df
        self.parameters_df = parameters_df

        # Determine which parameters to estimate
        if estimable_params is None:
            # Use all parameters marked for estimation
            self.estimable_params = parameters_df[
                parameters_df['estimate'] == 1
            ]['parameterId'].tolist()
        else:
            self.estimable_params = estimable_params

        # Create parameter info dict for bounds and scales
        self.param_info = {}
        for _, row in parameters_df.iterrows():
            self.param_info[row['parameterId']] = {
                'scale': row['parameterScale'],
                'lower': row['lowerBound'],
                'upper': row['upperBound'],
                'nominal': row['nominalValue']
            }

        # Extract unique timepoints for simulation
        self.timepoints = np.sort(measurements_df['time'].unique())

        # Initialize simulator (single simulator mode)
        self.simulator = None  # Created on first use
        self._last_params = {}  # Cache of last parameter values

        print(f"Initialized estimator with {len(self.estimable_params)} estimable parameters")
        print(f"Timepoints: {self.timepoints}")

    def _get_simulator(self, param_dict: Dict[str, float]) -> BNGLSimulator:
        """
        Create or update simulator with given parameter values.

        Args:
            param_dict: Dictionary of parameter values (may include non-model params)

        Returns:
            BNGLSimulator instance
        """
        # Filter out non-model parameters (scaling, baseline, and noise parameters)
        # These are PETAb-specific and not used in the BNGL ODE simulation
        # baseline_ parameters are only used in observable formulas, not in reactions
        # They will be available in param_dict for observable evaluation
        model_params = {
            k: v for k, v in param_dict.items()
            if not k.startswith('scale_')
            and not k.startswith('baseline_')
            and not k.startswith('sigma_')
        }

        # Single simulator mode
        if self.simulator is None:
            # Create simulator on first use
            self.simulator = BNGLSimulator(self.bngl_file, param_values=model_params)
            self._last_params = model_params.copy()
        else:
            # Update parameters that have changed
            for param_name, value in model_params.items():
                if param_name not in self._last_params or \
                   self._last_params[param_name] != value:
                    self.simulator.set_parameter(param_name, value)

            # Update cache
            self._last_params = model_params.copy()

        return self.simulator

    def _eval_observable_formula(self, formula: str,
                                 sim_result,
                                 param_dict: Dict[str, float],
                                 time_idx: int) -> float:
        """
        Evaluate an observable formula from PETAb observables table.

        Formula can reference:
        - Model observables (from simulation result)
        - Parameters (from param_dict)
        - Standard math operations

        Examples:
            "pSTAT1" -> sim_result['pSTAT1'][time_idx]
            "scale_pSTAT1 * pSTAT1 + baseline_pSTAT1" -> parametric transformation

        Args:
            formula: Observable formula string from PETAb
            sim_result: Simulation result (RoadRunner NamedArray or dict)
            param_dict: Parameter dictionary
            time_idx: Time index in simulation result

        Returns:
            Evaluated observable value
        """
        # Create namespace for evaluation
        namespace = {}

        # Add all parameters
        namespace.update(param_dict)

        # Add model observables at this timepoint
        # RoadRunner NamedArray has colnames attribute
        if hasattr(sim_result, 'colnames'):
            # RoadRunner NamedArray
            for col_name in sim_result.colnames:
                if col_name.lower() != 'time':
                    namespace[col_name] = sim_result[col_name][time_idx]
        else:
            # Regular dictionary
            for obs_name, obs_values in sim_result.items():
                if obs_name != 'time' and isinstance(obs_values, np.ndarray):
                    namespace[obs_name] = obs_values[time_idx]

        # Add numpy functions
        namespace['log'] = np.log
        namespace['log10'] = np.log10
        namespace['exp'] = np.exp
        namespace['sqrt'] = np.sqrt

        # Evaluate formula
        try:
            result = eval(formula, {"__builtins__": {}}, namespace)
            return float(result)
        except Exception as e:
            raise ValueError(f"Failed to evaluate observable formula '{formula}': {e}")

    def _parse_observable_parameters(self, obs_params: str) -> Dict[str, float]:
        """
        Parse observableParameters field from measurements table.

        Can be:
        - Empty/NaN: No additional parameters
        - "value1;value2;..." : Semicolon-separated values (e.g., "SD;N" for SEM calculation)

        Args:
            obs_params: observableParameters string

        Returns:
            Dictionary with parsed values (empty if none)
        """
        if pd.isna(obs_params) or obs_params == '':
            return {}

        # Split by semicolon
        parts = str(obs_params).split(';')

        # Try to parse as SD;N format (for SEM)
        if len(parts) == 2:
            try:
                sd = float(parts[0])
                n = int(parts[1])
                return {'SD': sd, 'N': n, 'SEM': sd / np.sqrt(n)}
            except ValueError:
                pass

        # Otherwise return raw parts as strings
        return {f'param{i}': p for i, p in enumerate(parts)}

    def vector_to_params(self, param_vector: np.ndarray) -> Dict[str, float]:
        """
        Convert parameter vector to dictionary.

        Handles log10 scaling as specified in PETAb.

        Args:
            param_vector: Array of parameter values (in estimation space)

        Returns:
            Dictionary mapping parameter names to values (in natural space)
        """
        param_dict = {}
        for i, param_name in enumerate(self.estimable_params):
            value = param_vector[i]
            scale = self.param_info[param_name]['scale']

            # Transform from estimation space to natural space
            if scale == 'log10':
                param_dict[param_name] = 10 ** value
            else:
                param_dict[param_name] = value

        # Add non-estimable parameters at nominal values
        for param_name, info in self.param_info.items():
            if param_name not in param_dict:
                param_dict[param_name] = info['nominal']

        return param_dict

    def params_to_vector(self, param_dict: Dict[str, float]) -> np.ndarray:
        """
        Convert parameter dictionary to vector.

        Handles log10 scaling as specified in PETAb.

        Args:
            param_dict: Dictionary of parameter values (in natural space)

        Returns:
            Array of parameter values (in estimation space)
        """
        param_vector = np.zeros(len(self.estimable_params))
        for i, param_name in enumerate(self.estimable_params):
            value = param_dict[param_name]
            scale = self.param_info[param_name]['scale']

            # Transform from natural space to estimation space
            if scale == 'log10':
                param_vector[i] = np.log10(value)
            else:
                param_vector[i] = value

        return param_vector

    def get_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get parameter bounds in estimation space.

        Returns:
            Tuple of (lower_bounds, upper_bounds) arrays
        """
        lower = np.zeros(len(self.estimable_params))
        upper = np.zeros(len(self.estimable_params))

        for i, param_name in enumerate(self.estimable_params):
            info = self.param_info[param_name]

            if info['scale'] == 'log10':
                lower[i] = np.log10(info['lower'])
                upper[i] = np.log10(info['upper'])
            else:
                lower[i] = info['lower']
                upper[i] = info['upper']

        return lower, upper

    def compute_nllh(self, param_vector: np.ndarray, verbose: bool = False) -> float:
        """
        Compute negative log-likelihood for given parameters.

        This is the objective function to minimize.

        Args:
            param_vector: Array of parameter values (in estimation space)
            verbose: Print debug information

        Returns:
            Negative log-likelihood value
        """
        # Convert to parameter dict
        param_dict = self.vector_to_params(param_vector)

        # Compute log-likelihood contributions
        nllh = 0.0
        n_datapoints = 0

        # Group measurements by condition
        for _, cond_row in self.conditions_df.iterrows():
            cond_id = cond_row['conditionId']

            # Get condition-specific parameters
            cond_params = {}
            for col in cond_row.index:
                if col != 'conditionId' and col != 'conditionName' and pd.notna(cond_row[col]):
                    cond_params[col] = cond_row[col]

            # Get measurements for this condition
            cond_measurements = self.measurements_df[
                self.measurements_df['simulationConditionId'] == cond_id
            ]

            if len(cond_measurements) == 0:
                continue

            # Get unique measurement times for this condition
            meas_times = np.sort(cond_measurements['time'].unique())

            # Calculate simulation timepoints
            t_end = float(meas_times[-1])  # Convert to native Python float
            # Use enough steps to have good resolution (at least 2 points per minute)
            n_steps = int(max(int(t_end * 2), 100))

            # Get simulator
            simulator = self._get_simulator(param_dict)

            # Run simulation with condition-specific parameters set via reference_values
            try:
                result = simulator.simulate(
                    t_end=t_end,
                    n_steps=n_steps,
                    reference_values=cond_params,
                    reset=True
                )
            except Exception as e:
                if verbose:
                    print(f"Simulation failed for condition {cond_id}: {e}")
                return np.inf

            # Compute likelihood for each measurement
            for _, meas_row in cond_measurements.iterrows():
                obs_id = meas_row['observableId']
                time = meas_row['time']
                measurement = meas_row['measurement']

                # Find the timepoint index
                time_idx = np.argmin(np.abs(result['time'] - time))

                # Get observable formula
                obs_row = self.observables_df[
                    self.observables_df['observableId'] == obs_id
                ]

                if len(obs_row) == 0:
                    if verbose:
                        print(f"Warning: Observable {obs_id} not found in observables table")
                    return np.inf

                obs_formula = obs_row.iloc[0]['observableFormula']
                noise_formula = obs_row.iloc[0]['noiseFormula']

                # Evaluate observable formula to get prediction
                try:
                    prediction = self._eval_observable_formula(
                        obs_formula, result, param_dict, time_idx
                    )
                except Exception as e:
                    if verbose:
                        print(f"Failed to evaluate observable {obs_id}: {e}")
                    return np.inf

                # Get noise standard deviation
                # First check if SEM provided in observableParameters
                obs_params_str = meas_row.get('observableParameters', '')
                obs_params = self._parse_observable_parameters(obs_params_str)

                if 'SEM' in obs_params:
                    # Use SEM from observableParameters (already calculated SD/sqrt(N))
                    noise_std = obs_params['SEM']
                else:
                    # Evaluate noise formula (e.g., "sigma_pSTAT1")
                    try:
                        noise_std = self._eval_observable_formula(
                            noise_formula, result, param_dict, time_idx
                        )
                    except Exception as e:
                        if verbose:
                            print(f"Failed to evaluate noise formula {noise_formula}: {e}")
                        return np.inf

                # Gaussian log-likelihood: -0.5 * ((y - y_pred) / sigma)^2 - log(sigma) - 0.5*log(2*pi)
                # NLL (negative log-likelihood) is the negative of this
                residual = measurement - prediction
                nllh += 0.5 * (residual / noise_std) ** 2
                nllh += np.log(noise_std)
                nllh += 0.5 * np.log(2 * np.pi)

                n_datapoints += 1

        if verbose:
            print(f"NLLH = {nllh:.4f} ({n_datapoints} datapoints)")

        return nllh

    def simulate_conditions(self, param_vector: np.ndarray,
                          n_steps: int = 200) -> Dict[str, dict]:
        """
        Simulate all experimental conditions with given parameters.

        Args:
            param_vector: Array of parameter values (in estimation space)
            n_steps: Number of time steps for simulation (higher = smoother curves)

        Returns:
            Dictionary mapping condition IDs to simulation results
        """
        param_dict = self.vector_to_params(param_vector)
        results = {}

        for _, cond_row in self.conditions_df.iterrows():
            cond_id = cond_row['conditionId']

            # Get condition parameters
            cond_params = {}
            for col in cond_row.index:
                if col != 'conditionId' and col != 'conditionName' and pd.notna(cond_row[col]):
                    cond_params[col] = cond_row[col]

            # Get simulator
            simulator = self._get_simulator(param_dict)

            # Simulate with condition-specific parameters set via reference_values
            result = simulator.simulate(
                t_end=float(self.timepoints[-1]),  # Convert to native Python float
                n_steps=n_steps,
                reference_values=cond_params,
                reset=True
            )

            results[cond_id] = result

        return results

    def plot_fit(self, param_vector: np.ndarray,
                figsize: Optional[Tuple[int, int]] = None,
                save_path: Optional[str] = None):
        """
        Plot model fit vs experimental data.

        Args:
            param_vector: Array of parameter values (in estimation space)
            figsize: Figure size (optional, auto-calculated if not provided)
            save_path: Path to save figure (optional)
        """
        # Get parameter dict for observable evaluation
        param_dict = self.vector_to_params(param_vector)

        # Simulate all conditions
        results = self.simulate_conditions(param_vector, n_steps=200)

        # Get unique observables
        unique_obs = self.measurements_df['observableId'].unique()

        # Determine grid layout
        n_conditions = len(self.conditions_df)
        n_obs = len(unique_obs)

        if figsize is None:
            figsize = (5 * n_conditions, 4 * n_obs)

        fig, axes = plt.subplots(n_obs, n_conditions, figsize=figsize, squeeze=False)

        for obs_idx, obs_id in enumerate(unique_obs):
            for cond_idx, (_, cond_row) in enumerate(self.conditions_df.iterrows()):
                cond_id = cond_row['conditionId']
                ax = axes[obs_idx, cond_idx]

                result = results[cond_id]

                # Get observable formula
                obs_formula = self.observables_df[
                    self.observables_df['observableId'] == obs_id
                ]['observableFormula'].iloc[0]

                # Evaluate observable at all timepoints
                obs_values = []
                for t_idx in range(len(result['time'])):
                    obs_val = self._eval_observable_formula(
                        obs_formula, result, param_dict, t_idx
                    )
                    obs_values.append(obs_val)

                # Plot simulation
                ax.plot(result['time'], obs_values, '-', linewidth=2, label='Model')

                # Plot measurements
                cond_data = self.measurements_df[
                    (self.measurements_df['simulationConditionId'] == cond_id) &
                    (self.measurements_df['observableId'] == obs_id)
                ]

                if len(cond_data) > 0:
                    # Check for SEM in observableParameters
                    has_sem = False
                    sems = []
                    for _, row in cond_data.iterrows():
                        obs_params = self._parse_observable_parameters(
                            row.get('observableParameters', '')
                        )
                        if 'SEM' in obs_params:
                            sems.append(obs_params['SEM'])
                            has_sem = True
                        else:
                            sems.append(None)

                    if has_sem:
                        ax.errorbar(cond_data['time'], cond_data['measurement'],
                                  yerr=sems, fmt='o', markersize=6,
                                  capsize=4, label='Data ± SEM')
                    else:
                        ax.plot(cond_data['time'], cond_data['measurement'],
                              'o', markersize=6, label='Data')

                # Labels
                ax.set_xlabel('Time (min)', fontweight='bold')
                ax.set_ylabel(obs_id, fontweight='bold')
                ax.set_title(f"{cond_id}", fontweight='bold')
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved plot to {save_path}")

        plt.show()
