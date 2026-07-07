"""
LAT Aggregation: Steady-State Sampling and Analysis
=====================================================
Samples steady-state snapshots from NFsim simulations of the LAT-Grb2-SOS1
aggregation model, then analyzes the resulting cluster size distributions.

Each simulation (equilibration + sampling) runs as an independent worker
process via Python's multiprocessing.Pool, so all NUM_SIMULATIONS runs
proceed in parallel up to the number of available CPU cores.

Usage:
    Place this script in the same directory as LAT.bngl and run:
        python lat_sample_and_analyze.py

Output:
    <RESULTS_DIR>/equil{i}/         -- equilibration and snapshot output per run
    frac_LAT.png                    -- scatter plot of LAT cluster size distribution
    fg_LAT.png                      -- histogram of largest cluster size per run
"""

import bionetgen
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

# =============================================================================
# Hyperparameters
# =============================================================================

RESULTS_DIR          = "results"   # Root directory for all simulation output
BNGL_FILE            = "LAT"       # Name of the model file (without .bngl)
NUM_SIMULATIONS      = 30          # Number of independent equilibration runs
NUM_SAMPLES          = 100         # Number of steady-state snapshots per run
TIME_BETWEEN_SAMPLES = 10          # Simulation time between successive snapshots
SIM_TIMEOUT          = 300         # Max seconds to allow a single NFsim call (0 = no limit)
MAX_RETRIES          = 20          # Max consecutive failures before aborting a run
TOTAL_LAT            = 2e3         # Total LAT molecules (for normalizing fractions)
NUM_WORKERS          = mp.cpu_count()  # Parallel workers; defaults to all available cores


# =============================================================================
# Worker function: one full simulation (equilibration + sampling) for run i
# =============================================================================
# This function is executed independently for each simulation index i by the
# multiprocessing pool. Because each run writes to its own subdirectory, there
# are no shared file-system conflicts between workers.

def run_simulation(i):
    """
    Equilibrate from LAT.bngl, then take NUM_SAMPLES steady-state snapshots
    for simulation index i. Returns a status string for logging.
    """

    # ------------------------------------------------------------------
    # Step 1: Equilibration
    # ------------------------------------------------------------------
    # Run the model to steady state from the initial conditions in LAT.bngl.
    # Output goes to RESULTS_DIR/equil{i}/ so each worker has its own space.
    print(f"[Run {i}] Equilibrating...")
    bionetgen.run(f"{BNGL_FILE}.bngl", out=f"{RESULTS_DIR}/equil{i}")

    # Load the equilibrated model as a mutable BioNetGen model object
    # so we can advance it incrementally during sampling
    model = bionetgen.bngmodel(f"{RESULTS_DIR}/equil{i}/{BNGL_FILE}.bngl")

    # ------------------------------------------------------------------
    # Step 2: Steady-State Sampling
    # ------------------------------------------------------------------
    # Advance the simulation by TIME_BETWEEN_SAMPLES time units, record the
    # species snapshot, and repeat NUM_SAMPLES times. This produces a
    # Markov-chain-style sample of the steady-state distribution with a
    # decorrelation time of TIME_BETWEEN_SAMPLES.

    for samp in range(NUM_SAMPLES):

        # Advance the simulation by one sampling interval
        model.add_action(action_type="simulate", action_args={
            "method":  "'nf'",
            "t_end":   str(TIME_BETWEEN_SAMPLES),
            "n_steps": str(TIME_BETWEEN_SAMPLES)
        })
        # Persist the current model state (rates, parameters, etc.) to disk
        model.add_action(action_type='writeModel', action_args={'overwrite': '1'})

        # Build a pruned species block containing only species with count > 0.
        # This is critical for NFsim: carrying zero-count species forward wastes
        # memory and can cause unexpected behavior in rare configurations.
        current_species = model.species
        new_block = bionetgen.modelapi.blocks.SpeciesBlock()
        for spec in current_species:
            if int(spec.count) > 0:
                new_block.add_species(spec.pattern, spec.count)
        model.add_species_block(new_block)

        # Write this snapshot as its own .bngl file for reproducibility
        snapshot_path = f"{RESULTS_DIR}/equil{i}/{BNGL_FILE}_samp{samp}.bngl"
        with open(snapshot_path, "w") as f:
            f.write(str(model))

        # --- Run the snapshot (with retry on failure) ---
        # NFsim can occasionally error on rare species graph configurations;
        # retrying is safe because each snapshot is a fresh independent call.
        run_dir = f"{RESULTS_DIR}/equil{i}/run{samp}"
        attempts = 0
        while True:
            try:
                bionetgen.run(snapshot_path, out=run_dir, timeout=SIM_TIMEOUT)
                break
            except Exception as e:
                attempts += 1
                print(f"[Run {i}] Sample {samp} run failed (attempt {attempts}): {e}")
                if attempts > MAX_RETRIES:
                    return f"[Run {i}] FAILED to run sample {samp} after {MAX_RETRIES} attempts — aborting."
                continue

        # --- Reload model from completed run (with retry on failure) ---
        # Model reloading can also fail if the output file is malformed;
        # retry independently of the run retry above.
        attempts = 0
        while True:
            try:
                model = bionetgen.bngmodel(f"{run_dir}/{BNGL_FILE}_samp{samp}.bngl")
                break
            except Exception as e:
                attempts += 1
                print(f"[Run {i}] Sample {samp} model load failed (attempt {attempts}): {e}")
                if attempts > MAX_RETRIES:
                    return f"[Run {i}] FAILED to load model at sample {samp} after {MAX_RETRIES} attempts — aborting."
                continue

        print(f"[Run {i}] Sample {samp + 1}/{NUM_SAMPLES} complete")

    return f"[Run {i}] Finished successfully."


# =============================================================================
# Analysis function: cluster size distribution across all runs
# =============================================================================

def analyze_results():
    """
    Parse .species output files from every snapshot across all runs and
    compute the LAT cluster size distribution and largest-cluster histogram.
    """

    print("\n" + "=" * 60)
    print("Step 3: Analyzing cluster size distributions...")
    print("=" * 60)

    # analysis[:,0] = cluster size i
    # analysis[:,1] = accumulated LAT molecule count in size-i clusters
    analysis = np.empty((0, 2))

    # fg: largest cluster per snapshot — proxy for gel/superaggregate formation
    fg = []

    for sim in range(NUM_SIMULATIONS):
        for samp in range(NUM_SAMPLES):
            species_file = f"{RESULTS_DIR}/equil{sim}/run{samp}/{BNGL_FILE}_samp{samp}.species"
            largest_cluster = 0

            with open(species_file, 'r') as f:
                for line in f:
                    if line[0] == '#':
                        continue  # Skip comment/header lines

                    species_name = line.strip().split(' ')[0]

                    # Count LAT molecules in this species by counting 'L'
                    # occurrences in the species string (one 'L' per LAT molecule)
                    num_LAT = species_name.count('L')
                    if num_LAT == 0:
                        continue  # Skip non-LAT species (free Grb2, SOS1, etc.)

                    if num_LAT > largest_cluster:
                        largest_cluster = num_LAT

                    # Accumulate total LAT in size-i aggregates
                    if num_LAT in analysis[:, 0]:
                        analysis[analysis[:, 0] == num_LAT, 1] += num_LAT
                    else:
                        analysis = np.vstack([analysis, [num_LAT, num_LAT]])

            fg.append(largest_cluster)

    # Normalize to fraction of total LAT across all snapshots
    analysis[:, 1] /= (NUM_SIMULATIONS * NUM_SAMPLES * TOTAL_LAT)
    print(f"  Sanity check — sum of fractions: {np.sum(analysis[:, 1]):.4f} (should be ~1.0)")

    # ------------------------------------------------------------------
    # Plot 1: Cluster size distribution
    # ------------------------------------------------------------------
    # Log scale on x-axis: cluster sizes can span several orders of magnitude
    plt.figure()
    plt.scatter(np.log10(analysis[:, 0]), analysis[:, 1], color='k', s=5)
    plt.xlabel(r'$\log_{10}$(LAT aggregate size $i$)')
    plt.ylabel(r'Fraction of LAT in aggregates of size $i$')
    plt.title('LAT Cluster Size Distribution')
    plt.tight_layout()
    plt.savefig('frac_LAT.png', dpi=150)
    plt.close()
    print("  Saved: frac_LAT.png")

    # ------------------------------------------------------------------
    # Plot 2: Largest cluster histogram (gel fraction proxy)
    # ------------------------------------------------------------------
    # A bimodal or heavy-tailed distribution here signals that the system
    # spends time in a gel-like superaggregate phase
    plt.figure()
    mean_fg = np.mean(fg)
    plt.hist(fg, bins=20, density=True, color='steelblue', edgecolor='white', alpha=0.85)
    plt.vlines(x=mean_fg, ymin=0, ymax=plt.gca().get_ylim()[1],
               color='k', zorder=10, label=fr'Avg. $f_g$ = {mean_fg:.1f}')
    plt.xlabel(r'Largest cluster size $f_g$')
    plt.ylabel('Probability density')
    plt.title('Distribution of Largest LAT Cluster per Snapshot')
    plt.legend()
    plt.tight_layout()
    plt.savefig('fg_LAT.png', dpi=150)
    plt.close()
    print("  Saved: fg_LAT.png")


# =============================================================================
# Entry point
# =============================================================================
# The if __name__ == '__main__' guard is required for multiprocessing on
# macOS and Windows, where new processes are spawned (not forked) and would
# otherwise re-execute this top-level code in each worker.

if __name__ == '__main__':

    print("=" * 60)
    print(f"Steps 1–2: Equilibrating and sampling {NUM_SIMULATIONS} runs")
    print(f"           using {NUM_WORKERS} parallel workers...")
    print("=" * 60)

    # Distribute simulation runs across workers. pool.map blocks until all
    # NUM_SIMULATIONS calls to run_simulation() have completed.
    with mp.Pool(processes=NUM_WORKERS) as pool:
        results = pool.map(run_simulation, range(NUM_SIMULATIONS))

    # Print per-run status messages from each worker
    print("\nWorker results:")
    for r in results:
        print(f"  {r}")

    # Run analysis sequentially once all simulation output is on disk
    analyze_results()

    print("\nDone.")