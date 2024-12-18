# KSHV_mathematical_analysis_tidy
Mathematical analysis for replication and segregation of KSHV

## Estimating replication and segregation efficiency
To calculate estimates of replication and segregation efficiency, run the following scripts for each experimental condition:

- **analyze_fixed_8TR.R** for estimates from the fixed images of 8TR cells 
  - results in Figure 2, Figure S4A, Figure S5A, Figure S14A
- **analyze_fixed_KSHV.R** for estimates from the fixed images of KSHV cells 
  - results in Figure 5, Figure S4C, Figure S5B, Figure S14B
- **analyze_live_KSHV.R** for estimates from the images of live KSHV cells 
  - results in Figure 6, Figure S4D, Figure S5C, Figure S14C&D

These estimates rely on functions in the following scripts.

**functions_run_pipeline.R** Contains four functions:

- `load_data()` to read in and format the data
- `run_pipeline()` to estimate the number of episomes per cell via Gibbs sampling
and find ML estimates of replication and segregation efficiencies
- `make_plots()` to make diagnostic plots from the analysis outputs, some of which are included in the supplement
- `figures()` to make the figure included in the main text

**functions_inference.R** Includes functions called from `run_pipeline()` and `make_plots()` required 
to implement Gibbs sampling, compute likelihoods, and quantify uncertainty. The main functions are:

- `likelihood()` Computes the likelihood of observing an observed daughter cell pair given $X_0$, $P_r$, and $P_s$
- `calculate_maximum_likelihood_unknownX0()` runs a grid search to find the maximum likelihood estimates of Pr and Ps marginalized over values of $X_0$
- `calculate_CI()` calculates the joint 95% confidence interval based on the results of `calculate_maximum_likelihood_unknownX0()`
- `run_grid_search()` combines the prior two functions and plots the results of the grid search
- `run_gibbs()` implements Gibbs sampling to estimate the number of episomes per LANA dot
- `convergence_results()` calculate convergence heuristics for Gibbs sampling

Supplemental figures S4, S5, and S14 are generated with the **generate_supplemental_figures.R** script.


## Simulations

Simulations for four cell-growth scenarios can be run with the following scripts: 

- **simulate_constant.R** simulates a constant-sized cell population without selection 
  - Figures 7, S8, S9, S10
- **simulate_constant_selection.R** simulates a constant-sized cell population under selection
  - Results not included in manuscript
- **simulate_expo_selection.R** simulates an exponentially-growing cell population of immortal cells under selection
  - Figure 8
- **simulate_PEL_growth.R** simulates a KSHV-dependent tumor under therapy that reduces replication or segregation efficiency
  - Figures 9, S11, S12

Functions for each of these scenarios are defined in **functions_simulations.R**, along with a few plotting functions. The main simulation functions are:

- `makeChildren()` Simulates the fate of episomes during cell division according to a given replication and segregation efficiency
- `simStepFlex()` Simulates one step of the cell population dynamics by sampling the type of event (either cell birth or death) and the time advance
- `extinction()` Simulates a dividing cell population that grows until it reaches a designated size, around which it fluctuates. If there is no selection against cells without episomes, simulations end when there are no more episomes in the population or when the simulation reaches a designated stop time. If there is selection against cells without episomes, the simulation is run for 700 generations. The distribution of episomes in the population is recorded at each time step.
- `exponential_growth()` Simulates a dividing cell population that grows exponentially until it reaches a designated size. The distribution of episomes in the population is recorded at designated population sizes.
- `PEL_simulations()` Uses the `exponential_growth()` function to simulate a tumor that grows until it reaches a designated size with a baseline replication and segregation efficiency. Then, replication and/or segregation efficiency is reduced and the tumor is simulated until it reaches a designated size, dies off, or until a specified time.

## Benchmarking methods

Bias of the ML estimates can be assessed using synthetic data by running the **benchmark_MLE.R** script (Figures S15, S16). Synthetic data is generated with `simulate_multiple_cells()`, which simulates one division for all cells in a specified population size.

The sensitivity of parameter estimates from MCMC to the choice of prior for $n_k$ can be assessed by running the **benchmark_prior_sensitivity.R** script (Figure S13).


