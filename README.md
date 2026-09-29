# KSHV_math
Mathematical analysis for replication and segregation of KSHV

## Getting started

Open `KSHV_mathematical_analysis.Rproj` in RStudio. All scripts locate files with
the `here` package relative to the project root, so they should be run from the
project (not from inside `scripts/`).

Scripts write to a `results/` directory that is not tracked in git. 

### Environment

Package versions are pinned with [renv](https://rstudio.github.io/renv/). The
project `.Rprofile` activates renv automatically when the project is opened; to
install the pinned versions:

```r
renv::restore()
```

The lockfile targets R 4.4.0. Note that `rstan` compiles from source and needs a
working C++ toolchain, so the first restore can take a while.

For reference, the packages used directly by the analysis scripts are:

`bayesplot`, `cluster`, `cowplot`, `doParallel`, `factoextra`, `fitdistrplus`,
`foreach`, `furrr`, `ggbeeswarm`, `ggdist`, `ggExtra`, `ggh4x`, `ggnewscale`,
`ggrepel`, `ggthemes`, `here`, `Hmisc`, `MASS`, `patchwork`, `purrr`,
`RColorBrewer`, `rstan`, `scales`, `scico`, `tidyverse`, `tune`, `ggh4x`, and `ggbeeswarm` 

`rstan` is used only for convergence diagnostics (`Rhat`, `ess_bulk`,
`ess_tail`) in `functions_inference.R` and `functions_run_pipeline.R`.


## Data

All data needed to reproduce the analysis are in `data/derived/`. 

### Imaging data (SUM159 cells)

Six Excel files, one dividing/non-dividing pair per experimental condition.
These are read by `load_data()` in `functions_run_pipeline.R`.

| File | Rows | Mother cells |
| --- | --- | --- |
| `fixed_8TR_dividing_cells.xlsx` | 91 | 40 |
| `fixed_8TR_non_dividing_cells.xlsx` | 98 | — |
| `fixed_KSHV_dividing_cells.xlsx` | 36 | 33 |
| `fixed_KSHV_non_dividing_cells.xlsx` | 55 | — |
| `live_KSHV_dividing_cells.xlsx` | 28 | 23 |
| `live_KSHV_non_dividing_cells.xlsx` | 31 | — |

**Dividing-cell files** (`*_dividing_cells.xlsx`) supply the daughter-cell data.
One row per LANA cluster, with the two daughters of each division stored
side by side:

- `mother_cell_id` — identifier for the division event, linking the daughter pair
- `Cluster_daughter1` / `Cluster_daughter2` — cluster index within each daughter
- `Min episome in cluster_daughter1` / `_daughter2` — minimum number of episomes
  the cluster is known to contain
- `Total cluster intensity_daughter1` / `_daughter2` — summed LANA fluorescence
  intensity for the cluster, the observable used to infer episome number

A mother cell with a cluster count in only one daughter has `NA` in the other
daughter's columns; these are dropped during loading. 

**Non-dividing-cell files** (`*_non_dividing_cells.xlsx`) supply the
data used to inform the distribution of episomes per mother cell prior to division:

- `Image`, `Cell` — identifiers, combined into `cell_id` during loading
- `Cluster` — cluster index within the cell (`0` marks a cell with no clusters)
- `Min episome in cluster` — minimum episomes in the cluster
- `Total cluster intensity` — summed LANA fluorescence intensity

### Longitudinal data (Brk.219 cells)

Two CSV files used only by `simulate_BRK219_experiments.R`.

`brk219_full_LANA_dots.csv` — 547 rows, one per imaged nucleus.

- `day` — days post-infection; measurements at days 0, 2, 4, …, 20
- `LANA_dots` — number of LANA dots counted in that nucleus (range 0–58)

`brk219_cell_growth.csv` — 33 rows, one per counting timepoint.

- `day` — day of the passaging experiment, 0 through 45
- `live_cells` — live cell count, normalized to 1 at day 0
- `dead_cells` — dead cell count; `NA` for most timepoints (13 of 33 recorded)


## Estimating replication and segregation efficiency
To calculate estimates of replication and segregation efficiency, run the following scripts for each experimental condition:

- **analyze_fixed_8TR.R** for estimates from the fixed images of 8TR cells 
  - results in Figure 2, Figure S4A, Figure S5A, Figure S15A
- **analyze_fixed_KSHV.R** for estimates from the fixed images of KSHV cells 
  - results in Figure 5, Figure S4C, Figure S5B, Figure S15B
- **analyze_live_KSHV.R** for estimates from the images of live KSHV cells 
  - results in Figure 6, Figure S4D, Figure S5C, Figure S15C&D

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

`run_pipeline()` caches expensive intermediate results as `.RData` files in the
results folder (`MCMC_samples_per_cluster.RData`, `MCMC_samples_per_cell.RData`,
`MCMC_convergence.RData`, `MLE_with_uncertainty.RData`,
`MLE_Pr_with_uncertainty.RData`) and reloads them on subsequent runs. Run with `overwrite = TRUE` to force a re-run.

Supplemental figures S4, S5, and S14 are generated with the **generate_supplemental_figures.R** script.

Histograms for Figures 1 and 4 are generated with **generate_fig1_4_histograms.R**.
This script has no external inputs — its data are hardcoded — and writes directly
to `results/`.

## Simulations

Simulations for four cell-growth scenarios can be run with the following scripts: 

- **simulate_constant.R** simulates a constant-sized cell population without selection 
  - Figures 7, S8, S9, S10
- **simulate_constant_selection.R** simulates a constant-sized cell population under selection
  - Results not included in manuscript
- **simulate_expo_selection.R** simulates an exponentially-growing cell population of immortal cells under selection
  - Figure 8
- **simulate_PEL_growth.R** simulates a KSHV-dependent tumor under therapy that reduces replication or segregation efficiency
  - Figures 9, S12, S13

Functions for each of these scenarios are defined in **functions_simulations.R**, along with a few plotting functions. The main simulation functions are:

- `makeChildren()` Simulates the fate of episomes during cell division according to a given replication and segregation efficiency
- `simStepFlex()` Simulates one step of the cell population dynamics by sampling the type of event (either cell birth or death) and the time advance
- `extinction()` Simulates a dividing cell population that grows until it reaches a designated size, around which it fluctuates. If there is no selection against cells without episomes, simulations end when there are no more episomes in the population or when the simulation reaches a designated stop time. If there is selection against cells without episomes, the simulation is run for 700 generations. The distribution of episomes in the population is recorded at each time step.
- `exponential_growth()` Simulates a dividing cell population that grows exponentially until it reaches a designated size. The distribution of episomes in the population is recorded at designated population sizes.
- `PEL_simulations()` Uses the `exponential_growth()` function to simulate a tumor that grows until it reaches a designated size with a baseline replication and segregation efficiency. Then, replication and/or segregation efficiency is reduced and the tumor is simulated until it reaches a designated size, dies off, or until a specified time.

## Validation of predicted longitudinal episome dynamics
We evaluated the consistency of model parameters informed by images in SUM159 cells using experiments tracking longitudinal episome dynamics in Brk.219 cells. The interpretation of these longitudinal LANA dot measurements and comparison to SUM159-informed model predictions (Figure S11) is contained in **simulate_BRK219_experiments.R** simulates. It relies on the simulation functions described above.

The script fits piecewise exponential growth rates to the cell counts, converts
LANA-dot measurements from days to generations, and fits the decay in dot number
per generation with a negative binomial regression (`MASS::glm.nb`). Replication
efficiency is recovered from the fitted slope as `Pr = 1 + b`, and the fitted
intercept gives the mean episome count at day 0 used to seed the simulations. 

Three local helpers are defined in this script rather than in the shared function
files:

- `simulate_cell_growth_variable()` simulates cell growth with a separate rate per
  passaging interval
- `sample_initial_epi()` draws per-cell episome counts from a negative binomial and
  bins them into the initial-conditions vector the simulator expects
- `sim_passage_wrapper_vary_b()` runs the passaging simulation with a growth rate
  that varies by interval

### Run order dependency

**`analyze_fixed_KSHV.R` must be run before `simulate_BRK219_experiments.R`.**

The BRK219 script reads the SUM159 confidence intervals produced by the fixed
KSHV analysis:

```r
SUM159_confidence_intervals <- read_csv(here("results", "fixed_KSHV", "MLE_parameter_estimates.csv"))
```

That file is written by `analyze_fixed_KSHV.R` at line 86. Because `results/` is
gitignored, it will not exist in a fresh clone, and the BRK219 script will fail
at this read with a missing-file error. `Pr_CI` and `Ps_CI` are parsed from it
and used to bracket the simulated confidence bounds.

The BRK219 script also sources `functions_simulations.R`, `functions_inference.R`,
and `functions_run_pipeline.R`, and reads both `brk219_*.csv` data files directly.


## Benchmarking methods

Bias of the ML estimates can be assessed using synthetic data by running the **benchmark_MLE.R** script (Figures S16, S17). Synthetic data are generated with `simulate_multiple_cells()`, which simulates one division for all cells in a specified population size.

The sensitivity of parameter estimates from MCMC to the choice of prior for $n_k$ can be assessed by running the **benchmark_prior_sensitivity.R** script (Figure S14).

## Output directories

`results/` is gitignored and starts empty. These scripts create their output
subdirectory on startup if it does not already exist:

| Script | Directory |
| --- | --- |
| `analyze_fixed_8TR.R`, `analyze_fixed_KSHV.R`, `analyze_live_KSHV.R` | `results/fixed_8TR/`, `results/fixed_KSHV/`, `results/live_KSHV/` (via `run_pipeline()`) |
| `generate_supplemental_figures.R` | `results/supplemental_figures/` |
| `generate_fig1_4_histograms.R` | `results/` |
| `simulate_constant.R` | `results/simulations_constant/` |
| `simulate_constant_selection.R` | `results/simulations_constant_selection/` |
| `simulate_expo_selection.R` | `results/simulations_expo_selection/` |
| `simulate_PEL_growth.R` | `results/PEL_simulations_with_selection/` |
| `simulate_BRK219_experiments.R` | `results/brk219/` |
| `benchmark_prior_sensitivity.R` | `results/supplemental_figures/`, `results/supplemental_figures_updated_pdf_n_prior/` |
| `benchmark_MLE.R` | `results/benchmarking/` |



