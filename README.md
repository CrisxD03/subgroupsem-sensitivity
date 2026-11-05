# SubgroupSEM Sensitivity Study

This repository contains the code and results for the sensitivity study in our replication report of Kiefer et al. (2022), "Subgroup Discovery in Structural Equation Models."

The paper shows SubgroupSEM recovering planted subgroups from one synthetic dataset with one sample size, one effect size and one interestingness measure. We test what happens when you change those.

## What we test

Four hypotheses:

- **H1 (effect size):** Does recovery hold when the planted effects are smaller?
- **H2 (sample size):** Does recovery hold at smaller or larger populations?
- **H3 (interestingness measure):** Do LR, Wald and user-defined agree on which groups are recovered?
- **H4 (SEM tree):** Does SubgroupSEM recover at least as many groups as the SEM tree?

We run 16 conditions (4 effect sizes x 4 sample sizes), 20 datasets each, 320 datasets total.

## How to reproduce

### Requirements

- R 4.4.0+
- CRAN packages: lavaan (0.7-2), OpenMx (2.21.13), semtree (0.9.23), reticulate, withr, jsonlite, ggplot2, testthat (for the tests)
- subgroupsem from GitHub at commit fd24703:
  ```r
  remotes::install_github("chkiefer/subgroupsem@fd24703")
  ```
- Python 3.11 with pysubgroup 0.7.8 and setuptools < 70. You don't install these yourself: the code calls `reticulate::py_require()`, which downloads them on first use.

### Steps

Run everything from the project root.

1. Run the full study (takes several hours, writes one .rds per dataset to `results/sensitivity/runs/`):
   ```
   Rscript scripts/20_run_grid.R
   ```
   This uses 6 parallel workers. You can change that in the `settings` list at the top of `scripts/20_run_grid.R`. Finished datasets are skipped on rerun, so you can stop and resume. `Rscript scripts/20_run_grid.R timing` runs one dataset per sample size, to see how long a run takes on your machine.
2. Summarise and evaluate hypotheses:
   ```
   Rscript scripts/30_summarise.R
   ```
   This produces `results/sensitivity/summary/by_condition.csv`, `runs.csv`, `hypotheses.txt` and the figures. It needs the .rds files from step 1 and stops if any dataset is missing or failed.

To run the unit tests:
```
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### What produces what

| Report table | Data source | Script |
|---|---|---|
| Per-group recovery | `summary/by_condition.csv`, rows "SubgroupSEM (LR)", columns `rate_G1` to `rate_G3` (times 20) and `full_k` for "All 3" | `30_summarise.R` |
| H3 disagreements | `summary/by_condition.csv`, columns `rate_G1` to `rate_G3` for LR, Wald and UD (times 20). The disagreeing pairs are listed under H3 in `summary/hypotheses.txt` | `30_summarise.R` |
| H4: LR vs SEM tree | `summary/by_condition.csv`, column `mean_recovered` for "SubgroupSEM (LR)" and "SEM tree" | `30_summarise.R` |

### Pre-computed results

If you do not want to rerun the study, the `results/sensitivity/summary/` folder has the aggregated CSVs. The raw .rds files are not included. The full run was done on R 4.4.0, macOS Sonoma 14.6.1 (see `results/sensitivity/sessionInfo.txt`). The settings of that run are in `results/sensitivity/settings.json` and every dataset with its seed is in `results/sensitivity/design.csv`.

## Project structure

```
R/                    Core functions (simulate, measure, score, search, tree)
scripts/              Runner and summariser scripts
tests/testthat/       Unit tests for the core functions
results/sensitivity/  Pre-computed summary results
```

## Seeds

Every dataset has its own seed: 2022 + 1000 * i_N + r, where i_N is the position of the sample size in (500, 2500, 5000, 10000), so 1 to 4, and r is the run number (1-20). For example N = 500, r = 1 uses seed 3023. The seed does not depend on the effect size, so the four effect sizes of one run share the same random draws.

The base model is the Introductory Example from the paper (one factor, five indicators, parameters from Table 1).

## Reference

Kiefer, C., Lemmerich, F., Langenberg, B., Mayer, A. (2022). Subgroup discovery in structural equation models. Psychological Methods. https://doi.org/10.1037/met0000524
