# Development checks

Scripts used to validate rrcsurv against the reference implementations of the paper's simulations.
They are not part of the package build (`.Rbuildignore`) and need the analysis repository
[measurementERROR_COX_ERR_RCC](https://github.com/fyGuo/measurementERROR_COX_ERR_RCC) cloned next to
`rrcsurv/` (its `sim1_repro/` and `sim2_repro/` folders), so they run only in that layout, from `rrcsurv/`.

| Script | What it checks | Result |
|---|---|---|
| `check_engine_vs_sim1.R` | Cox / ERR engine against Simulation 1's partial likelihood code | |
| `make_sim1_fixture.R` | builds `tests/testthat/fixtures/sim1_cells.rds` (Simulation 1 references) | |
| `check_variance_mc.R` | Monte Carlo coverage of the sandwich variance, 500 reps x 4 cells | `check_variance_mc.log` |
| `sim2_data.R` | shared loader and converter of Simulation 2 data to rrc() input | |
| `check_sim2.R` | rrc() against sim2_repro's fits, 5 cells x 3 reps x 6 fits at full size | `check_sim2.log` |
| `check_sim2_mc.R` | Monte Carlo at full Simulation 2 size, 300 reps x 2 cells | `check_sim2_mc.log` |
| `make_sim2_fixture.R` | builds `tests/testthat/fixtures/sim2_cells.rds` (Simulation 2 references) | |
