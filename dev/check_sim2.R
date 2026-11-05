# Phase 5 check: rrcsurv::rrc() against Simulation 2's own estimators (sim2_repro, the code behind
# Table 2) at the full Simulation 2 sizes (n1 = 6000, n2 = 2000 with up to 8 pairs each, 60 years).
# For each (cell, replicate), sim2_repro's one_rep() gives the reference fits; the same data are then
# regenerated from the same seed, converted to person-time format and refitted with rrc().
# run from rrcsurv/:  Rscript dev/check_sim2.R [cells] [reps] [cores]
args <- commandArgs(TRUE)
cells <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else c(1, 2, 5, 8, 11)
reps <- if (length(args) >= 2) seq_len(as.integer(args[2])) else 1:2
cores <- if (length(args) >= 3) as.integer(args[3]) else 5

Sys.setenv(SIM2_SCENARIO = "beta1")
source("dev/sim2_data.R")
sim <- globalenv()

one_check <- function(ci, rep) {
  cell <- sim$CELLS[sim$CELLS$cell == ci, ]
  t0 <- proc.time()[3]
  ref <- sim$one_rep(cell, rep)       # sim2_repro's fits (it also fits Cox and ERR on the true x)
  t_old <- proc.time()[3] - t0

  # regenerate the same data from the same seed
  s <- sim2_draw(cell, rep)
  ex <- s$ex; sv <- s$sv; pairs <- s$pairs; pos <- s$pos
  B <- sim$BETA
  dA <- to_rrc(ex, sv, pairs)
  # the log calibration uses only pairs with c, C > 0; the others count as unvalidated years
  dM <- if (all(pos)) dA else to_rrc(ex, sv, pairs[pos, ])
  # naive: every year at its surrogate value
  d0 <- dA$dat; d0$dose <- ex$C[cbind(d0$id, d0$start + 1)]

  fm <- FM2
  t0 <- proc.time()[3]
  new <- list()
  fit1 <- function(d, model, cal, floor = NULL, naive = FALSE, start) {
    f <- suppressWarnings(rrc(fm, data = d$dat, id = "id", validation = d$val, calibration = cal,
                              model = model, min_pairs = sim$MIN_RS_PAIRS, naive = naive,
                              surrogate_floor = floor, start = start,
                              variance = if (naive) "none" else "sandwich"))
    f
  }
  for (m in c("cox", "err")) {
    f <- fit1(list(dat = d0, val = dA$val), m, c ~ dose, naive = TRUE, start = B[[m]])
    new[[paste0("naive_", m)]] <- c(b = unname(f$naive$coef), se = sqrt(f$naive$vcov[1, 1]),
                                    conv = f$naive$converged)
    f <- fit1(dA, m, c ~ dose, start = B[[m]])
    new[[paste0(m, "_A")]] <- c(b = unname(coef(f)), se = sqrt(vcov(f)[1, 1]), conv = f$fit$converged)
    f <- fit1(dM, m, log(c) ~ log(dose), floor = if (cell$error == "add") 1e-8, start = B[[m]])
    new[[paste0(m, "_M")]] <- c(b = unname(coef(f)), se = sqrt(vcov(f)[1, 1]), conv = f$fit$converged)
  }
  t_new <- proc.time()[3] - t0
  do.call(rbind, lapply(names(new), function(nm) data.frame(
    cell = ci, model = cell$model, error = cell$error, rho = cell$rho, rep = rep, fit = nm,
    b_old = ref[[paste0("b_", nm)]], b_new = new[[nm]][["b"]],
    se_old = ref[[paste0("se_", nm)]], se_new = new[[nm]][["se"]],
    conv_old = ref[[paste0("conv_", nm)]], conv_new = as.logical(new[[nm]][["conv"]]),
    t_old = t_old, t_new = t_new)))
}

jobs <- expand.grid(rep = reps, cell = cells)
t0 <- proc.time()[3]
res <- parallel::mclapply(seq_len(nrow(jobs)), function(i) one_check(jobs$cell[i], jobs$rep[i]),
                          mc.cores = cores)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) print(res[bad])
res <- do.call(rbind, res[!bad])
res$d_b <- abs(res$b_new - res$b_old); res$d_se <- abs(res$se_new - res$se_old)
saveRDS(res, "dev/check_sim2.rds")

both <- res$conv_old & res$conv_new
cat(sprintf("%d fits (%d cells x %d reps x 6), %.1f min\n", nrow(res), length(cells), length(reps),
            (proc.time()[3] - t0) / 60))
cat(sprintf("converged in both: %d; old only: %d; new only: %d; neither: %d\n", sum(both),
            sum(res$conv_old & !res$conv_new), sum(!res$conv_old & res$conv_new),
            sum(!res$conv_old & !res$conv_new)))
cat(sprintf("max |b_new - b_old| (both converged): %.2e;  max |se_new - se_old|: %.2e\n",
            max(res$d_b[both]), max(res$d_se[both], na.rm = TRUE)))
print(aggregate(cbind(d_b, d_se) ~ fit, res[both, ], max), digits = 3)
if (any(!both)) {
  cat("\nfits not converged in at least one implementation:\n")
  print(res[!both, c("cell", "model", "error", "rho", "rep", "fit", "b_old", "b_new", "conv_old", "conv_new")],
        digits = 5, row.names = FALSE)
}
tm <- unique(res[, c("cell", "rep", "t_old", "t_new")])
cat(sprintf("\ntime per replicate (6 fits): old %.0f s, rrcsurv %.0f s (median)\n",
            median(tm$t_old), median(tm$t_new)))
