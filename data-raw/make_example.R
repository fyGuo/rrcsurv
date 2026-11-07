#### Build the example data sets rrc_main and rrc_validation ####
# run from rrcsurv/:  Rscript data-raw/make_example.R
#
# A cohort of 1500 workers followed for up to 25 years. Each year a worker receives a dose c; the
# exposure records hold only a surrogate C. As in the paper's simulations, the surrogate is generated
# first (lognormal, correlated within a worker) and the true dose from it, with multiplicative error
# that depends on smoking:
#   log c = -0.3 + 0.6 log C - 0.2 smoke + e,  e ~ N(0, 0.4^2), independent across years.
# The hazard depends on the true cumulative dose of the previous years, x(t), and on smoking:
#   Cox: lambda(t) = 0.004 exp(0.25 x(t) + 0.5 smoke), with events at the end of a year.
# An internal validation study measured the true dose of 1 to 5 of the years of 400 workers. Those
# years carry the true dose in the main study (surrogate = FALSE).

set.seed(20261107)
n <- 1500; T <- 25; n_val <- 400
smoke <- rbinom(n, 1, 0.3)
level <- rnorm(n, log(0.4), 0.5)                       # worker's typical log surrogate dose
drop_out <- pmin(T, 1 + rgeom(n, 0.03))                # last year under follow-up if no event

rows <- vector("list", n)
for (i in seq_len(n)) {
  C <- exp(level[i] + rnorm(T, 0, 0.4))
  c_true <- exp(-0.3 + 0.6 * log(C) - 0.2 * smoke[i] + rnorm(T, 0, 0.4))
  x_prev <- c(0, cumsum(c_true))[1:T]                  # dose of years 1..j-1 at year j
  h <- 0.004 * exp(0.25 * x_prev + 0.5 * smoke[i])
  ev <- rbinom(T, 1, 1 - exp(-h))
  last <- min(drop_out[i], which(ev == 1)[1], na.rm = TRUE)
  j <- seq_len(last)
  rows[[i]] <- data.frame(id = i, start = j - 1, stop = j, event = as.integer(j == last & ev[last] == 1),
                          dose = signif(C[j], 4), dose_true = signif(c_true[j], 4), smoke = smoke[i])
}
main <- do.call(rbind, rows)

# validation: 400 workers, 1 to 5 of their follow-up years each
vid <- sample(n, n_val)
val <- do.call(rbind, lapply(vid, function(i) {
  yrs <- main$stop[main$id == i]
  k <- min(length(yrs), sample(1:5, 1))
  data.frame(id = i, year = sort(yrs[sample.int(length(yrs), k)]))
}))
key_m <- paste(main$id, main$stop); key_v <- paste(val$id, val$year)
m <- match(key_v, key_m)
rrc_validation <- data.frame(id = val$id, year = val$year, c = main$dose_true[m], dose = main$dose[m],
                             smoke = main$smoke[m])
main$surrogate <- !(key_m %in% key_v)
main$dose[!main$surrogate] <- main$dose_true[!main$surrogate]
rrc_main <- main[, c("id", "start", "stop", "event", "dose", "surrogate", "smoke", "dose_true")]
rownames(rrc_main) <- rownames(rrc_validation) <- NULL

cat(sprintf("%d rows, %d subjects, %d events; %d validation pairs from %d subjects\n", nrow(rrc_main),
            n, sum(rrc_main$event), nrow(rrc_validation), length(unique(rrc_validation$id))))
save(rrc_main, file = "data/rrc_main.rda", compress = "xz")
save(rrc_validation, file = "data/rrc_validation.rda", compress = "xz")
