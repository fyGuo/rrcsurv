#### Example data ####

#' Example cohort with an error-prone cumulative exposure
#'
#' Simulated data for the examples and the vignette: a cohort of 1500 workers followed for up to
#' 25 years, with a yearly dose recorded with multiplicative error, and an internal validation
#' study in which the true dose of 1 to 5 years of 400 workers was measured. Built by
#' `data-raw/make_example.R` in the package source.
#'
#' The recorded (surrogate) dose C comes first, lognormal and correlated within a worker; the true
#' dose c follows the measurement error model
#' `log c = -0.3 + 0.6 log C - 0.2 smoke + e`, `e ~ N(0, 0.4^2)`, independent across years. The
#' hazard is `0.004 exp(0.25 x(t) + 0.5 smoke)`, where x(t) is the true dose received before year
#' t starts; events are recorded at the end of a year. True values: 0.25 for the exposure and 0.5
#' for smoking (Cox model).
#'
#' @format `rrc_main`: the main study in person-time format, one row per worker and year of
#'   follow-up (22,947 rows, 482 events):
#'   \describe{
#'     \item{id}{worker.}
#'     \item{start, stop}{the year, (start, stop].}
#'     \item{event}{1 if the worker's event occurred in this year.}
#'     \item{dose}{the dose of the year as known to the analyst: the surrogate C, or the true
#'       dose c in a validated year.}
#'     \item{surrogate}{`TRUE` if `dose` is the surrogate, `FALSE` in validated years.}
#'     \item{smoke}{smoker (1) or not (0).}
#'     \item{dose_true}{the true dose c of every year; known only because the data are simulated.}
#'   }
#' @source Simulated.
#' @seealso [rrc()], [cumfit()].
"rrc_main"

#' @rdname rrc_main
#' @format `rrc_validation`: the validation study, one row per validated worker-year (1,159
#'   pairs from 400 workers):
#'   \describe{
#'     \item{id}{worker, as in `rrc_main`.}
#'     \item{year}{when the dose was received: the end of its year, the `stop` of its row in
#'       `rrc_main` (use as `dose_time`).}
#'     \item{c}{the true dose.}
#'     \item{dose}{the surrogate dose C, under the exposure's name in `rrc_main`.}
#'     \item{smoke}{as in `rrc_main`.}
#'   }
"rrc_validation"
