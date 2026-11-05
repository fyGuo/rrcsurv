#### Risk sets from counting-process (person-time) data ####

#' Build the risk-set structure of counting-process data
#'
#' Each row covers the interval (start, stop]. The distinct event times are the `stop` values of rows
#' with `event == 1`; the risk set at event time t holds every row with start < t <= stop (the row
#' that describes its subject at t), which allows delayed entry and gaps. Ties use the Breslow
#' convention: every case at t shares the full risk set at t.
#'
#' @param start,stop numeric interval limits, `start < stop`.
#' @param event 0/1 event indicator; an event happens at `stop`.
#' @return A list of class `rrc_rs`:
#'   \describe{
#'     \item{time}{sorted distinct event times}
#'     \item{d}{number of cases at each event time}
#'     \item{row}{for each (event time, at-risk row) pair, the row index in the data}
#'     \item{uid}{for each pair, the index of its event time in `time`}
#'     \item{case}{for each case, the index of its own pair in `row`/`uid`}
#'   }
#' @keywords internal
rs_build <- function(start, stop, event) {
  n <- length(stop)
  if (length(start) != n || length(event) != n) stop("start, stop and event must have the same length")
  if (anyNA(start) || anyNA(stop) || anyNA(event)) stop("start, stop and event must not be missing")
  if (any(start >= stop)) stop("every row needs start < stop")
  if (!all(event %in% c(0, 1))) stop("event must be 0/1")
  ev <- which(event == 1)
  if (length(ev) == 0) stop("no events")

  time <- sort(unique(stop[ev]))
  d <- tabulate(match(stop[ev], time), length(time))
  # event times covered by row r are time[lo:hi]: start < time <= stop
  lo <- findInterval(start, time) + 1L
  hi <- findInterval(stop, time)
  len <- pmax(hi - lo + 1L, 0L)
  row <- rep.int(seq_len(n), len)
  uid <- sequence(len, from = lo)
  o <- order(uid, row)
  row <- row[o]; uid <- uid[o]

  # a case's own pair: its row at its event time (its stop)
  cu <- match(stop[ev], time)
  key <- (uid - 1) * n + row
  case <- match((cu - 1) * n + ev, key)
  structure(list(time = time, d = d, row = row, uid = uid, case = case, n_rows = n),
            class = "rrc_rs")
}
