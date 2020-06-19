#' Sample from LogisticNormal over full range of variable
#'
#' Sample mean and credible interval from LogisticNormal for a given model (stmobj)
#' and full observed range of variable (est_var),
#' holding all variables but est_var as median/majority.
#'
#' @param stmobj Fitted stm model.
#' @param est_var Variable for which to sample over full observed range.
#' @param formula Formula object (the prevalence specification used to fit the stm).
#' @param metadata Metadata that was used to fit the stm.
#' @param nsims Number of draws from the LogisticNormal used to obtain empirical mean and quantiles.
#' @param ci_lower Lower bound of credible interval.
#' @param ci_upper Upper bound of credible interval.
#' @param seed Seed.
#' @export
sample_props_logisticn <- function(stmobj, est_var, formula, metadata, nsims = 1000,
                                   ci_lower = 0.025, ci_upper = 0.975, seed = NULL) {
  gamma <- stmobj$mu$gamma
  Sigma <- stmobj$sigma
  if(is.numeric(metadata[,est_var])) {
    range_est_var <- seq(min(metadata[,est_var]), max(metadata[,est_var]), length.out = 500)
  } else {
    range_est_var <- unique(metadata[,est_var])
  }
  xmat <- make_median_xmat(est_var, formula, metadata, range_est_var)
  mu <- xmat %*% gamma
  set.seed(seed)
  est <- sim_theta(mu, Sigma, nsims = nsims, ci_lower = ci_lower, ci_upper = ci_upper)
  res <- lapply(est, function(x) setNames(cbind(range_est_var, x), c(est_var, names(x))))
  return(res)
}
