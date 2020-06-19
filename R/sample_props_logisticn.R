# sample mean and credible interval from LogisticNormal for a given model (stmobj)
## and full observed range of variable (est_var),
## holding all variables but est_var as median/majority
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
