# implements method of composition:
## (1) perform nsim regressions (beta regression or quasibinomial glm, depending on specified type)
## (2) sample from resulting distributions of regression coefficents
sample_coefs <- function(stmobj, formula, type, metadata, nsims = 25, seed = NULL) {
  response <- as.character(formula)[2]
  topic_n <- eval(parse(text=response))
  topic_nam <- paste0("Topic", topic_n)
  f <- paste(topic_nam, "~", as.character(formula)[3])
  res <- list()
  set.seed(seed)
  formals(glm)$family <- quasibinomial(link = "logit")
  fit_reg <- if (type == "beta") betareg::betareg else if (type == "quasibinomial") glm
  for (k in 1:length(topic_n)) {
    theta_sim <- do.call(rbind,
                         stm::thetaPosterior(stmobj, nsims = nsims, type = "Global"))[,topic_n[k]]
    theta_sim <- lapply(split(1:(NROW(theta_sim)), 1:nsims),
                        function(i) setNames(data.frame(theta_sim[i]), topic_nam[k]))
    res[[k]] <- lapply(theta_sim, function(x) {
      res_tmp <- fit_reg(formula = as.formula(f[k]), data = cbind(x, metadata))
      sample_normal(res_tmp, type)}
    )
  }
  return(setNames(res, topic_nam))
}

# obtain theta proportions and credible intervals (given set of previously sampled coefficients)
## for full range of variable est_var, holding all variables but est_var as median/majority
predict_props <- function(beta_coefs, est_var, formula,
                          metadata, ci_lower = 0.025, ci_upper = 0.975) {
  response <- as.character(formula)[2]
  topic_n <- eval(parse(text=response))
  topic_nam <- paste0("Topic", topic_n)
  if(is.numeric(metadata[,est_var])) {
    range_est_var <- seq(min(metadata[,est_var]), max(metadata[,est_var]), length.out = 500)
  } else {
    range_est_var <- unique(metadata[,est_var])
  }
  xmat <- make_median_xmat(est_var, formula, metadata, range_est_var)
  res <- list()
  for (k in topic_nam){
    fit_vals <- do.call(cbind, lapply(beta_coefs[[k]], function(x) sigmoid(xmat %*% t(x))))
    mu <- quanteda::rowMeans(fit_vals)
    ci <- apply(fit_vals, 1, function(x) quantile(x, probs = c(ci_lower, ci_upper)))
    res[[k]] <- data.frame(range_est_var, mu, ci[1,], ci[2,])
    names(res[[k]]) <- c(est_var, "proportion", "ci_lower", "ci_upper")
  }
  return(setNames(res, topic_nam))
}
