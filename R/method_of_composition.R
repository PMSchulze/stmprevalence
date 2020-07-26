#' Perform Method of Composition employing either a quasibinomial GLM or beta regression
#'
#' For a fitted stm, first sample topic proportions from the approximate prosterior,
#' then perform a regression of these proportions on prevalence covariates,
#' and lastly sample from the resulting distribution of regression coefficents.
#' The obtained values of this process, which is repeated nsims times,
#' are samples of the marginal posterior of regression coefficients.
#' This procedure is known as the method of composition in the social sciences.
#'
#' @param stmobj Fitted stm model.
#' @param formula Formula (subset of the prevalence specification used to fit the stm).
#' @param type Regression to perform: Either 'beta' or 'quasibinomial'.
#' @param metadata Metadata that was used to fit the stm.
#' @param nsims Number of repetitions.
#' @param seed Seed.
#' @return A list of lists: For each topic, nsims regression outputs are returned.
#' @export
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


#' Generate mean and credible intervals of topic proportions over full range of a specified variable
#'
#' Given the samples of coefficients obtained using sample_coefs, specify one variable and predict
#' mean and credible intervals (with respect to mean prediction) of topic proportions
#' over the full observed range of variable est_var,
#' while holding all other variables as median/majority. This function is typically used to visualize
#' the results of the method of composition.
#'
#' @param beta_coefs sampled coefficients from method of composition using sample_coefs
#' @param est_var Variable for which to sample over full observed range.
#' @param formula Formula (subset of the prevalence specification used to fit the stm).
#' @param metadata Metadata that was used to fit the stm.
#' @param ci_lower Lower bound of credible interval.
#' @param ci_upper Upper bound of credible interval.
#' @return A list of dataframes: For each topic, the empirical mean and credible intervals are returned.
#' @export
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


#' Apply fully Bayesian version of method of composition
#'
#' For a fitted stm, first sample topic proportions from the approximate prosterior,
#' and then perform a Bayesian beta regression of these proportions on prevalence covariates.
#' The procedure is repeated nsims times.
#'
#' @param stmobj Fitted stm model.
#' @param formula Formula (subset of the prevalence specification used to fit the stm).
#' @param metadata Metadata that was used to fit the stm.
#' @param nsims Number of repetitions.
#' @param seed Seed.
#' @return A list of lists: For each topic, nsims regression outputs are returned.
#' @export
beta_bayes <- function(stmobj, formula, metadata, nsims = 100, seed = 123) {
  response <- as.character(formula)[2]
  topic_n <- eval(parse(text=response))
  topic_nam <- paste0("Topic", topic_n)
  f <- paste(topic_nam, "~", as.character(formula)[3])
  res <- list()
  set.seed(seed)
  for (k in 1:length(topic_n)) {
    theta_sim <- do.call(rbind,
                         stm::thetaPosterior(stmobj, nsims = nsims, type = "Global"))[,topic_n[k]]
    theta_sim <- lapply(split(1:(NROW(theta_sim)), 1:nsims),
                        function(i) setNames(data.frame(theta_sim[i]), topic_nam[k]))
    res[[k]] <- lapply(theta_sim, function(x) {
      rstanarm::stan_betareg(
        formula = as.formula(f[k]),
        link = "logit", link.phi = "log",
        data = cbind(x, metadata),
        algorithm = "optimizing", seed = seed
      )})
  }
  return(setNames(res, topic_nam))
}


#' Generate mean and credible intervals of topic proportions over full range of a specified variable
#'
#' Given the MAP estimates obtained using beta_bayes, specify one variable and predict
#' mean and credible intervals of topic proportions over the full observed range of variable est_var,
#' while holding all other variables as median/majority. This function is typically used to visualize
#' the results of the method of composition.
#'
#' @param bayes_out Output of method of composition with Bayesian regression using beta_bayes
#' @param est_var Variable for which to sample over full observed range.
#' @param formula Formula (subset of the prevalence specification used to fit the stm).
#' @param metadata Metadata that was used to fit the stm.
#' @param ci_lower Lower bound of credible interval.
#' @param ci_upper Upper bound of credible interval.
#' @return A list of dataframes: For each topic, the empirical mean and credible intervals are returned.
#' @export
posterior_predict_props <- function(bayes_out, est_var, formula, metadata, ci_lower, ci_upper) {
  response <- as.character(formula)[2]
  topic_n <- eval(parse(text=response))
  topic_nam <- paste0("Topic", topic_n)
  # --------------------------------------------------------------------------------------------
  if(is.numeric(metadata[ ,est_var])) {
    range_est_var <- seq(min(metadata[,est_var]), max(metadata[,est_var]), length.out = 500)
  } else {
    range_est_var <- unique(metadata[,est_var])
  }
  dat_new <- lapply(metadata[, -which(names(metadata) == est_var)],
                    function(x) if(is.numeric(x)) median(x) else majority(x))
  xmat <- data.frame(dat_new, range_est_var)
  names(xmat) <- c(names(dat_new),est_var)
  levels(xmat$Partei) <- levels(metadata$Partei)
  levels(xmat$Bundesland) <- levels(metadata$Bundesland)
  # ----------------------------------------------------------------------------------------------
  nm <- c(est_var, "proportion", "ci_lower", "ci_upper")
  res <- list()
  for (k in topic_nam){
    preds <- do.call(rbind, lapply(bayes_out[[k]],
                                   function(x) rstanarm::posterior_predict(x, xmat, draws = 1000)))
    mu <- colMeans(preds)
    ci_l <- apply(preds, 2, quantile, ci_lower)
    ci_u <- apply(preds, 2, quantile, ci_upper)
    res[[k]] <- setNames(data.frame(range_est_var, mu, ci_l, ci_u), nm)
  }
  # ----------------------------------------------------------------------------------------------
  return(setNames(res, topic_nam))
}
