# sigmoid function; inverse of logit-link
sigmoid <- function(x) exp(x)/(1+exp(x))

# obtain majority vote for categorial variable
majority <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# create design matrix where all columns but est_var fixed as median/majority values
make_median_xmat <- function(est_var, formula, metadata, range_est_var) {
  dat <- metadata[, -which(names(metadata) == est_var)]
  dat_new <- lapply(dat, function(x) if(is.numeric(x)) median(x) else majority(x))
  dat_fit <- data.frame(dat_new, range_est_var)
  names(dat_fit) <- c(names(dat_new),est_var)
  f <- paste("~",as.character(formula)[3])
  return(stm::makeDesignMatrix(as.formula(f), metadata, dat_fit))
}

# helper function to sample from normal distribution of regression coefficients;
## allowed types are beta regression or quasibinomial glm
sample_normal <- function(mod, type) {
  mu <- var <- NULL
  if (type == "beta"){
    mu <- mod$coefficients$mean
    var <- mod$vcov[1:length(mu), 1:length(mu)]
  } else if (type == "quasibinomial") {
    mu <- mod$coefficients
    var <- vcov(mod)
  } else {
    stop("Error: Please set type='beta' or type='quasibinomial'")
  }
  return(mvtnorm::rmvnorm(1, mean = mu, sigma = var))
}

# helper function to simulate theta (either mean or quantile p) from LogisticNormal
## for single document, given mean mu_d of gaussian
sim_theta_d <- function(mu_d, Sigma, nsims, p = "mean") {
  eta_sim <- cbind(mvtnorm::rmvnorm(nsims, mu_d, Sigma),0)
  theta_sim <- exp(eta_sim - matrixStats::rowLogSumExps(eta_sim))
  colnames(theta_sim) <- c()
  return(apply(theta_sim, 2, function(x) if(p=="mean") mean(x) else quantile(x, probs = p)))
}

# helper function to simulate mean as well as quantiles ci_lower, ci_upper of LogisticNormal,
## for all documents
sim_theta <- function(mu, Sigma, nsims, ci_lower, ci_upper) {
  topic_n <- ncol(mu)+1
  mean_emp <- apply(mu, 1, sim_theta_d, Sigma = Sigma, p = "mean", nsims = nsims)
  ci_lower <- apply(mu, 1, sim_theta_d, Sigma = Sigma, p = ci_lower, nsims = nsims)
  ci_upper <- apply(mu, 1, sim_theta_d, Sigma = Sigma, p = ci_upper, nsims = nsims)
  nm <- c("proportion", "ci_lower", "ci_upper")
  res <- lapply(1:topic_n,
                function(i) setNames(data.frame(mean_emp[i,], ci_lower[i,], ci_upper[i,]), nm))
  return(setNames(res, paste0("Topic", 1:topic_n)))
}
