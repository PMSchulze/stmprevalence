
<!-- README.md is generated from README.Rmd. Please edit that file -->

## stmprevalence: Inspect prevalence covariates

This package provides additional functionalities to examine the link
between prevalence covariates and topic proportions. Namely, these are:

  - Implementation of the method of composition, using either a
    quasibinomial GLM or a beta regression. The method of composition is
    implemented in the [stm](http://www.structuraltopicmodel.com/)
    package through its `estimateEffect` function, employing a simple
    OLS regression. This violates the assumption of (sampled) topic
    proportions being restricted to `(0,1)`, which is the motivation for
    our extension.
  - Direct assessment of the prevalence output produced by the stm,
    i.e. the MAP estimates for `Gamma` and `Sigma`, by sampling from a
    LogisticNormal distribution.

Both approaches can be used to visualize the empirical mean as well as
credible intervals of topic proportions for the full observed range of a
specified variable, while holding all other variables as mean/majority
vote.

### Installation

You can install `stmprevalence` from github with:

``` r
# install.packages("devtools")
devtools::install_github("PMSchulze/stmprevalence")
```

### Examples
