#' @title The zero-truncated negative binomial family
#' @description The `fam_ztnbinom2` family implements a truncated negative binomial model in which the linear predictor controls the mean parameter of the negative binomial distribution, and *not* the mean of the distribution itself.
#' 
#' Requires positive integer count data.
#' 
#' @details
#' The zero-truncated negative binomial distribution is a two-parameter distribution for modeling count data that does not include zeros, whre the (non-truncated) negative binomial distribution is parametrized such that \eqn{Var(Y) = \mu + \mu^2 / \theta}, where \eqn{\mu} is the mean and \eqn{\theta} is the dispersion parameter. 
#' 
#' The fitted values for this family represent the negative binomial mean parameter, \eqn{\mu}, of the distribution, and *not* the mean of the distribution itself, which is given by \eqn{E[Y] = \mu / (1 - P(y = 0))}, where \eqn{P(y = 0) = (\theta / (\mu + \theta))^\theta}.
#' 
#' @return 
#' A family inheriting from glass `general.family`
#' 
#' @name fam_ztnbinom2
#' @rdname fam_ztnbinom2
#' @export fam_ztnbinom2
#' @examples 
#' \dontrun{
#' library(gamFactory)
#' library(sdmTMB)
#' library(actuar)
#' 
#' # Example 1: Generate data from and fitting a zero-truncated NB GLM
#' dat <- gamSim(1, n = 2000, dist = "poisson", scale = .1)
#' dat$eta = 0 + 0.5*dat$x0 + 0.5*dat$x1 - 0.5*dat$x2 + 0.5*dat$x3
#' dat$mu <- exp(dat$eta)
#' dat$y <- rztnbinom(n = 2000, size = 1/2, prob = 1/2/(dat$mu + 1/2))
#' 
#' table(dat$y)
#' 
#' 
#' ## Fit misspecified NB model
#' fit_wrong <- gam(y ~ x0 + x1 + x2 + x3,
#' data = dat, 
#' method = "REML",
#' family = nb())
#' 
#' ## Reference fit using sdmTMB's truncated_nbinom family for fitting truncated NB models
#' fit_gold <- sdmTMB(y ~ x0 + x1 + x2 + x3,
#' data = dat,
#' spatial = "off", 
#' reml = TRUE,
#' family = truncated_nbinom2())
#' 
#' ## Create specific zero-truncated NB family and fit it
#' fit <- gam(list(y ~ x0 + x1 + x2 + x3, ~1), 
#' data = dat, 
#' method = "REML",
#' family = fam_ztnbinom2())
#' 
#' fit_wrong$coefficients # Estimated parameters should be off
#' coef(fit_gold); exp(fit_gold$sd_report$par.fixed)
#' fit$coefficients # Estimated parameters match those from sdmTMB and be close to true parameters
#' 
#' 
#' # Example 2: Generate data from and fit a zero-truncated NB GAM
#' dat <- gamSim(1, n = 2000, dist = "poisson", scale = .1)
#' dat$eta = 0.1*(dat$f0 + dat$f1 + dat$x2 + dat$x3)
#' dat$mu <- exp(dat$eta)
#' dat$y <- rztnbinom(n = 2000, size = 1/2, prob = 1/2/(dat$mu + 1/2))
#' 
#' 
#' ## Reference fit using sdmTMB's truncated_nbinom family for fitting truncated NB models
#' fit_gold <- sdmTMB(y ~ s(x0) + s(x1) + s(x2) + s(x3),
#' data = dat,
#' spatial = "off", 
#' reml = TRUE,
#' family = truncated_nbinom2())
#' 
#' ## Create specific zero-truncated NB family and fit it
#' fit <- gam(list(y ~ s(x0) + s(x1) + s(x2) + s(x3), ~ 1), 
#' data = dat
#' method = "REML",
#' family = fam_ztnbinom2())
#' 
#' 
#' err <- abs(fit$fitted.values - exp(fit_gold$fitted.values[1:2000,1]))
#' summary(err)
#' # Error here should be small
#' }
#' 
fam_ztnbinom2 <- function(){
  
  bundle <- bundle_ztnbinom2()
  
  fam <- build_family(bundle)()
  
  return(fam)
  
}


