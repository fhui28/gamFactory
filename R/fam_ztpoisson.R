#' @title The zero-truncated Poisson family
#' @description The `fam_ztpoisson` family implements a truncated Poisson model in which the linear predictor controls the rate (mean) parameter of the Poisson distribution, and *not* the mean of the distribution itself. 
#' 
#' Requires positive integer count data.
#' 
#' @details
#' The zero-truncated Poisson distribution is a one-parameter distribution for modeling count data that does not include zeros. The fitted values for this family represent the Poisson rate (mean) parameter, \eqn{\mu}, of the distribution, and *not* the mean of the distribution itself, which is given by \eqn{E[Y] = \mu / (1 - \exp(-\mu))}.
#' 
#' @return 
#' A family inheriting from glass `general.family`
#' 
#' @name fam_ztpoisson
#' @rdname fam_ztpoisson
#' @export fam_ztpoisson
#' @examples 
#' 
#' \dontrun{
#' library(gamFactory)
#' library(actuar)
#' 
#' # Example 1: Generate data from and fitting a zero-truncated Poisson GLM
#' dat <- gamSim(1, n = 2000, dist = "poisson", scale = .1)
#' dat$eta = 0 + 0.5*dat$x0 + 0.5*dat$x1 - 0.5*dat$x2 + 0.5*dat$x3
#' dat$mu <- exp(dat$eta)
#' dat$y <- rztpois(n = 2000, lambda = dat$mu)
#' 
#' table(dat$y)
#' 
#' 
#' ## Fit misspecified Poisson model
#' fit_wrong <- gam(y ~ x0 + x1 + x2 + x3,
#' data = dat, 
#' method = "REML",
#' family = poisson())
#' 
#' ## Reference fit using mgcv's ziplss family for fitting hurdle Poisson models
#' fit_gold <- gam(list(y ~ x0 + x1 + x2 + x3, ~ 1),
#' data = rbind(dat, 0), # Adding a zero as ziplss can sometimes hang in datas without zeros 
#' method = "REML",
#' family = ziplss())
#' 
#' ## Create specific zero-truncated poisson family and fit it
#' fit <- gam(list(y ~ x0 + x1 + x2 + x3), 
#' data = dat, 
#' method = "REML",
#' family = fam_ztpoisson())
#' 
#' fit_wrong$coefficients # Estimated parameters should be off
#' fit_gold$coefficients[1:5] 
#' fit$coefficients # Estimated parameters match ziplss GAM and be close to true parameters
#' 
#' 
#' # Example 2: Generate data from a Poisson GAM and fit a zero-truncated Poisson GAM
#' dat <- gamSim(1, n = 2000, dist = "poisson", scale = .1)
#' dat$eta = 0.1*(dat$f0 + dat$f1 + dat$f2 + dat$f3)
#' dat$mu <- exp(dat$eta)
#' dat$y <- rztpois(n = 2000, lambda = dat$mu)
#' 
#' 
#' ## Reference fit using mgcv's ziplss family for fitting hurdle Poisson models
#' fit_gold <- gam(list(y ~ s(x0) + s(x1) + s(x2) + s(x3), ~ 1),
#' data = rbind(dat, 0), # Adding a zero as ziplss can sometimes hang in datas without zeros 
#' family = ziplss())
#' 
#' ## Create specific zero-truncated poisson family and fit it
#' fit <- gam(list(y ~ s(x0) + s(x1) + s(x2) + s(x3)), 
#' data = dat,
#' method = "REML",
#' family = fam_ztpoisson())
#' 
#' 
#' err <- abs(fit$fitted.values - exp(fit_gold$fitted.values[1:2000,1]))
#' summary(err)
#' # Error here should be relatively small
#' }
#' 

fam_ztpoisson <- function(){
  
  bundle <- bundle_ztpoisson()
  
  fam <- build_family(bundle)()
  
  return(fam)
  
}

