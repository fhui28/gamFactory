#' Bundle for zero-truncated negative binomial model
#' 
#' @name bundle_ztnbinom2
#' @rdname bundle_ztnbinom2
#' @importFrom actuar rztnbinom
#' @export bundle_ztnbinom2
#'
bundle_ztnbinom2 <- function() {
    out <- list(np = 2,
                available_deriv = 4,
                llk = gamFactory::llk_ztnbinom2,
                links = list(c("log"), c("identity")), 
                nam = "ztnbinom2",
                bundle_nam = as.character(match.call()[[1]]),
                residuals = function(object, type = c("deviance", "pearson", "response")) {
                    type <- match.arg(type)
                    r <- .resid_ztnbinom2(object, type)
                    return( r )
                },
                rd = function(mu, wt, scale) {
                    return( actuar::rztnbinom(nrow(mu), 
                                              size = mu[,2,drop=TRUE], 
                                              prob = mu[,2,drop=TRUE] / (mu[,1,drop=TRUE] + mu[,2,drop=TRUE])) )
                },
                initialize = function(y, nobs, E, x, family, offset, jj, unscaled){
                    if(any(y <= 0)) {
                        stop("All response variables must be greater than zero for ztnbinom2.")
                        }
                    if(all.equal(y, round(y)) != TRUE) {
                        stop("Non-integer response variables are not allowed with ztnbinom2.")
                        }
                    
                    n <- rep(1, nobs)
                    ## should E be used unscaled or not?..
                    use.unscaled <- if(!is.null(attr(E, "use.unscaled"))) { 
                        TRUE 
                    } 
                    else { 
                        FALSE 
                    }
                    
                    jj <- attr(x, "lpi")
                    start <- rep(0, ncol(x))
                    x1 <- x[, jj[[2]], drop = FALSE]
                    e1 <- E[, jj[[2]], drop = FALSE]
                    yt1 <- as.numeric(as.logical(y))
                    if (use.unscaled) {
                        qrx <- qr(rbind(x1, e1))
                        x1 <- rbind(x1, e1)
                        startji <- qr.coef(qr(x1), c(yt1, rep(0, nrow(E))))
                        startji[!is.finite(startji)] <- 0
                    } else startji <- penreg(x1, e1, yt1)
                    start[jj[[2]]] <- startji
                    p <- drop(x1[1:nobs, , drop = FALSE] %*% startji)
                    ind <- y == 0 & p < 0.5
                    w <- rep(1, nobs)
                    w[ind] <- 0.1
                    yt1 <- family$linfo[[1]]$linkfun(log(abs(y) + (y == 0) * 0.2))
                    yt1 <- yt1 * w
                    x1 <- w * x[, jj[[1]], drop = FALSE]
                    e1 <- E[, jj[[1]], drop = FALSE]
                    if (use.unscaled) {
                        x1 <- rbind(x1, e1)
                        startji <- qr.coef(qr(x1), c(yt1, rep(0, nrow(E))))
                        startji[!is.finite(startji)] <- 0
                        } 
                    else startji <- penreg(x1, e1, yt1)
                    start[jj[[1]]] <- startji
                    
                    return(start)
                    }
                )
    
    # Fixing the environment of all functions
    for(ii in 1:length(out)){
        if( class(out[[ii]]) == "function" ){
            environment(out[[ii]]) <- environment()
        }
    }
    
    
    return( out )
}


#' @noMd
#' @noRd
.resid_ztnbinom2 <- function(object, type) { 
    y <- drop(object$y)
    mu <- object$fitted.values
    theta <- object$theta #' Makes an assumption that this is in the object
    nb_prob <- theta / (mu + theta) 
    Ey <- theta * (1-nb_prob) / (nb_prob * (1-nb_prob^theta)) # mean of zero-truncated negative binomial distribution
    wts <- object$prior.weights
    
    if(type == "deviance") {
        rsd <- 2 * (wts * (y * log(y/mu) + y * log((mu + theta) / (y + theta)) + log(1 - (theta / (theta + mu))^theta) - log(1 - (theta / (theta + y))^theta) ))
        s <- sign(y - Ey)
        rsd <- sqrt(pmax(rsd, 0)) * s
        } 
    else {
        rsd <- y - Ey
        if(type == "pearson") {
            getvar <- (theta * (1-nb_prob) * (1 - (1 + theta * (1-nb_prob)) * nb_prob^theta)) / (nb_prob * (1-nb_prob^theta))^2
            rsd <- rsd * sqrt(wts / getvar)
            rm(getvar)
            }
        }
    return(rsd)
    }



