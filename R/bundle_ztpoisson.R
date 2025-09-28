#' Bundle for zero-truncated Poisson regression model
#' 
#' @name bundle_ztpoisson
#' @rdname bundle_ztpoisson
#' @importFrom actuar rztpois
#' @export bundle_ztpoisson
#'
bundle_ztpoisson <- function(){
    out <- list(np = 1,
                available_deriv = 4,
                llk = gamFactory::llk_ztpoisson,
                links = list(c("log", "sqrt")), 
                nam = "ztpoisson",
                bundle_nam = as.character(match.call()[[1]]),
                residuals = function(object, type = c("deviance", "pearson", "response")) {
                    type <- match.arg(type)
                    r <- .resid_ztpoisson(object, type)
                    return( r )
                },
                rd = function(mu, wt, scale) {
                    return( actuar::rztpois(nrow(mu), lambda = mu) )
                },
                initialize = function(y, nobs, E, x, family, offset, jj, unscaled){
                    if(any(y <= 0)) {
                        stop("All response variables must be greater than zero for ztpoisson.")
                    }
                    if(all.equal(y, round(y)) != TRUE) {
                        stop("Non-integer response variables are not allowed with ztpoisson.")
                    }
                    n <- rep(1, nobs)
                    ## should E be used unscaled or not?..
                    use.unscaled <- if(!is.null(attr(E, "use.unscaled"))) { 
                        TRUE 
                    } 
                    else { 
                        FALSE 
                    }
                    
                    lpi <- attr(x, "lpi")
                    start <- rep(0, ncol(x))
                    
                    #### Get scale phi = sigma = E(y) (if xi == 0)
                    yt1 <- if (family$link[[1]]=="identity"){ 
                        y 
                    } else {
                        family$linfo[[1]]$linkfun(as.double(y) + 1e-3)
                    }
                    x1 <- x[ , lpi[[1]], drop=FALSE]
                    e1 <- E[ , lpi[[1]], drop=FALSE] ## square root of total penalty
                    
                    if (use.unscaled) {
                        qrx <- qr( rbind(x1, e1) )
                        x1 <- rbind(x1, e1)
                        startji <- qr.coef(qr(x1), c(yt1,rep(0,nrow(E))))
                        startji[!is.finite(startji)] <- 0       
                    } else {
                        startji <- penreg(x1, e1, yt1)
                    }
                    start[lpi[[1]]] <- startji
                    
                    return( start )
                }
    )
    
    # Fixing the environment of all functions
    for(ii in 1:length(out)){
        if(inherits(out[[ii]], "function")) {
            environment(out[[ii]]) <- environment()
        }
    }
    
    return( out )
}



#' @noMd
#' @noRd
.resid_ztpoisson <- function(object, type) { 
    y <- drop(object$y)
    mu <- object$fitted.values
    Ey <- mu / (1 - exp(-mu))  # mean of zero-truncated Poisson
    wts <- object$prior.weights
    
    if(type == "deviance") {
        rsd <- 2 * (wts * (y * log(y/mu) + log(exp(mu - 1)) - log(exp(y) - 1)))
        s <- sign(y - Ey)
        rsd <- sqrt(pmax(rsd, 0)) * s
        } 
    else {
        rsd <- y - Ey
        if(type == "pearson") {
            rsd <- rsd * sqrt(wts / (Ey * (1 + mu - Ey)))
            }
        }
    return(rsd)
}
