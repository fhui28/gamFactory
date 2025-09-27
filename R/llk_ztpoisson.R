#' Log-likelihood of a zero-truncated Poisson model
#' 
#' @description XXX.
#' @param param XXX.
#' @param deriv XXX.
#' @rdname llk_ztpoisson
#' @export llk_ztpoisson
#' @examples 
#' \dontrun{
#' # Please see the help file for [fam_ztpoisson()]. 
#' }
#' 
llk_ztpoisson <- function(y, 
                          param, 
                          deriv = 0, 
                          ...) {
    
    if (is.list(param) ) 
        param <- do.call("cbind", param)
    if (is.vector(param)) 
        param <- matrix(param, ncol = 1)
    if (ncol(param) != 1) 
        stop("Wrong number of parameters provided")
    
    p <- ncol( param )
    mu <- param[ , 1, drop = TRUE] # rate parameter is Poisson distribution
    n <- length(y)
    
    if (length(mu) == 1) {
        mu <- rep(mu, n)
    }
    
    d0 <- y * log(mu) - log(exp(mu) - 1) - lfactorial(y)
    
    out <- list()
    out$d0 <- d0
    
    if( deriv > 0 ) {
        d1 <- y / mu - exp(mu) / (exp(mu) - 1)
        out[["d1"]] <- list(d1)
        
        if( deriv > 1 ){
            d11 <- - y / mu^2 + exp(mu) / (exp(mu) - 1)^2
            out[["d2"]] <- list(d11) 
            
            if( deriv > 2 ){
                d111 <- 2 * y / mu^3 - exp(mu) * (exp(mu) + 1) / (exp(mu) - 1)^3
                out[["d3"]] <- list(d111) 
                
                if( deriv > 3){
                    d1111 <- - 6 * y / mu^4 + exp(mu) * (exp(2*mu) + 4 * exp(mu) + 1) / (exp(mu) - 1)^4
                    out[["d4"]] <- list(d1111)
                }
            }
        }
    }
    
    return( out )
    }
