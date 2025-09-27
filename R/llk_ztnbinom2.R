#' Log-likelihood of a zero-truncated negative binomial model
#' 
#' @description Many of the derivatives are constructed courtesy of the [Deriv()] package, with close to zero checks made on this after the second derivative!
#' @param param XXX.
#' @param deriv XXX.
#' @rdname llk_ztnbinom2
#' @export llk_ztnbinom2
#' @examples 
#' \dontrun{
#' # Please see the help file for [fam_ztnbinom2()]. 
#' }
#' 
llk_ztnbinom2 <- function(y, 
                          param, 
                          deriv = 0, 
                          ...) {
    if (is.list(param) ) 
        param <- do.call("cbind", param)
    if (is.vector(param)) 
        param <- matrix(param, nrow = 1)
    if (ncol(param) != 2) 
        stop("Wrong number of parameters provided")
    
    p <- ncol( param )
    mu <- param[ , 1, drop = TRUE]
    theta <- param[ , 2, drop = TRUE] # Parametrized such that the variance of a (non-zero-truncated) negative binomial distribution is mu + mu^2 / theta
    n <- length(y)
    
    if (length(mu) == 1) {
        mu <- rep(mu, n)
        theta <- rep(theta, n)
        }
    
    d0 <- lgamma(y + theta) - lgamma(theta) - lfactorial(y) + theta * log(theta) - theta * log(mu + theta) + y * log(mu) - y * log(mu + theta) - log(1 - (theta / (mu + theta))^theta)
    
    out <- list()
    out$d0 <- d0
    
    if( deriv > 0 ) {
        d1 <- y / mu - y / (mu + theta) - (theta^(theta + 1) / (mu + theta)^(theta + 1)) / (1 - (theta / (mu + theta))^theta)
        d2 <- digamma(y + theta) - digamma(theta) + log(theta) + 1 - log(mu + theta) - theta / (mu + theta) - y / (mu + theta) 
        d2_part2 <- function(mu, theta) {
            .e1 <- mu + theta
            .e2 <- theta / .e1
            .e3 <- .e2^theta
            
            ((log(theta) - log(.e1)) * .e3 + theta * (1 - .e2) * .e2^(therta - 1) / .e1) / (1 - .e3)
            }
        d2 <- d2 + d2_part2(mu, theta)
        rm(d2_part2)
        out[["d1"]] <- list(d1, d2)
        
        if( deriv > 1 ){ 
            d11 <- -y / mu^2 + y / (mu + theta)^2
            d11_part2 <- function(mu, theta) {
                .e1 <- theta + mu
                .e2 <- 1 + theta
                .e3 <- theta / .e1
                .e5 <- 1 - .e3^theta
                .e6 <- 3 + theta 
                
                (theta^.e2 * .e2 * .e1^(theta - 2 * .e2) + theta^.e6 * .e3^(theta - 1)/(.e5 * .e1^.e6))/.e5           
                }
            d11 <- d11 + d11_part2(mu, theta)
            rm(d11_part2)
            
            d12 <- y / (mu + theta)^2
            d12_part2 <- function(mu, theta) {
                .e1 <- mu + theta
                .e2 <- 1 + theta
                .e3 <- theta/.e1
                .e4 <- .e3^theta
                .e5 <- .e1^.e2
                .e6 <- 1 - .e4
                .e7 <- log(.e1)
                .e8 <- log(theta)
                
                -(theta^.e2 * (((.e8 - .e7) * .e4 + theta * (1 - .e3) * .e3^(theta - 1)/.e1)/.e6 + .e2/theta + .e8 - (.e2 * .e1^theta + .e5 * .e7)/.e1^(2 * .e2 - .e2))/(.e6 * .e5))      
                }
            d12 <- d12 + d12_part2(mu, theta)
            rm(d12_part2)
            
            d22 <- trigamma(y + theta) - trigamma(theta) + 1/theta - 1/(mu + theta) - (theta * (y - mu)) / (mu + theta)^2
            d22_part2 <- function(mu, theta) {
                .e1 <- mu + theta
                .e2 <- theta/.e1
                .e3 <- .e2^theta
                .e4 <- 1 - .e2
                .e5 <- log(.e1)
                .e6 <- log(theta)
                .e7 <- theta - 1
                .e8 <- .e2^.e7
                .e10 <- (.e6 - .e5) * .e3 + theta * .e4 * .e8/.e1
                .e11 <- 1 - .e3
                .e12 <- 1/.e1
                
                ((.e10/.e11 + .e6 - .e5) * .e10 + (1 + theta * ((.e4 * .e7/theta + .e6 - (.e12 + .e5)) * .e4 - (2 - .e2)/.e1)) * .e8/.e1 + (1/theta - .e12) * .e3) / .e11
                }
            d22 <- d22 + d22_part2(mu, theta)
            rm(d22_part2)
            out[["d2"]] <- list(d11, d12, d22) 
            
            if( deriv > 2 ){
                d111 <- 2 * y / mu^3 - 2 * y / (mu + theta)^3
                d111_part2 <- function(mu, theta) {
                    .e1 <- theta + mu
                    .e2 <- theta/.e1
                    .e3 <- 1 + theta
                    .e5 <- 1 - .e2^theta
                    .e6 <- theta - 1
                    .e7 <- 3 + theta
                    .e8 <- .e2^.e6
                    .e9 <- 2 * .e3
                    .e10 <- .e5 * .e1^.e7
                    .e11 <- 2 + theta
                    .e12 <- theta - .e9
                    .e13 <- theta^.e3
                    .e14 <- theta^2
                    .e14 * (theta^.e6 * .e3 * .e12 * .e1^(theta - (1 + .e9)) - (((.e13 * .e3 * .e1^.e12 + theta^.e7 * .e8/.e10)/(.e5 * .e1^2) + .e13 * 
                                                                                     (.e5 * .e7 * .e1^.e11 + .e14 * .e1^.e3 * .e8)/.e10^2) * 
                                                                                    .e8 + theta^.e11 * .e6 * .e2^(theta - 2)/(.e5 * .e1^(5 + theta)))) / .e5
                    }
                d111 <- d111 + d111_part2(mu, theta)
                rm(d111_part2)
                
                d112 <- -2 * y / (mu + theta)^3
                d112_part2 <- function(mu, theta) {
                    .e1 <- mu + theta
                    .e2 <- theta/.e1
                    .e3 <- 1 + theta
                    .e4 <- 3 + theta
                    .e5 <- .e2^theta
                    .e6 <- theta - 1
                    .e7 <- .e2^.e6
                    .e8 <- .e1^.e4
                    .e9 <- 1 - .e5
                    .e10 <- 2 * .e3
                    .e11 <- log(.e1)
                    .e12 <- log(theta)
                    .e13 <- theta - .e10
                    .e14 <- theta^.e3
                    .e15 <- theta^.e4
                    .e16 <- .e1^.e13
                    .e17 <- 1 - .e2
                    .e18 <- .e12 - .e11
                    .e19 <- .e9 * .e8
                    .e21 <- .e18 * .e5 + theta * .e17 * .e7/.e1
                    .e22 <- 2 + theta
                    
                    ((.e21 * (.e14 * .e3 * .e16 + .e15 * .e7/.e19) + (.e7 * (theta^.e22 * .e4 + .e15 * .e12) + .e15 * (.e17 * .e6 * .e2^(theta - 2)/.e1 + .e18 * .e7))/.e8)/.e9 + 
                            (.e3 * (.e14 * .e12 + theta^theta * .e3) + .e14) * .e16 + .e14 * (.e1^(theta - (1 + .e10)) * .e13 - .e11 * .e16) * .e3 - .e15 * 
                            ((.e4 * .e1^.e22 + .e11 * .e8) * .e9 - .e21 * .e8) * .e7/.e19^2) / .e9
                    }
                d112 <- d112 + d112_part2(mu, theta)
                rm(d112_part2)
                
                d122 <- -2 * y / (mu + theta)^3
                d122_part2 <- function(mu, theta) {
                    .e1 <- mu + theta
                    .e2 <- theta/.e1
                    .e3 <- 1 + theta
                    .e4 <- log(.e1)
                    .e5 <- .e2^theta
                    .e6 <- .e1^.e3
                    .e7 <- log(theta)
                    .e8 <- theta - 1
                    .e9 <- 1 - .e2
                    .e10 <- .e1^theta
                    .e11 <- .e2^.e8
                    .e12 <- 1 - .e5
                    .e13 <- .e7 - .e4
                    .e15 <- .e3 * .e10 + .e4 * .e6
                    .e17 <- .e13 * .e5 + theta * .e9 * .e11/.e1
                    .e18 <- .e17/.e12
                    .e19 <- .e3/theta
                    .e20 <- theta^.e3
                    .e24 <- .e18 + .e19 + .e7 - .e15/.e6
                    .e25 <- .e12 * .e6
                    
                    -((.e24 * (.e20 * .e7 + theta^theta * .e3) + .e20 * (((.e18 + .e7 - .e4) * .e17 + 
                                                                              ((1 - theta * (1 + 2 * .e9)/.e1) * .e11 + theta * (.e9 * .e8 * .e2^(theta - 2)/.e1 + .e13 * .e11) * .e9)/.e1 + 
                                                                              (1/theta - 1/.e1) * .e5)/.e12 + .e15^2/.e1^(2 * .e3) + (2 - .e19)/theta - 
                                                                             (.e15 * .e4 + .e3 * (.e4 * .e10 + theta * .e1^.e8) + 2 * .e10)/.e6))/.e25 - .e20 * (.e15 * .e12 - .e17 * .e6) * .e24/.e25^2)
                }
                d122 <- d122 + d122_part2(mu, theta)
                rm(d122_part2)
                
                d222 <- psigamma(y + theta, deriv = 2) - psigamma(theta, deriv = 2) - 1 / theta^2 + 1 / (mu + theta)^2 - (y * mu - mu^2 - theta * y + theta * mu) / (mu + theta)^3
                d222_part2 <- function(mu, theta) {
                    .e1 <- mu + theta
                    .e2 <- theta/.e1
                    .e3 <- 1 - .e2
                    .e4 <- theta - 1
                    .e5 <- log(.e1)
                    .e6 <- log(theta)
                    .e7 <- .e2^theta
                    .e8 <- .e2^.e4
                    .e9 <- .e3 * .e4
                    .e10 <- .e9/theta
                    .e11 <- 1/.e1
                    .e12 <- .e6 - .e5
                    .e14 <- .e12 * .e7 + theta * .e3 * .e8/.e1
                    .e16 <- .e10 + .e6 - (.e11 + .e5)
                    .e17 <- 1 - .e7
                    .e18 <- .e16 * .e3
                    .e19 <- (2 - .e2)/.e1
                    .e23 <- .e14/.e17 + .e6 - .e5
                    .e25 <- (1 + theta * (.e18 - .e19)) * .e8/.e1
                    .e26 <- (1/theta - .e11) * .e7
                    .e27 <- .e1^2
                    .e28 <- 2/.e1
                    
                    ((.e16 * (2 + theta * (.e18 - (3 - .e2)/.e1)) + theta * (((2 - ((.e9 + theta)/.e1 + .e10))/theta - (.e10 + 1 + .e6 - (.e28 + .e5))/.e1) * 
                                                                                 .e3 + (3 - 2 * .e2)/.e27) - .e19) * .e8/.e1 + (.e14 * .e12 + .e25 + .e26) * .e23 + 
                            .e14 * (2 * ((.e23 * .e14 + .e25 + .e26)/.e17) + 2/theta - .e28) + (1/.e27 - 1/theta^2) * .e7)/.e17
                    }
                d222 <- d222 + d222_part2(mu, theta)
                rm(d222_part2)
                out[["d3"]] <- list(d111, d112, d122, d222) 
                
                if( deriv > 3){
                    d1111 <- -6 * y / mu^4 + 6 * y / (mu + theta)^4 
                    d1111_part2 <- function(mu, theta) {
                        .e1 <- mu + theta
                        .e2 <- theta/.e1
                        .e3 <- 1 + theta
                        .e5 <- 1 - .e2^theta
                        .e6 <- theta - 1
                        .e7 <- 3 + theta
                        .e8 <- .e2^.e6
                        .e9 <- 2 * .e3
                        .e10 <- theta^2
                        .e11 <- .e1^.e7
                        .e12 <- theta^.e3
                        .e13 <- .e5 * .e11
                        .e14 <- theta - 2
                        .e15 <- .e1^.e3
                        .e16 <- 2 + theta
                        .e17 <- theta - .e9
                        .e18 <- .e1^2
                        .e19 <- .e2^.e14
                        .e22 <- .e5 * .e7 * .e1^.e16 + .e10 * .e15 * .e8
                        .e23 <- .e5 * .e18
                        .e24 <- 5 + theta
                        .e25 <- .e12 * .e3
                        .e26 <- theta^.e7
                        .e27 <- .e13^2
                        .e28 <- .e5 * .e1^.e24
                        .e31 <- theta - (1 + .e9)
                        .e33 <- .e25 * .e1^.e17 + .e26 * .e8/.e13
                        .e34 <- .e1^.e31
                        .e35 <- .e1^theta
                        .e37 <- .e33/.e23 + .e12 * .e22/.e27
                        
                        .e10 * (theta * (.e12 * ((.e5 * .e24 * .e1^(4 + theta) + .e10 * .e11 * .e8) * .e19/.e28^2 + theta * .e14 * .e2^(theta - 3)/(.e5 * .e1^(7 + theta))) * .e6 + 
                                             theta^.e14 * .e3 * .e1^(theta - (2 + .e9)) * .e31 * .e17 - 
                                             (theta * .e8 * (theta^.e6 * .e3 * .e34 * .e17 - (.e37 * .e8 + theta^.e16 * .e6 * .e19/.e28))/.e5 - .e37 * .e6 * .e19)/.e18) - 
                                    ((.e25 * .e34 * .e17 - .e26 * (.e22 * .e8/(.e5 * .e15) + theta * .e6 * .e19)/.e28)/.e23 + 
                                         .e12 * ((.e5 * .e16 * .e15 + .e10 * .e35 * .e8) * .e7 + .e10 * (.e3 * .e35 * .e8 - theta * .e1^.e6 * .e6 * .e19) - 
                                                     2 * (.e22^2/.e13))/.e27 - (2 * .e5 + .e10 * .e8/.e1) * .e1 * .e33/.e23^2) * .e8) / .e5
                    }    
                    d1111 <- d1111 + d1111_part2(mu, theta)
                    rm(d1111_part2)
                    
                    d1112 <- 6 * y / (mu + theta)^4 
                    d1112_part2 <- function(mu, theta) {
                        .e1 <- mu + theta
                        .e2 <- theta/.e1
                        .e3 <- 1 + theta
                        .e4 <- .e2^theta
                        .e5 <- theta - 1
                        .e6 <- 1 - .e4
                        .e7 <- 3 + theta
                        .e8 <- .e2^.e5
                        .e9 <- log(.e1)
                        .e10 <- 2 * .e3
                        .e11 <- 2 + theta
                        .e12 <- log(theta)
                        .e13 <- .e1^.e7
                        .e14 <- 1 - .e2
                        .e15 <- theta^.e3
                        .e16 <- .e6 * .e13
                        .e17 <- .e12 - .e9
                        .e18 <- theta - .e10
                        .e19 <- theta - 2
                        .e20 <- .e1^.e3
                        .e21 <- .e1^.e11
                        .e23 <- .e17 * .e4 + theta * .e14 * .e8/.e1
                        .e24 <- .e2^.e19
                        .e25 <- 5 + theta
                        .e26 <- .e1^.e25
                        .e27 <- .e1^.e18
                        .e29 <- theta - (1 + .e10)
                        .e30 <- theta^.e11
                        .e31 <- theta^.e5
                        .e32 <- .e6 * .e1^2
                        .e33 <- .e1^.e29
                        .e34 <- theta^.e7
                        .e35 <- .e16^2
                        .e38 <- .e6 * .e7 * .e21 + theta^2 * .e20 * .e8
                        .e39 <- .e6 * .e26
                        .e40 <- .e14 * .e5
                        .e43 <- .e15 * .e3 * .e27 + .e34 * .e8/.e16
                        .e44 <- ((.e7 * .e21 + .e9 * .e13) * .e6 - .e23 * .e13)/.e16
                        .e47 <- .e40 * .e24/.e1 + .e17 * .e8
                        .e49 <- .e43/.e32 + .e15 * .e38/.e35
                        .e50 <- .e31 * .e3
                        
                        theta * ((2 + theta * .e23/.e6) * (.e50 * .e33 * .e18 - (.e49 * .e8 + .e30 * .e5 * .e24/.e39)) + 
                                     theta * (((.e3 * (.e31 * .e12 + theta^.e19 * .e5) + .e31) * .e18 - .e50) * .e33 + 
                                                  .e30 * ((.e25 * .e1^(4 + theta) + .e9 * .e26) * .e6 - .e23 * .e26) * .e5 * .e24/.e39^2 + 
                                                  .e31 * (.e1^(theta - (2 + .e10)) * .e29 - .e9 * .e33) * .e3 * .e18 - 
                                                  ((((.e3 * (.e15 * .e12 + theta^theta * .e3) + .e15) * .e27 + .e15 * (.e33 * .e18 - .e9 * .e27) * .e3 + 
                                                         .e34 * ((.e40 + 3 + theta)/theta + 2 * .e12 - (.e44 + .e9)) * .e8/.e16)/.e32 + 
                                                        .e15 * (.e38 * (.e3/theta + .e12 - 2 * .e44) + (.e11 * .e20 + .e9 * .e21) * .e6 * .e7 + 
                                                                    (1 - (.e23 * .e7 + .e4)) * .e21 + theta * ((2 * .e20 + theta * (.e3 * .e1^theta + .e9 * .e20)) * .e8 + theta * .e47 * .e20))/.e35 - 
                                                        (2 * .e6 - .e23 * .e1) * .e1 * .e43/.e32^2) * .e8 + ((.e5 * (.e15 * .e11 + .e30 * .e12) + .e30) * .e24 + 
                                                                                                                 .e30 * (.e14 * .e19 * .e2^(theta - 3)/.e1 + .e17 * .e24) * .e5)/.e39 + 
                                                       .e47 * .e49))) / .e6
                    }
                    d1112 <- d1112 + d1112_part2(mu, theta)
                    rm(d1112_part2)
                    
                    d1222 <- 6 * y / (mu + theta)^4
                    d1222_part2 <- function(mu, theta) {
                        .e1 <- mu + theta
                        .e2 <- theta/.e1
                        .e3 <- 1 + theta
                        .e4 <- log(.e1)
                        .e5 <- theta - 1
                        .e6 <- log(theta)
                        .e7 <- 1 - .e2
                        .e8 <- .e2^theta
                        .e9 <- .e2^.e5
                        .e10 <- .e1^theta
                        .e11 <- .e1^.e3
                        .e12 <- .e6 - .e4
                        .e14 <- .e3 * .e10 + .e4 * .e11
                        .e16 <- .e12 * .e8 + theta * .e7 * .e9/.e1
                        .e17 <- 1 - .e8
                        .e18 <- theta - 2
                        .e19 <- .e1^.e5
                        .e20 <- .e2^.e18
                        .e21 <- .e7 * .e5
                        .e24 <- .e21 * .e20/.e1 + .e12 * .e9
                        .e25 <- theta * .e19
                        .e26 <- .e16/.e17
                        .e28 <- .e4 * .e10 + .e25
                        .e29 <- .e3/theta
                        .e30 <- 2 * .e7
                        .e31 <- theta^.e3
                        .e32 <- .e3 * .e28
                        .e35 <- 1/theta - 1/.e1
                        .e36 <- theta * .e24
                        .e37 <- (1 - 2 * .e2) * .e9
                        .e38 <- .e35 * .e8
                        .e39 <- 2 * .e3
                        .e40 <- theta^theta
                        .e42 <- .e14 * .e4 + .e32
                        .e46 <- .e26 + .e29 + .e6 - .e14/.e11
                        .e48 <- .e26 + .e6 - .e4
                        .e49 <- 1 + .e30
                        .e51 <- .e31 * .e6 + .e40 * .e3
                        .e54 <- (.e37 + 2 * .e36) * .e7 + (1 - theta * .e49/.e1) * .e9
                        .e56 <- .e14 * .e17 - .e16 * .e11
                        .e57 <- .e42 + 2 * .e10
                        .e58 <- 2 * (.e48 * .e16)
                        .e59 <- 2 * .e38
                        .e60 <- (.e54/.e1 + .e58 + .e59)/.e17
                        .e61 <- ((.e14 * (1/.e11 + 1/.e1^(.e39 - .e3)) - 2 * .e4) * .e14 - (2 * .e32 + 4 * .e10))/.e11
                        .e65 <- (.e37 + .e36) * .e7/.e1 + .e16 * .e12 + .e38
                        .e66 <- .e17 * .e11
                        .e68 <- .e1^.e39
                        .e69 <- 2 * ((2 - .e29)/theta)
                        .e70 <- 2 + .e30
                        .e71 <- theta^2
                        
                        -(((.e60 + .e61 + .e69) * .e51 + .e46 * (.e6 * .e51 + .e40 * ((1 + .e6) * .e3 + 2)) + 
                               .e31 * (((.e60 + 2/theta - 2/.e1) * .e16 + .e65 * .e48 + 
                                            ((.e24 + theta * (((.e7 * .e18 * .e2^(theta - 3)/.e1 + .e12 * .e20) * .e7 * .e5 + (1 - (2 * .e21 + theta)/.e1) * .e20)/.e1 + 
                                                                  .e24 * .e12 + .e35 * .e9)) * .e7 + 
                                                 .e24 * (1 - theta * (.e70 - .e2)/.e1) - ((.e70 - theta * (.e30 + 2 * .e49)/.e1) * .e9 + .e36 * .e7)/.e1)/.e1 + 
                                            (1/.e1^2 - 1/.e71) * .e8)/.e17 + .e14 * (2 * (.e57/.e68) - .e14 * (2 * (.e3 * .e1^(.e39 - 1)) + 2 * (.e4 * .e68))/.e1^(4 * .e3)) - 
                                           (((.e42 + 3 * .e10) * .e4 + .e14 * (1 - .e57/.e10)/.e1 + .e3 * (2 * .e19 + .e4 * .e28 + theta * (.e1^.e18 * .e5 + .e4 * .e19)) + 
                                                 2 * .e28 + .e25)/.e11 + (3 - 2 * .e29)/.e71)) - ((.e56 * .e51 + .e31 * (.e57 * .e17 - (.e65 * .e11 + 2 * (.e14 * .e16)))) * .e46 + 
                                                                                                      (.e46 * .e51 + .e31 * (((.e54 - 2 * (.e56 * .e46/.e10))/.e1 + .e58 + .e59)/.e17 + .e61 + .e69)) * .e56) / .e66) / .e66)
                    }
                    d1222 <- d1222 + d1222_part2(mu, theta)
                    rm(d1222_part2)
                    
                    d1122 <- 6 * y / (mu + theta)^4
                    d1122_part2 <- function(mu, theta) {
                        .e1 <- mu + theta
                        .e2 <- theta/.e1
                        .e3 <- 1 + theta
                        .e4 <- 3 + theta
                        .e5 <- log(.e1)
                        .e6 <- log(theta)
                        .e7 <- theta - 1
                        .e8 <- .e2^theta
                        .e9 <- .e2^.e7
                        .e10 <- .e1^.e4
                        .e11 <- 1 - .e2
                        .e12 <- 2 * .e3
                        .e13 <- 2 + theta
                        .e14 <- .e6 - .e5
                        .e15 <- 1 - .e8
                        .e16 <- theta - .e12
                        .e17 <- theta^.e3
                        .e18 <- theta^.e4
                        .e20 <- .e14 * .e8 + theta * .e11 * .e9/.e1
                        .e21 <- .e1^.e13
                        .e22 <- .e11 * .e7
                        .e23 <- .e1^.e16
                        .e24 <- .e15 * .e10
                        .e26 <- .e4 * .e21 + .e5 * .e10
                        .e28 <- theta - (1 + .e12)
                        .e29 <- .e1^.e28
                        .e30 <- theta - 2
                        .e31 <- theta^.e13
                        .e32 <- .e2^.e30
                        .e33 <- 1/.e1
                        .e34 <- theta^theta
                        .e36 <- .e26 * .e15 - .e20 * .e10
                        .e39 <- .e22 * .e32/.e1 + .e14 * .e9
                        .e41 <- .e29 * .e16 - .e5 * .e23
                        .e43 <- .e17 * .e6 + .e34 * .e3
                        .e45 <- .e31 * .e4 + .e18 * .e6
                        .e47 <- 1/theta - .e33
                        .e48 <- 2 * .e6
                        .e49 <- .e36/.e24
                        .e50 <- .e24^2
                        .e52 <- (.e22 + 3 + theta)/theta + .e48
                        .e57 <- .e20 * .e14 + (1 + theta * ((.e22/theta + .e6 - (.e33 + .e5)) * .e11 - (2 - .e2)/.e1)) * .e9/.e1 + .e47 * .e8
                        .e58 <- .e3 * .e43
                        .e60 <- .e9 * .e45 + .e18 * .e39
                        .e63 <- .e17 * .e3 * .e23 + .e18 * .e9/.e24
                        
                        ((((.e4 * (.e17 * .e13 + .e31 * .e6) + 2 * .e31 + .e6 * .e45) * .e9 + 2 * (.e39 * .e45) + 
                               .e18 * (((.e11 * .e30/theta + .e6 - (2/.e1 + .e5)) * .e7 + 1) * .e11 * .e32/.e1 + .e39 * .e14 + .e47 * .e9) - .e26 * .e60/.e10)/.e10 + 
                              .e57 * .e63 + .e20 * (2 * ((.e20 * .e63 + .e60/.e10)/.e15) + 2 * ((.e58 + .e17) * .e23) + 
                                                        .e17 * (2 * (.e41 * .e3) + theta^2 * ((.e52 - (.e49 + .e5))/.e24 - .e36/.e50) * .e9)))/.e15 +
                                (.e41 * .e43 + .e17 * ((.e1^(theta - (2 + .e12)) * .e28 - .e5 * .e29) * .e16 - (.e41 * .e5 + 2 * .e29))) * .e3 + 
                                (.e3 * (.e6 * .e43 + .e34 * ((1 + .e6) * .e3 + 3)) + .e17 * (.e3/theta + .e48)) * .e23 + 
                                (.e58 + 2 * .e17) * .e41 - .e18 * ((.e52 - (2 * .e49 + .e5)) * .e36 + 
                                                                       ((.e13 * .e1^.e3 + .e5 * .e21) * .e4 + .e26 * .e5 + 2 * .e21) * .e15 - (.e57 * .e10 + 2 * (.e26 * .e20))) * .e9/.e50)/.e15    
                    }
                    d1122 <- d1122 + d1122_part2(mu, theta)
                    rm(d1122_part2)
                    
                    d222 <- psigamma(y + theta, deriv = 2) - psigamma(theta, deriv = 2) - 1 / theta^2 + 1 / (mu + theta)^2 - (y * mu - mu^2 - theta * y + theta * mu) / (mu + theta)^3
                    d2222 <- psigamma(y + theta, deriv = 3) - psigamma(theta, deriv = 3) + 2 / theta^3 - 2 / (mu + theta)^3 - (2 * y * theta - 2 * mu * theta - 4 * y * mu + 4 * mu^2) / (mu + theta)^4
                    d2222_part2 <- function(mu, theta) {
                        .e1 <- mu + theta
                        .e2 <- theta/.e1
                        .e3 <- 1 - .e2
                        .e4 <- theta - 1
                        .e5 <- log(.e1)
                        .e6 <- log(theta)
                        .e7 <- .e3 * .e4
                        .e8 <- .e2^theta
                        .e9 <- .e7/theta
                        .e10 <- .e2^.e4
                        .e11 <- 1/.e1
                        .e12 <- .e6 - .e5
                        .e14 <- .e9 + .e6 - (.e11 + .e5)
                        .e16 <- .e12 * .e8 + theta * .e3 * .e10/.e1
                        .e17 <- .e14 * .e3
                        .e18 <- 1 - .e8
                        .e19 <- (2 - .e2)/.e1
                        .e20 <- (.e7 + theta)/.e1
                        .e22 <- .e1^2
                        .e23 <- 1 + theta * (.e17 - .e19)
                        .e25 <- 1/theta - .e11
                        .e28 <- .e23 * .e10/.e1
                        .e29 <- .e25 * .e8
                        .e30 <- (2 - (.e20 + .e9))/theta
                        .e31 <- 2/.e1
                        .e34 <- .e16/.e18 + .e6 - .e5
                        .e35 <- 2 * .e2
                        .e37 <- (.e30 - (.e9 + 1 + .e6 - (.e31 + .e5))/.e1) * .e3
                        .e38 <- 3 - .e35
                        .e42 <- .e16 * .e12 + .e28 + .e29
                        .e43 <- .e38/.e22
                        .e44 <- theta * (.e17 - (3 - .e2)/.e1)
                        .e45 <- theta^2
                        .e46 <- (.e34 * .e16 + .e28 + .e29)/.e18
                        .e48 <- 1/.e22
                        .e49 <- 1/.e45
                        .e50 <- theta * (.e37 + .e43)
                        .e51 <- ((.e17 + .e50 - (3 + .e44)/.e1) * .e10 + (.e7 * .e2^(theta - 2)/.e1 + .e12 * .e10) * .e23)/.e1
                        .e52 <- .e42 * .e34
                        .e53 <- .e16 * (2 * .e46 + 2/theta - .e31)
                        .e54 <- (.e49 - .e48) * .e8
                        .e55 <- .e1^3
                        .e56 <- 2 * .e20
                        .e57 <- 2 * .e9
                        .e58 <- 2 + .e44
                        
                        ((.e51 + .e42 * .e12 + 2 * (.e16 * .e25) - .e54) * .e34 + (((.e14 * .e58 + .e50 - .e19) * .e10/.e1 + .e52 + 
                                                                                        .e53 + (.e48 - .e49) * .e8 + 2 * (.e51 + .e52 + .e53 - .e54))/.e18 + 3/.e22 - 3/.e45) * .e16 + 
                                ((.e14 * (3 + theta * (.e17 - (4 - .e2)/.e1)) + theta * ((7 - 4 * .e2)/.e22 + 2 * .e37) - (5 - .e35)/.e1) * .e14 + .e37 + 
                                     (.e30 - (1 - .e11)/.e1) * .e58 + 2 * .e43 - theta * ((((2 - .e56)/.e1 + (3 - (.e56 + .e57))/theta)/theta + 
                                                                                               (2 * .e30 - (.e57 + 2 * .e6 + 3 - (2 * .e5 + 6/.e1))/.e1)/.e1) * .e3 + 
                                                                                              (2 * .e3 + 2 * .e38)/.e55)) * .e10/.e1 + .e42 * (3 * .e46 + 3/theta - 3/.e1) + 
                                (2/theta^3 - 2/.e55) * .e8) / .e18    
                    }
                    d2222 <- d2222 + d2222_part2(mu, theta)
                    rm(d2222_part2)
                    out[["d4"]] <- list(d1111, d1112, d1122, d1222, d2222)
                    }
                }
            }
        }
    
    return( out )
    }

