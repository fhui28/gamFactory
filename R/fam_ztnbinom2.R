#' The zero-truncated negative binomial family
#' 
#' @name fam_ztnbinom2
#' @rdname fam_ztnbinom2
#' @export fam_ztnbinom2
#' @examples 
#' \dontrun{
#' # TBD
#' }
#' 
fam_ztnbinom2 <- function(){
  
  bundle <- bundle_ztnbinom2()
  
  fam <- build_family(bundle)()
  
  return(fam)
  
}

