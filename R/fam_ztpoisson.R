#'
#' The zero-truncated Poisson family
#' 
#' @name fam_ztpoisson
#' @rdname fam_ztpoisson
#' @export fam_ztpoisson
#' @examples 
#' 
fam_ztpoisson <- function(){
  
  bundle <- bundle_ztpoisson()
  
  fam <- build_family(bundle)()
  
  return(fam)
  
}

