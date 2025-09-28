#'
#' Predictions from GAM models with non-linear effects
#' 
#' @name predict.gamnl
#' @rdname predict.gamnl
#' @param object XXX
#' @param newdata XXX
#' @param type XXX
#' @param se.fit XXX
#' @param ... Can be ignored
#' @importFrom mgcv Predict.matrix
#' @export
#'
predict.gamnl <- function(object, newdata, type="link", se.fit=FALSE, ...){
  
  out <- predict.gam(object = object, newdata = newdata, type = type, block.size = Inf, ...)
  
  return(out)
  
}