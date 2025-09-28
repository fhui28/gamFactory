#'
#' Predict using nested effects
#' 
#' @name Predict.matrix.nested
#' @rdname Predict.matrix.nested
#' @param object XXX
#' @param data XXX
#' @param ... Can be ignored.
#' @export
#'
Predict.matrix.nested <- function(object, data, ...){
  
  if(class(object)[1] == "si"){
    return( .predict.matrix.si(object, data, ...) ) 
  }
  if(class(object)[1] == "nexpsm"){
    return( .predict.matrix.nexpsm(object, data, ...) ) 
  }
  if(class(object)[1] == "mgks"){
    return( .predict.matrix.mgks(object, data, ...) ) 
  }
  
  stop("I do not know this effect type")
  
}