##' Select variables to keep or drop and a data frame in which variables are kept or dropped is returned.
##'
##' 
##' @title Select variables to keep or drop.
##' @param data a data frame in which variables are included.
##' @param keep a character vector specifying variables to keep.
##' @param drop a character vector specifying variables to drop.
##' @return a data frame in which variables are kept or dropped.
##' @examples
##' df <- data.frame(id=seq(1,10), str=letters[1:10], fac=factor(seq(1,10)), stringsAsFactors=FALSE)
##' selectVariable(df, keep=c("id", "str"))
##' selectVariable(df, drop=c("fac"))
##' @export
selectVariable <- function(data=NULL, keep=NULL, drop=NULL){
  idx <- NULL
  varlist <- names(data)
  if(!is.data.frame(data))
    stop("argument \"data\" must be data frame")

  if(length(drop)==0){
    if(length(keep)==0)
      stop("either drop or keep must be specified")
    # case of keep
    for(var in keep)
      if(!var %in% varlist)
        stop(var , " is not included in data")
    idx <- match(keep, varlist)
  } else {
    if(length(keep) > 0)
      stop("both drop and keep cannot be specified")
    #case of drop
    for(var in drop)
      if(!var %in% varlist)
        stop(var , " is not included in data")
    idx <- -match(drop, varlist)
  }
  return(data[,idx])
}
