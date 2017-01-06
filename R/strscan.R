##' Extract words from a string vector.
##'
##' This function split a string into several string components divided by delimiter, and devided components selected by position are returned.
##' @param src a source string vector to be split.
##' @param delimiter a delimiter string.
##' @param pos a subset indicator vector to select.
##' @return a string vector to be extracted.
##' @examples
##' bdate <- c("3/21/1941", "9/21/1919", "6/1/1951")
##' datestr <- paste(strscan(bdate, "/", 3), strscan(bdate, "/", 1), strscan(bdate, "/", 2),
##'                  sep="-")
##' datestr
##' as.Date(bdate, '%m/%d/%Y')
##' @export
strscan <- function(src, delimiter, pos){
  seplst <- strsplit(src, delimiter)
  numelem <- sapply(seplst, length)
  maxelem <- max(numelem)
  if(min(numelem) != maxelem)
    for(i in 1:length(seplst))
      length(seplst[[i]]) <- maxelem
  
  ret <- unlist(seplst)
  dim(ret) <- c(maxelem, length(ret)/maxelem)
  return(t(ret)[,pos])
}
