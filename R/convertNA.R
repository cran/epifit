##' Convert a character pattern or numeric value into NA and vice versa.
##'
##' Convert a character pattern or numeric value into NA and vice versa.
##' @param data a data.frame to summarize.
##' @param na.character a character vector specifying missing character.
##' @param na.numeric a numeric vector specifying missing value.
##' @param fromNA a logical value specifying replacement that NA is replaced with the first element of na.character or na.numeric.
##' @return a data.frame with NA replacement.
##' @seealso
##' \code{\link{countNA}}
##' @examples
##' dat <- data.frame(a=c("","2","3"),b=c("4", NA, "."), c=c(-1,1,3),
##'                   d=c(NA,3,2), stringsAsFactors=FALSE)
##' dat2 <- convertNA(dat, na.character=c("", "."))
##' dat3 <- convertNA(dat, na.character=".", na.numeric=-1, fromNA=TRUE)
##' dat
##' dat2
##' dat3
##' @export
convertNA <- function(data=NULL, na.character=c(""), na.numeric=NA, fromNA=FALSE){

  if (is.null(data) || !is.data.frame(data)) 
    stop("data is not specified or not data.frame")
  
  for(col in 1:ncol(data)){
    
    if(is.numeric(data[,col])){ # numeric

      if(fromNA){

        for(row in 1:nrow(data)){
          if(is.na(data[row, col])){
            data[row, col] <- na.numeric[1]
          }
        }
        
      } else {

        for(row in 1:nrow(data)){
          if(data[row, col] %in% na.numeric){
            data[row, col] <- NA
          }
        }
        
      }
    
    } else { # factor or character
      
      if(fromNA){

        for(row in 1:nrow(data)){
          if(is.na(data[row, col])){
            data[row, col] <- na.character[1]
          }
        }
        
      } else {
        
        for(row in 1:nrow(data)){
          if(as.character(data[row, col]) %in% na.character){
            data[row, col] <- NA
          }
        }
        
      }
    }
  }
  return(data)
}
