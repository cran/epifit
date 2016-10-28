##' Count NAs in variables.
##'
##' Count NAs, and calculate NA proportion in data.frame.
##' @param data a data.frame to summarize.
##' @return a matrix with total data, NA count and NA proportion.
##' @seealso
##' \code{\link{convertNA}}
##' @examples
##' df <- data.frame(id=1:1000, cov1=rnorm(1000), cov2=runif(1000))
##' df$cov1 <- ifelse(df$cov1 < 0, NA, df$cov1)
##' df$cov2 <- ifelse(df$cov2 < 0.2, NA, df$cov2)
##' countNA(df)
##' @export
countNA <- function(data=NULL){

  if (is.null(data) || !is.data.frame(data)) 
    stop("data is not specified or not data.frame")
  
  n <- nrow(data)
  result <- matrix(n, ncol(data), 4)
  rownames(result) <- colnames(data)
  colnames(result) <- c("observed", "missing", "total", "percent(%)")
  
  for(i in 1:ncol(data)){
    result[i,2] <- sum(as.integer(is.na(data[,i])))
  }
  result[,4] <- result[,2]/result[,3]*100
  result[,1] <- result[,3] - result[,2]
  return(result)
}
