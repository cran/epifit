##' Calculate the difference between two date in terms of unit of time.
##'
##' This function calculate the difference between two date in terms of unit of time, and age can be obtained when \sQuote{year} is specified as unit argument.
##' @param birthday a character or character vector specifying birthday or base date.
##' @param targetdate a character specifying current or target date.
##' @param unit a character specifying unit for calculating the difference between the two dates. Values of "year", "month" and "day" are supported.
##' @return a vector of age
##' @examples calcAge("1963-2-3")
##' @examples calcAge("1970-1-1", unit="day")
##' @export
calcAge <- function(birthday, targetdate=Sys.Date(), unit="year"){
  sapply(birthday,
         function(x){
           tryCatch(
             {length(seq(as.Date(x), as.Date(targetdate), unit)) - 1},
             error=function(e){NA})
         },
         USE.NAMES=FALSE
         )
}
