##' Convert factor variable to numeric or character variable.
##'
##' This function convert factor variable to numeric or character variable. In \sQuote{automatic} mode, all variables not specified by numeric or character arguments are converted into numeric or character automatically. Explicitly specified variables are converted as specified. In \sQuote{only_specified} mode, variables specified by numeric argument is converted into numeric variables, and the character is conveted into character variables. When converted into numeric variable, NAs can be produced for incompatible data with warnings. In \sQuote{all_character}, and \sQuote{all_numeric} modes, all factor variables are converted into character or numeric variables, respectively.
##' @title Convert factor variable to numeric or character variable.
##' @param data a data.frame which contains factor variables.
##' @param numeric a character vector specifying variable names to be converted into numeric variables.
##' @param character a character vector specifying variable names to be converted into character variables.
##' @param mode a character vector specifying converting mode. See details and the default is \sQuote{only_specified}.
##' @return a converted data.frame.
##' @seealso \code{\link{showContents}}, \code{\link{listNumericIncompatibility}}
##' @examples
##' a <- factor(rnorm(5))
##' b <- c("a", "b", "c", "d", "e")
##' c <- c("1", "2", "3", "4", NA)
##' d <- c("1", "2", "3", "4", ".")
##' dat <- data.frame(a,b,c,d)
##' dat2 <- convertFromFactor(dat)
##' dat3 <- convertFromFactor(dat, numeric=c("d"))
##' dat4 <- convertFromFactor(dat, mode="all_character")
##' @export
convertFromFactor <- function(data=NULL, numeric=c(""), character=c(""), mode=c("automatic", "only_specified", "all_character", "all_numeric")){
  
  if(!is.data.frame(data))
    stop("data argument must be data.frame")

  mode <- match.arg(mode)

  varlist <- colnames(data)

  numflag <- numeric %in% varlist
  if(sum(numflag) != length(numeric) && numeric != ""){
    warning("non-existent variable (", paste0(numeric[!numflag], collapse=", "), ") is specified in numeric argument")
  }
  numeric <- numeric[numflag]
  
  charflag <- character %in% varlist
  if(sum(charflag) != length(character) && character != ""){
    warning("non-existent variable (",  paste0(character[!charflag], collapse=", "),") is specified in character argument")
  }
  character <- character[charflag]

  if(mode == "automatic"){

    automatic <- varlist[!varlist %in% c(numeric, character)]

  } else if(mode == "all_character"){

    if(length(numeric) > 0)
      warning("Variables specified in numeric argument are ignored")
    automatic <- NULL
    numeric <- NULL
    character <- varlist
    
  } else {

    if(length(numeric) > 0)
      warning("Variables specified in character argument are ignored")
    automatic <- NULL
    numeric <- varlist
    character <- NULL
    
  }

  for(var in automatic){
    content <- data[[var]]
    if(is.factor(content)){
      tryCatch(
        {
          data[[var]] <- as.numeric(as.character(content))
        },
        warning=function(e){
          data[[var]] <<- as.character(content)
        }
        )
    }
  }
  
  for(var in numeric){
    content <- data[[var]]
    if(is.factor(content)){
      data[[var]] <- as.numeric(as.character(content))
    }
  }

  for(var in character){
    content <- data[[var]]
    if(is.factor(content)){
      data[[var]] <- as.character(content)
    }
  }

  return(data)
}
