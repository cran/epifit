% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pullOneValue.R
\name{pullOneValue}
\alias{pullOneValue}
\title{Pull one set of values from variables included in data which are thought to include essentially the same information.}
\usage{
pullOneValue(data = NULL, varlist = c(""), check = TRUE,
  choose = c("first", "last", "lower", "upper"))
}
\arguments{
\item{data}{a data.frame in which variables are included}

\item{varlist}{a character vector of variable list which is assumed to contain the same information.}

\item{check}{a logical value specifying whether to scan all variables in varlist for incompatible values or not.}

\item{choose}{a character value specifying which to choose among incompatible values. To choose \sQuote{first} or \sQuote{last} depends on the order of the data. \sQuote{lower} and \sQuote{upper} cannot be used for character variables. This option works only when check option is enabled (set to TRUE).}
}
\value{
a vector including one set of values obtained from varlist.
}
\description{
When some part of the data is missing, the missing information may be recovered from another source of information. This function scans all variables which are thought to include essentially the same information, and pulls value from a variable which is not missing.
}
\examples{
dat <- data.frame(a1=c(NA,2,3), a2=c(1,NA,2), a3=c(1,2,NA), b=c(10,11,20))
dat
pullOneValue(dat, c("a1", "a2", "a3"))
}

