% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countNA.R
\name{countNA}
\alias{countNA}
\title{Count NAs in variables.}
\usage{
countNA(data = NULL)
}
\arguments{
\item{data}{a data.frame to summarize.}
}
\value{
a matrix with total data, NA count and NA proportion.
}
\description{
Count NAs, and calculate NA proportion in data.frame.
}
\examples{
df <- data.frame(id=1:1000, cov1=rnorm(1000), cov2=runif(1000))
df$cov1 <- ifelse(df$cov1 < 0, NA, df$cov1)
df$cov2 <- ifelse(df$cov2 < 0.2, NA, df$cov2)
countNA(df)
}
\seealso{
\code{\link{convertNA}}
}

