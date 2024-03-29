% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.r
\docType{data}
\name{mushroom}
\alias{mushroom}
\title{mushroom data example}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}, \code{symbolic_tbl}) with 23 rows and 3 columns.
}
\usage{
data(mushroom)
}
\description{
mushroom interval data example.
}
\examples{
data(mushroom)
ggInterval_scatter(mushroom, aes(x = Cap.Widths, y = Stipe.Lengths))
}
\references{
Billard L. and  Diday E. (2006).Symbolic data analysis:
Conceptual statistics and data mining. Wiley, Chichester.
}
\keyword{datasets}
