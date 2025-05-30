% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_names.r
\name{match_names}
\alias{match_names}
\title{Match names}
\usage{
match_names(x, where)
}
\description{
Find an elements of \code{x} in the first space-separated element of \code{where}.
}
\details{
The typical use of this function is to match sequence names.
In fasta files, sequence names often consist of the sequence name itself and an additional information about
the sequence.
}
\examples{
\dontrun{
match_names("foo", "foo bar")
match_names(c("foo", "bar"), c("bar", "foo"))
match_names(c("foo", "bar") # throws error: foo not found
match_names("foo", c("foo", "foo")) # throws error: foo found twice
}
}
\keyword{internal}
