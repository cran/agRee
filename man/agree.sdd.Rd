\name{agree.sdd}
\title{Smallest Detectable Difference}
\alias{agree.sdd}
\description{
  Obtain confidence interval and point estimate of the
  smallest detectable difference (SDD).
}
\usage{
  agree.sdd(ratings, conf.level=0.95, NAaction=c("fail", "omit"))
}
\arguments{
  \item{ratings}{a matrix of observations with one subject per row and one
	rater per column.}
  \item{conf.level}{confidence level of the interval. The default
	is 0.95.}
  \item{NAaction}{a character string specifying what should happen
	when the data contain \code{NA}s. It must be one of "fail"
	and "omit" and may be abbreviated. The default is "fail" that causes
	the function to print an error message and terminate if there are
	any incomplete observations. If it is "omit", then the entire row(s)
	containing incomplete observation(s) will be deleted.}
}
\details{
  The calculation is based on one-way random-effects ANOVA and the details can be found in Baumgartner et al. (2015).
}
\value{
    Point estimate of the SDD and lower and upper bounds of the confidence
    interval. 
}
\references{
  Richard Baumgartner, Dai Feng and Aniket Joshi (2015) 
  Determination of smallest detectable difference for PET tracers
  using test-retest data: application in receptor occupancy studies 
  (\emph{under review})
} 
\examples{
  data(petVT)
  agree.sdd(petVT$cerebellum)
}
\keyword{htest}



