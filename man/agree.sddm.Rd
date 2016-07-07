\name{agree.sddm}
\title{Mean Normalized Smallest Detectable Difference}
\alias{agree.sddm}
\description{
  Obtain confidence interval and point estimate of the
  mean normalized smallest detectable difference (SDDm).
}
\usage{
  agree.sddm(ratings, conf.level=0.95, method=c("vst", "delta"),
             NAaction=c("fail", "omit"))
}
\arguments{
  \item{ratings}{a matrix of observations with one subject per row and one
	rater per column.}
  \item{conf.level}{confidence level of the interval. The default
	is 0.95.}
  \item{method}{a character string specifying the method used to
	obtain confidence interval of the WSCV, based on what the SDDm is calculated.
    It must be one of "vst" and "delta" and may be abbreviated.
	The default is "vst".}
  \item{NAaction}{a character string specifying what should happen
	when the data contain \code{NA}s. It must be one of "fail"
	and "omit" and may be abbreviated. The default is "fail" that causes
	the function to print an error message and terminate if there are
	any incomplete observations. If it is "omit", then the entire row(s)
	containing incomplete observation(s) will be deleted.}
}
\details{
   The calculation is based on the relationship with the WSCV and the details can be found in Baumgartner et al. (2015).
}
\value{
    Point estimate of the SDDm and lower and upper bounds of the confidence
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
  agree.sddm(petVT$cerebellum)
}
\keyword{htest}



