\name{agree.wscv}
\title{Within-subject Coefficient of Variation}
\alias{agree.wscv}
\alias{wscv.delta}
\alias{wscv.vst}
\description{
  Obtain confidence interval and point estimate of the
  within-subject coefficient of variation (WSCV).
}
\usage{
  agree.wscv(ratings, conf.level=0.95, method=c("vst", "delta"),
             NAaction=c("fail", "omit"))
}
\arguments{
  \item{ratings}{a matrix of observations with one subject per row and one
	rater per column.}
  \item{conf.level}{confidence level of the interval. The default
	is 0.95.}
  \item{method}{a character string specifying the method used to
	obtain confidence interval of the WSCV.
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
  The point estimate is based on what proposed in Quan and Shih (1996).
  To obtain confidence interval, the methods available include the delta
  method proposed in Quan and Shih (1996) and the variance stabilizing
  transformation in Shoukri et al. (2006).
}
\value{
    Point estimate of the WSCV and lower and upper bounds of the confidence
    interval. 
}
\references{
  Hui Quan and Weichung J. Shih (1996)
  Assessing reproducibility by the within-subject coefficient of
  variation with random effects models.
  \emph{Biometrics}
  \bold{52} 1195-1203
  
  Mohamed M Shoukri, Nasser Elkum and Stephen D Walter (2006)
  Interval estimation and optimal design for the within-subject
  coefficient of variation for continuous and binary variables.
  \emph{BMC Medical Research Methodology}
  \bold{6} 24
} 
\examples{
  data(lesionBurden)
  agree.wscv(lesionBurden.M)
}
\keyword{htest}



