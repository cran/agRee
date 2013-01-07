\name{agree.ccc}
\title{Concordance Correlation Coefficient}
\alias{agree.ccc}
\alias{ccc.mvn.mcmc}
\alias{ccc.nonpara.jackknife}
\description{
  Obtain confidence interval and point estimate of the
  concordance correlation coefficient (CCC) proposed in Lin (1989).
}
\usage{
  agree.ccc(ratings, conf.level=0.95,
            method=c("jackknifeZ", "jackknife",
                     "mvn.jeffreys", "mvn.conjugate"),
            nsim=10000)
}
\arguments{
  \item{ratings}{a matrix of observations with one subject per row and one
	rater per column.}
  \item{conf.level}{confidence level of the interval. The default
	is 0.95.}
  \item{method}{a character string specifying the method used to
	obtain the estimate of the CCC.
    It must be one of "jackknifeZ", "jackknife", "mvn.jeffreys" and
	"mvn.conjugate" 
    and may be abbreviated.	The default is "jackknifeZ".}
  \item{nsim}{number of iterations used in the Bayesian approach with
	Jeffreys or conjugate prior.}
}
\details{
  To obtain point estimate and confidence interval, the methods available include the
  jackknife method with and without Z-transformation and the Bayesian
  approach with Jeffreys and conjugate priors for the multivariate
  normal distribution.
}
\value{
    Point estimate and lower and upper bounds of the confidence
    interval of the CCC. 
}
\references{

  Dai Feng, Vladimir Svetnik, Alexandre Coimbra and Richard Baumgartner (2012)
  A comparison of confidence interval methods for the concordance
  correlation coefficient and intraclass correlation coefficient with
  small number of raters.
  \emph{Journal of Biopharmaceutical Statistics (accepted)}

  Lawrence I-Kuei Lin (1989)
  A concordance correlation coefficient to evaluate reproducibility.
  \emph{Biometrics} \bold{45} 255-268

} 
\seealso{
  \code{\link[epiR]{epi.ccc}}, \code{\link[cccrm]{cccvc}},
  \code{\link[Agreement]{agreement}}.
}
\examples{
  data(judgeRatings)
  agree.ccc(judgeRatings[,2:3])
}
\keyword{htest}



