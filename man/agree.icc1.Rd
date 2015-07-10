\name{agree.icc1}
\title{Intraclass correlation coefficient for one-way random anova model}
\alias{agree.icc1}
\alias{icc1.sf}
\description{
  Obtain confidence interval and point estimate of the
  intraclass correlation coefficient for one-way random anova model (ICC1).
}
\usage{
  agree.icc1(ratings, conf.level=0.95, method=c("sf"),
             NAaction=c("fail", "omit"))
}
\arguments{
  \item{ratings}{a matrix of observations with one subject per row and one
	rater per column.}
  \item{conf.level}{confidence level of the interval. The default
	is 0.95.}
  \item{method}{a character string specifying the method used to
	obtain confidence interval of the ICC1.
    Now only the "sf" method is supported.}
  \item{NAaction}{a character string specifying what should happen
	when the data contain \code{NA}s. It must be one of "fail"
	and "omit" and may be abbreviated. The default is "fail" that causes
	the function to print an error message and terminate if there are
	any incomplete observations. If it is "omit", then the entire row(s)
	containing incomplete observation(s) will be deleted.}
}
\details{
  The point estimate and confidence interval are
  based on a one-way random anova model as proposed in Shrout and Fleiss (1979).
}
\value{
    Point estimate of the ICC1 and lower and upper bounds of the confidence
    interval. 
}
\references{
  Patrick E Shrout and Joseph L Fleiss (1979).
  Intraclass correlations: uses in assessing rater reliability.
  \emph{Psychological Bulletin}
  \bold{86} 420-428
} 
\examples{
  data(lesionBurden)
  agree.icc1(lesionBurden.M)
}
\keyword{htest}



