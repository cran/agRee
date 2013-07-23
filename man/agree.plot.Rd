\name{agree.plot}
\title{Visualize the Agreement of Ratings among Different Raters}
\alias{agree.plot}
\description{
  Draw Bland-Altman plot(s) and scatter plot(s) with identity line.
}
\usage{
  agree.plot(ratings, NAaction=c("fail", "omit"))
}
\arguments{
  \item{ratings}{a matrix of ratings from different raters, one
	rater per column.}
  \item{NAaction}{a character string specifying what should happen
	when the data contain \code{NA}s. It must be one of "fail"
	and "omit" and may be abbreviated. The default is "fail" that causes
	the function to print an error message and terminate if there are
	any incomplete observations. If it is "omit", then the entire row(s)
	containing incomplete observation(s) will be deleted.}
}
\details{
  The function produces a matrix of plots. The upper panel consists
  of scatter plot(s) with identity line. The lower panel consists
  of the Bland-Altman plot(s) with confidence bounds and bias using dotted
  line in red color and the horizontal line passing through the origin
  in black.
}
\note{
  The confidence bounds are mean of the difference between two raters
  plus or minus twice of the SD of difference.
}
\value{
  \code{NULL}
}
\references{
  J. Martin Bland and Douglas G. Altman (1986)
  Statistical methods for assessing agreement between
  two methods of clinical measurement.
  \emph{Lancet} \bold{1} 307-310
} 
\examples{
  data(judgeRatings)
  agree.plot(judgeRatings)
}
\keyword{hplot}



