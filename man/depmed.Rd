% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depmed.R
\name{depmed}
\alias{depmed}
\title{Calculates the depth median.}
\usage{
depmed(G)
}
\arguments{
\item{G}{List containing the grid information produced by
\code{\link{gridfun}} and the halfspace location depths
produced by \code{\link{hldepth}}.}
}
\value{
An vector with a length equal to the number of dimension
    of the array in G, containing the coordinates of the depth
    median.
}
\description{
Calculates the depth median.
}
\details{
Calculates the depth median in a specified grid array with given
halfspace location depth at each grid location.
}
\examples{
## Attention: calculation is currently time-consuming.
\dontrun{

# A 3-dimensional example data set D1
n <- 200
x1 <- rnorm(n, 0, 1)
y1 <- rnorm(n, 0, 1)
z1 <- rnorm(n, 0, 1)
D1 <- data.frame(cbind(x1, y1, z1))
colnames(D1) <- c("x", "y", "z")

# Specification of the grid and calculation of the halfspace location depth at each grid location.
G <- gridfun(D1, grid.size=20)
G$H <- hldepth(D1, G, verbose=TRUE)
dm <- depmed(G) ## Calculation of the depth median
}
}
\references{
Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
    bagplot: a bivariate boxplot. The American Statistician,
    53(4), 382-387.
}
\seealso{
For more information, please refer to the package's documentation and the tutorial: \url{https://software.klausjung-lab.de/}.
}
\author{
Jochen Kruppa, Klaus Jung
}
