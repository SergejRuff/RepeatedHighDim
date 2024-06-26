% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RepeatedHighDim.R
\name{RHighDim}
\alias{RHighDim}
\title{Detection of global group effect}
\usage{
RHighDim(X1, X2, paired = TRUE)
}
\arguments{
\item{X1}{Matrix of expression levels in first group. Rows represent features, columns represent samples.}

\item{X2}{Matrix of expression levels in second group. Rows represent features, columns represent samples.}

\item{paired}{FALSE if samples are unpaired, TRUE if samples are paired.}
}
\value{
An object that contains the test results. Contents can be displayed by the summary function.
}
\description{
Detection of global group effect
}
\details{
Global test for a set of molecular features (e.g. genes, proteins,...) between two experimental groups. Paired or unpaired design is allowed.
}
\examples{
### Global comparison of a set of 100 genes between two experimental groups.
X1 = matrix(rnorm(1000, 0, 1), 10, 100)
X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
RHD = RHighDim(X1, X2, paired=FALSE)
summary_RHD(RHD)
}
\references{
Brunner, E (2009) Repeated measures under non-sphericity. \emph{Proceedings of the 6th St. Petersburg Workshop on Simulation}, 605-609.

Jung K, Becker B, Brunner B and Beissbarth T (2011) Comparison of Global Tests for Functional Gene Sets in Two-Group Designs and Selection of Potentially Effect-causing Genes. \emph{Bioinformatics}, \strong{27}, 1377-1383. \doi{10.1093/bioinformatics/btr152}
}
\seealso{
For more information, please refer to the package's documentation and the tutorial: \url{https://software.klausjung-lab.de/}.
}
\author{
Klaus Jung
}
