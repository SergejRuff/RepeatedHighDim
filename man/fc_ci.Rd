% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_ci.R
\name{fc_ci}
\alias{fc_ci}
\title{Calculation of adjusted confidence intervals}
\usage{
fc_ci(fit, alpha = 0.05, method = "raw")
}
\arguments{
\item{fit}{Object as returned from the function eBayes of the limma package}

\item{alpha}{1 - confidence level (e.g., if confidence level is 0.95, alpha is 0.05)}

\item{method}{Either 'raw' for unadjusted confidence intervals, or 'BH' for Bejamini Hochberg-adjusted confidence intervals, or 'BY' for Benjamini Yekutieli-adjusted confidence intervals}
}
\value{
A results matrix with one row per gene, and one column for the p-value, the log fold change, the lower limit of the CI, and the upper limit of the CI
}
\description{
Calculation of adjusted confidence intervals
}
\details{
Calculation of unadjusted and adjusted confidence intervals for the log fold change
}
\examples{
### Artificial microarray data
d = 1000 ### Number of genes
n = 10 ### Sample per group
fc = rlnorm(d, 0, 0.1)
mu1 = rlnorm(d, 0, 1) ### Mean vector group 1
mu2 = mu1 * fc ### Mean vector group 2
sd1 = rnorm(d, 1, 0.2)
sd2 = rnorm(d, 1, 0.2)
X1 = matrix(NA, d, n) ### Expression levels group 1
X2 = matrix(NA, d, n) ### Expression levels group 2
for (i in 1:n) {
  X1[,i] = rnorm(d, mu1, sd=sd1)
  X2[,i] = rnorm(d, mu2, sd=sd2)
}
X = cbind(X1, X2)
heatmap(X)

### Differential expression analysis with limma
if(check_limma()){
group = gl(2, n)
design = model.matrix(~ group)
fit1 = limma::lmFit(X, design)
fit = limma::eBayes(fit1)

### Calculation of confidence intervals
CI = fc_ci(fit=fit, alpha=0.05, method="raw")
head(CI)
CI = fc_ci(fit=fit, alpha=0.05, method="BH")
head(CI)
CI = fc_ci(fit=fit, alpha=0.05, method="BY")
head(CI)

fc_plot(CI, xlim=c(-0.5, 3), ylim=-log10(c(1, 0.0001)), updown="up")
fc_plot(CI, xlim=c(-3, 0.5), ylim=-log10(c(1, 0.0001)), updown="down")
fc_plot(CI, xlim=c(-3, 3), ylim=-log10(c(1, 0.0001)), updown="all")
}
}
\references{
Dudoit, S., Shaffer, J. P., & Boldrick, J. C. (2003). Multiple hypothesis testing in microarray experiments. \emph{Statistical Science}, \strong{18(1)}, 71-103. \url{https://projecteuclid.org/journals/statistical-science/volume-18/issue-1/Multiple-Hypothesis-Testing-in-Microarray-Experiments/10.1214/ss/1056397487.full}

Jung, K., Friede, T., & Beißbarth, T. (2011). Reporting FDR analogous confidence intervals for the log fold change of differentially expressed genes. \emph{BMC bioinformatics}, \strong{12}, 1-9. \url{https://link.springer.com/article/10.1186/1471-2105-12-288}
}
\seealso{
For more information, please refer to the package's documentation and the tutorial: \url{https://software.klausjung-lab.de/}.
}
\author{
Klaus Jung
}
