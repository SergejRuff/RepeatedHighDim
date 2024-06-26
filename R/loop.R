#' Calculates the fence and the loop of a gemplot (i.e. the outer gemstone).
#'
#' The fence inflates the the bag relative to the depth median by the
#' factor inflation. Data points outside the bag and inside the fence
#' the loop or outer gemstone are flagged as outliers. Data points
#' outside the fence are marked as outliers. In the case of a
#' 3-dimensional data set, the loop can be visualized by an outer
#' gemstone around the inner gemstone or bag.
#' @title Calculates the fence and the loop
#' @param D Data set with rows representing the individuals and
#'     columns representing the features. In the case of three
#'     dimensions, the colnames of D must be c("x", "y", "z").
#' @param B List containing the information about the coordinates of
#'     the bag and the convex hull that forms the bag (determined by
#'     \code{\link{bag}}).
#' @param inflation A numeric value > 0 that specifies the inflation
#'     factor of the bag relative to the median (default = 3).
#' @param dm The coordinates of the depth median as produced by
#'     \code{\link{depmed}}.
#' @return A list containing the following elements:
#' \describe{
#' \item{\emph{coords.loop}}{Coordinates of the data points that are inside the convex hull around the loop.}
#' \item{\emph{hull.loop}}{A data matrix that contains the indices of the margin data points of the loop that cover the convex hull by triangles. Each row represnts one triangle. The indices correspond to the rows of coords.loop.}
#' \item{\emph{coords.fence}}{Coordinates of the grid points that are inside the fence but outside the bag.}
#' \item{\emph{hull.fence}}{A data matrix that contains the indices of the margin grid points of the fence that cover the convex hull around the fence by triangles. Each row represnts one triangle. The indices correspond to the rows of coords.fence.}
#' \item{\emph{outliers}}{A vector of length equal to the sample size. Data points that are inside the fence are labelled by 0 and values outside the fence (i.e. outliers) are labelled by 1.}
#' }
#'
#' @references
#' Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The bagplot: a bivariate boxplot. \emph{The American Statistician}, \strong{53(4)}, 382-387. \doi{10.1080/00031305.1999.10474494}

#'
#' Kruppa, J., & Jung, K. (2017). Automated multigroup outlier identification in molecular high-throughput data using bagplots and gemplots. \emph{BMC bioinformatics}, \strong{18(1)}, 1-10. \url{https://link.springer.com/article/10.1186/s12859-017-1645-5}
#' @author Jochen Kruppa, Klaus Jung
#' @importFrom  rgl material3d bg3d points3d text3d spheres3d axes3d
#' @export
#' @examples
#' ## Attention: calculation is currently time-consuming.
#' \dontrun{
#'
#' # Two 3-dimensional example data sets D1 and D2
#' n <- 200
#' x1 <- rnorm(n, 0, 1)
#' y1 <- rnorm(n, 0, 1)
#' z1 <- rnorm(n, 0, 1)
#' D1 <- data.frame(cbind(x1, y1, z1))
#' x2 <- rnorm(n, 1, 1)
#' y2 <- rnorm(n, 1, 1)
#' z2 <- rnorm(n, 1, 1)
#' D2 <- data.frame(cbind(x2, y2, z2))
#' colnames(D1) <- c("x", "y", "z")
#' colnames(D2) <- c("x", "y", "z")
#'
#' # Placing outliers in D1 and D2
#' D1[17,] = c(4, 5, 6)
#' D2[99,] = -c(3, 4, 5)
#'
#' # Grid size and graphic parameters
#' grid.size <- 20
#' red <- rgb(200, 100, 100, alpha = 100, maxColorValue = 255)
#' blue <- rgb(100, 100, 200, alpha = 100, maxColorValue = 255)
#' yel <- rgb(255, 255, 102, alpha = 100, maxColorValue = 255)
#' white <- rgb(255, 255, 255, alpha = 100, maxColorValue = 255)
#' require(rgl)
#' material3d(color=c(red, blue, yel, white),
#'  alpha=c(0.5, 0.5, 0.5, 0.5), smooth=FALSE, specular="black")
#'
#' # Calucation and visualization of gemplot for D1
#' G <- gridfun(D1, grid.size=20)
#' G$H <- hldepth(D1, G, verbose=TRUE)
#' dm <- depmed(G)
#' B <- bag(D1, G)
#' L <- loop(D1, B, dm=dm)
#' bg3d(color = "gray39" )
#' points3d(D1[L$outliers==0,1], D1[L$outliers==0,2], D1[L$outliers==0,3], col="green")
#' text3d(D1[L$outliers==1,1], D1[L$outliers==1,2], D1[L$outliers==1,3],
#' as.character(which(L$outliers==1)), col=yel)
#' spheres3d(dm[1], dm[2], dm[3], col=yel, radius=0.1)
#' material3d(1,alpha=0.4)
#' gem(B$coords, B$hull, red)
#' gem(L$coords.loop, L$hull.loop, red)
#' axes3d(col="white")
#'
#' # Calucation and visualization of gemplot for D2
#' G <- gridfun(D2, grid.size=20)
#' G$H <- hldepth(D2, G, verbose=TRUE)
#' dm <- depmed(G)
#' B <- bag(D2, G)
#' L <- loop(D2, B, dm=dm)
#' points3d(D2[L$outliers==0,1], D2[L$outliers==0,2], D2[L$outliers==0,3], col="green")
#' text3d(D2[L$outliers==1,1], D2[L$outliers==1,2], D2[L$outliers==1,3],
#' as.character(which(L$outliers==1)), col=yel)
#' spheres3d(dm[1], dm[2], dm[3], col=yel, radius=0.1)
#' gem(B$coords, B$hull, blue)
#' gem(L$coords.loop, L$hull.loop, blue)
#' }
#'@seealso
#'For more information, please refer to the package's documentation and the tutorial: \url{https://software.klausjung-lab.de/}.
loop <- function (D, B, inflation = 3, dm)
{
  n <- dim(D)[1]
  d = dim(D)[2]
  if (d==3) dm = matrix(dm, 1, 3)
  index.F <- sort(intersect(as.vector(B$hull), as.vector(B$hull)))
  FENCE <- B$coords[index.F, ]
  MED.MAT <- t(matrix(dm, d, dim(FENCE)[1]))
  FENCE <- MED.MAT + inflation * (FENCE - MED.MAT)
  colnames(FENCE) <- colnames(D)
  convH <- convhulln(FENCE)
  outliers <- rep(0, n)
  for (i in 1:n) {
    Z <- rbind(FENCE, D[i, ])
    convH.Z <- convhulln(Z)
    if (!is.na(match(dim(FENCE)[1] + 1, convH.Z))) {
      outliers[i] <- 1
    }
  }
  LOOP <- D[which(outliers == 0), ]
  convH2 <- convhulln(LOOP)
  return(list(coords.loop = LOOP, hull.loop = convH2, coords.fence = FENCE,
              hull.fence = convH, outliers = outliers))
}
