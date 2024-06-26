#' Plots a gemstone to an interactive graphics device.
#'
#' Only applicable to 3-dimensional data sets. Transparent colors are
#' recommended for outer gemstone of the gemplot. Further graphical
#' parameters can be set using material3d() of the rgl-package.
#'
#' @title Plots a gemstone to an interactive graphics device
#' @param coords Matrix with coordinates of the grid or of data
#'     points that belong to the gemstone, calculated by either
#'     \code{\link{bag}} or \code{\link{loop}}. Each row represents a
#'     grid point and each column represents one dimension.
#' @param hull Matrix with indices of triangles that cover a convex
#'     hull arround the gemstone. Each row represents one triangle
#'     and the indices refer to the rows of coords.
#' @param clr Specifies the color of the gemstone.
#' @return NULL
#' @author Jochen Kruppa, Klaus Jung
#' @references
#' Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The bagplot: a bivariate boxplot. \emph{The American Statistician}, \strong{53(4)}, 382-387. \doi{10.1080/00031305.1999.10474494}

#'
#' Kruppa, J., & Jung, K. (2017). Automated multigroup outlier identification in molecular high-throughput data using bagplots and gemplots. \emph{BMC bioinformatics}, \strong{18(1)}, 1-10. \url{https://link.springer.com/article/10.1186/s12859-017-1645-5}
#' @importFrom rgl triangles3d material3d bg3d points3d text3d spheres3d axes3d
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
#' alpha=c(0.5, 0.5, 0.5, 0.5), smooth=FALSE, specular="black")
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
#'
#' # Example of outlier detection with four principal components.
#' # Attention: calculation is currently time-consuming.
#'
#' set.seed(123)
#' n <- 200
#' x1 <- rnorm(n, 0, 1)
#' x2 <- rnorm(n, 0, 1)
#' x3 <- rnorm(n, 0, 1)
#' x4 <- rnorm(n, 0, 1)
#' D <- data.frame(cbind(x1, x2, x3, x4))
#' D[67,] = c(7, 0, 0, 0)
#'
#' date()
#' G = gridfun(D, 20, 4)
#' G$H = hldepth(D, G, verbose=TRUE)
#' dm = depmed(G)
#' B = bag(D, G)
#' L = loop(D, B, dm=dm)
#' which(L$outliers==1)
#' date()
#' }
#'@seealso
#'For more information, please refer to the package's documentation and the tutorial: \url{https://software.klausjung-lab.de/}.
gem <- function (coords, hull, clr)
{
  for (i in 1:dim(hull)[1]) {
    x <- coords[hull[i, ], 1]
    y <- coords[hull[i, ], 2]
    z <- coords[hull[i, ], 3]
    triangles3d(x, y, z, col = clr)
  }
}
