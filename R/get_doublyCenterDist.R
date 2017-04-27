#' my program to get doubly center distance matrix
#'
#' @param x, distance matrix or data matrix
#' @param index, power to the Euclidean norm
#' @return the P-values of SPC and aSPC tests
#' @references Xu Z., Pan W. An adaptive and powerful test for two groups of variables with high dimension
#' @importFrom energy dcor
#' @importFrom stats dist
#' @keywords internal

get_doublyCenterDist = function (x, index = 1)
{
  if (!(class(x) == "dist"))
    x <- dist(x)

  x <- as.matrix(x)
  n <- nrow(x)

  if (!(all(is.finite(c(x)))))
    stop("Data contains missing or infinite values")
  if (index < 0 || index > 2) {
    warning("index must be in [0,2), using default index=1")
    index = 1
  }
  Akl <- function(x) {
    d <- as.matrix(x)^index
    m <- rowMeans(d)
    M <- mean(d)
    a <- sweep(d, 1, m)
    b <- sweep(a, 2, m)
    return(b + M)
  }
  A <- Akl(x)
  dVarX <- sqrt(mean(A * A))
  return(list(mat = A, dVar = dVarX))
}
