#' my program to get distance correlation
#'
#' @param A, doubly centered matrix from get_doublyCenterDist(x, index)$mat
#' @param B, doubly centered matrix from get_doublyCenterDist(x, index)$mat
#' @param dVarX, doubly centered matrix from get_doublyCenterDist(x, index)$dVar
#' @param dVarY, doubly centered matrix from get_doublyCenterDist(x, index)$dVar
#' @return the P-values of SPC and aSPC tests
#' @references Xu Z., Pan W. An adaptive and powerful test for two groups of variables with high dimension
#' @importFrom energy dcor
#' @keywords internal

My_DCOR = function (A, B, dVarX, dVarY)
{

  dCov <- sqrt(mean(A * B))
  V <- sqrt(dVarX * dVarY)
  if (V > 0) dCor <- dCov/V else dCor <- 0

  return(dCor = dCor)
}
