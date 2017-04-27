#' calculate the correlation matrix by dcor() in package energy
#'
#' @param df1, first matrix
#' @param df2, second matrix
#' @return a p by q correlation matrix
#' @importFrom energy dcor
#' @keywords internal

dcor_mat <- function(df1, df2){
  p <- ncol(df1)
  q <- ncol(df2)

  each_row <- function(y) {
    apply(df2, 2, function(x) dcor(x,y))
  }
  cor_mat <- t(apply(df1, 2, each_row))
  return(cor_mat)
}



