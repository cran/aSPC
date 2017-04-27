#' calculate the correlation matrix
#'
#' @param ls_distCenterMat_X, first list of doubly centered distance matrices for p columns in df1
#' @param ls_distCenterMat_Y, second list of doubly centered distance matrices for p columns in df1
#' @param ls_distvariance_X, first list of distance variance in df1
#' @param ls_distvariance_Y, second list of distance variance in df2
#'
#' @return a p by q correlation matrix
#' @export
#' @importFrom energy dcor
#' @keywords internal

dcor_list2 <- function(ls_distCenterMat_X, ls_distCenterMat_Y,
                       ls_distvariance_X, ls_distvariance_Y){

  each_row <- function(y) {
    sapply(1:length(ls_distCenterMat_Y), function(x) My_DCOR(ls_distCenterMat_Y[[x]],
                                                             ls_distCenterMat_X[[y]],
                                                             ls_distvariance_Y[[x]],
                                                             ls_distvariance_X[[y]]))
  }
  cor_mat <- t(sapply(1:length(ls_distCenterMat_X), each_row))
  return(cor_mat)
}
