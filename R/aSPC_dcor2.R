#' An Adaptive Sum of Powered Correlation Test (aSPC) with dcor
#'
#' @param df1, first matrix
#' @param df2, second matrix
#' @param pow, power integer candidates, default c(1:8, Inf)
#' @param B, number of permutations to calculate a P-value
#' @return the P-values of SPC and aSPC tests
#' @references Xu Z., Pan W. An adaptive and powerful test for two groups of variables with high dimension
#'
#'
#'
#' @importFrom energy dcor
#' @keywords internal



aSPC_dcor2 = function(df1, df2, pow = pow, B = B){
  # test parameter
  # df1 = X; df2 = Y; pow = c(1:8, Inf); B=2


  X = scale(df1)
  Y = scale(df2)

  n = nrow(X) ## number of subjects
  ### X and Y has to be standardize before input
  ls_X = lapply(1:ncol(X), function(x) get_doublyCenterDist(X[,x]))
  ls_distCenterMat_X = lapply(1:ncol(X), function(x) ls_X[[x]]$mat )
  ls_distvariance_X = lapply(1:ncol(X), function(x) ls_X[[x]]$dVar )

  ls_Y = lapply(1:ncol(Y), function(x) get_doublyCenterDist(Y[,x]))
  ls_distCenterMat_Y = lapply(1:ncol(Y), function(x) ls_Y[[x]]$mat )
  ls_distvariance_Y = lapply(1:ncol(Y), function(x) ls_Y[[x]]$dVar )

  # a = proc.time()
  mat_obs = dcor_list2(ls_distCenterMat_X, ls_distCenterMat_Y,
                      ls_distvariance_X, ls_distvariance_Y)
  # proc.time() - a

  ## obs statistics
  T_obs = rep(NA,length(pow))


  for(k in 1:length(pow)){
    if(pow[k]<Inf) T_obs[k] = sum(mat_obs^pow[k]) else T_obs[k] = max(abs(mat_obs))
  }

  pPerm0 = rep(NA,length(pow))

  T0s = matrix(nrow = B, ncol = length(pow))
  for(i in 1: B){
    index = sample(nrow(X))
    ls_distCenterMat_X_sample = lapply(1:length(ls_distCenterMat_X), function(x) ls_distCenterMat_X[[x]][index, index])
    mat0 = dcor_list2(ls_distCenterMat_X_sample, ls_distCenterMat_Y,
                      ls_distvariance_X, ls_distvariance_Y)
    for(j in 1:length(pow)){
      if(pow[j] < Inf) T0s[i,j] = sum(mat0^pow[j]) else T0s[i,j] = max(abs(mat0))
    }
    #if(show_b == T) print(i)
  }

  for(j in 1:length(pow)){
    pPerm0[j] = round((sum(abs(T_obs[j])<=abs(T0s[1:(B-1),j]))+1)/(B), digits=8)
    P0s = (B-rank(abs(T0s[,j]))+1)/(B)
    if (j==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]

  }

  Paspu<-(sum(minp0<=min(pPerm0))+1)/(B+1)
  pvs <- c(pPerm0, Paspu)
  names(pvs) = c(paste0("SPC.",pow),"aSPC")
  return(pvs)


}
