#' An adaptive sum of powered correlation test (aSPC)
#' for association between two random vectors
#'
#' @param df1, first sample matrix
#' @param df2, second sample matrix
#' @param pow, power integer candidates, default c(1:8, Inf)
#' @param B, number of permutations to calculate a P-value. Default is 100.
#' @param Z.transform, whether to do Fisher's z-transformation on Pearson correlation, default is TRUE.
#' @param method, one of "pearson", "spearman", or "dcor". Default is "pearson".
#' @return the P-values of SPC and aSPC tests
#' @references Xu Z., Pan W. 2017. Adaptive testing for association between two random vectors in moderate to high dimensions. Submitted to Genetic Epidemiology
#' @references Kim J., Zhang Y., Pan W. Powerful and Adaptive Testing for Multi-trait and Multi-SNP Associa-tions with GWAS and Sequencing Data. Genetics, 2016, 203(2): 715-731.
#' @examples
#' library(mvtnorm)
#' sigma = diag(0.9, 10) + 0.1
#' n = 50 # sample size
#' Z = rmvnorm(n=n, mean=rep(0,10), sigma=sigma)
#' X = rmvnorm(n=n, mean=rep(0,15), sigma=diag(1, 15))
#' Y = rmvnorm(n=n, mean=rep(0,15), sigma=diag(1, 15))
#' X = as.data.frame(cbind(Z[,1:5], X))
#' Y = as.data.frame(cbind(Z[,6:10], Y))
#' set.seed(123) # to ensure we can replicate the permutation P-value

#' p = 2; q = 2; n=50
#' X = rmvnorm(n=n, mean=rep(0,p), sigma=diag(1, p))
#' Y = rmvnorm(n=n, mean=rep(0,q), sigma=diag(1, q))
#' a = proc.time()
#' aSPC(X, Y, pow = c(1:8, Inf), B = 1000, method = "pearson")
#' proc.time() - a
#'
#' #' a = proc.time()
#' aSPC(X, Y, pow = c(1:8, Inf), B = 1000, method = "spearman")
#' proc.time() - a
#'
#' a = proc.time()
#' aSPC(X, Y, pow = c(1:8, Inf), B = 500, method = "dcor")
#' proc.time() - a

#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats cor

aSPC = function(df1, df2, pow = c(1:6, Inf), B = 100, Z.transform = TRUE, method = "pearson"){
  if(method == "dcor") return(aSPC_dcor2(df1, df2, pow = pow, B = B) )
  # test
  # pow = c(1:6, Inf); B = 100; Z.transform = TRUE; method = "dcor"

  X = scale(df1)
  Y = scale(df2)
  a = 10
  n = nrow(X) ## number of subjects
  ### X and Y has to be standardize before input
  if(method == "pearson") mat_obs = t(X)%*%Y/n else if(method=="dcor")
    mat_obs = dcor_mat(X, Y) else
    mat_obs = cor(X, Y, method = method)
  if(Z.transform == TRUE & method == "pearson") mat_obs = a*sqrt(n-3)*1/2*log((1+mat_obs)/(1-mat_obs))
  ## obs statistics
  T_obs = rep(NA,length(pow))


  for(k in 1:length(pow)){
    if(pow[k]<Inf) T_obs[k] = sum(mat_obs^pow[k]) else T_obs[k] = max(abs(mat_obs))
  }

  pPerm0 = rep(NA,length(pow))
  T0s = numeric(B)
  s <- sample(1:10^5,1)

  for(j in 1:length(pow)){
    set.seed(s) # to ensure the same samples are drawn for each pow
    for(i in 1:B){
      X_sample = X[sample(nrow(X)),]
      if(method == "pearson") mat0 = t(X_sample) %*% Y/n else if(method=="dcor")
        mat0 = dcor_mat(X_sample, Y) else
        mat0 = cor(X_sample, Y, method = method)
      if(Z.transform == TRUE & method == "pearson") mat0 = a*sqrt(n-3)*1/2*log((1+mat0)/(1-mat0))

      if(pow[j] < Inf) T0s[i] = sum(mat0^pow[j]) else T0s[i] = max(abs(mat0))
    }
    pPerm0[j] = round((sum(abs(T_obs[j])<=abs(T0s[1:(B-1)]))+1)/(B), digits=8)
    P0s = (B-rank(abs(T0s))+1)/(B)
    if (j==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
  }


  Paspu<-(sum(minp0<=min(pPerm0))+1)/(B+1)
  pvs <- c(pPerm0, Paspu)
  names(pvs) = c(paste0("SPC.",pow),"aSPC")
  return(pvs)
}

