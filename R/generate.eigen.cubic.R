#' @title Grid Points and Eigendecomposition of Reproducing Kernal for Cubic Splines
#' @description Generates eigendecomposition from grid points from an interval.
#' @param N The number of grid points. It should be an integer. The default value is 1000.
#' @param a The lower limit of the interval used for grid points. The default value is 0.
#' @param b The upper limit of the interval used for grid points. The default value is 1.
#' @export
#' @import assist
#' @return It returns (and saves) a list of following components:
#' \item{e}{A list of two elements of "values" and "vectors", which refer respectively, the eigenvalues and eigenfunctions of reproducing kernel for cubic splines at the pre-selected grid points.}
#' \item{xg}{A vector of grid points.}
#' @examples
#' \dontrun{
#' eigen_res <- generate.eigen.cubic(N=1000,a=2,b=3)
#' }
#'


generate.eigen.cubic <- function(N=1000,a=0,b=1){
  xg <- seq(a,b,length.out = N)
  Sigma <- cubic(xg)
  e <- eigen(Sigma)
  eigenMN <- list(e=e,xg=xg)
  save(eigenMN,file=paste("eigenM_",a,"&",b,"_",N,".rda",sep=""))
  return(eigenMN)
}
