#' Compute the Bhattacharyya Coefficient (BC) between two 1d Gaussian 
#' distributions.
#' 
#' @description
#' As the BC between two multivariate Gaussian distributions with 
#' diagonal covariance matrices (independent dimensions) is the product of
#' the BC across dimensions. This function can be used to vectorize the 
#' computation of BC between multivariate Gaussian distributions 
#'
#' @param mu1 Mean of the first distribution
#' @param mu2 Mean of the second distribution
#' @param var1 Variance of the first distribution
#' @param var2 Variance of the second distribution
#'
#' @return The BC between the 2 distributions
#' @export
#'
#' @examples bhattacharyya_1d(1, 0, 2, 1)
bhattacharyya_1d <- function(mu1, mu2, var1, var2) {
#  stopifnot("variance must be strictly positive" = {var1 >0; var2>0})
  sqrt(2*sqrt(var1*var2)/(var1+var2))*exp(-.25*(mu1-mu2)**2/(var1+var2))
}



kullback_1d <- function(mu1, mu2, var1, var2) {
#  stopifnot("variance must be strictly positive" = {var1 >0; var2>0})
  log(sqrt(var2/var1)) + .5*var1/var2 + .5*(mu1-mu2)**2/var2 - .5
}


coord2std <- function(coord, center, scale = 1e4) {
  R = 6377726
  (pi*R*(coord - center)/180)/scale
}

std2coord <- function(std, center, scale = 1e4) {
  R = 6377726
  scale*std*180/(pi*R) + center
}