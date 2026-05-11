#' Isotropytest based on the variogram
#'
#' Returns robust and non-robust isotropy test based on the variogram using a subsampling or a blockpermutation approach for test decision
#'
#' @param data A data frame with three column containing the data (first column) and the grid (second and third column).
#' @param lagmat  A matrix of spatial lags, where each row represent a lag vector specified as \eqn{(lag.x, lag.y)} for which the semivariogram value is to be estimated.
#' @param A A contrast matrix, where each row defines a contrast of the estimated semivariogram evaluated at the lags specified in \code{lagmat}.
#' @param estimator A character string or vector of character strings specifying the estimators to be used. Possible values are \code{"matheron"}, \code{"genton"}, \code{"mcd"}. If \code{"all"} all three estimators are used. See details. Default is \code{"all"}.
#' @param method A character string specifying the mehtod, which is used for estimation of the covariance matrix as well as for p-value calculation. Possible values are \code{"subsampling"} or \code{"blockpermutation"}. See details.
#' @param window.dims A length-two vector specifying to the width and height of the moving windows -given as the number of columns and rows, respectively- used in the subsampling approach.
#' @param var.robust Logical. If \code{TRUE} a robust covariance estimator is used for estimation of the covariance matrix else the classical sample covariance matrix estimator is used. See details. Default is \code{TRUE}.
#' @param edge.sub Logical. If \code{TRUE} a correction for edge effects is used in the supsampling approach as suggested in Guan et. al. (2004). Default is \code{TRUE}.
#' @param c.regularisation Numeric. A value between 0 and 1. Needed for the regularisation of the covariance estimation, if this is singular. See Gierse and Fried (2026); default: 0.1.
#' @param B Numeric. Used permutations in the blockpermutations approach. Default is 1000. description
#'
#' @details
#' This functions implements tests for isotropy based on the idea of Guan et al. (2004). All test are explained in Gierse and Fried (2026).
#'
#' The idea of the test is to compare variogram estimations for lag vectors with the same distance but different directions using a contrast test.
#' For the p-value calculation a subsampling approach with small overlapping subsamples or a blockpermutation approach with non-overlapping block is used.
#'
#' Implemented are three different estimators the Matheron variogram estimator (Matheron, 1962), the robust Genton estimator (Genton, 1998) and the robust MCD.diff estimator (Gierse and Fried, 2025).
#' More information about the implementation can be found in documentations of the functions
#' \code{\link{variogram_est}} and \code{\link{variogram_est_general}}.
#'
#' @return An object of class "isotropyRob" which is a list with an list for each estimator containing the following elements
#' \itemize{
#' \item teststatistic: value of the teststatistic
#' \item p.value: calculated p-value
#' \item cov: a list with the following two elements:
#' \itemize{
#' \item cov: a matrix containing the estimated covariance matrix of the estimated semivariogram vector for all lags of \code{lagmat}
#' \item gamma.sub: a matrix with the semivariogram estimations for each subsample/permutated data set
#' }
#' \item regularization: information if the regularized version of the covariance estimation is used in the teststatistic}
#'
#' @references
#' Genton, M. (1998). Highly robust variogram estimation. \emph{Mathematical Geology}, 30, 213-221. \doi{https://doi.org/10.1023/A:1021728614555}
#'
#' Gierse, J., & Fried, R. (2025). Nonparametric directional variogram estimation in the presence of outlier blocks. \emph{Statistical Papers}, 66(134). \doi{https://doi.org/10.1007/s00362-025-01754-2}
#'
#' Gierse, J. & Fried. R. (2026). EINFÜGEN
#'
#' Guan, Y., Sherman, M. & Calvin, J. A. (2004). A nonparametric test for spatial isotropy using subsampling. \emph{J. Americ. Statist. Assoc.}, 99, 810-821. \doi{https://doi.org/10.1198/016214504000001150}
#'
#' Matheron, G. (1962). Traité de géostatistique appliquée, Tome I. \emph{Mémoires du Bureau de Recherches Géologiques et Miniéres}, no. 14, Editions Technip, 333
#'
#' @seealso \code{\link{isotropy_subsampling}}, \code{\link{isotropy_blockpermutation}}
#'
#' @examples
#' \donttest{
#'  ## Simulation of an anisotrop GRF with an quadratic outlier block
#'  dat <-  simulate_grf(gridsize = c(24, 24), param.variogram = c(1, 6),
#'                       aniso.param = c(pi/4, 2), out.type = "block",
#'                       block.type = "square", amount = 0.1,
#'                       param.outlier = c(0,5), n.it = 1)
#'
#'  dat <- cbind(dat$data, dat$grid)
#'
#'  ## Test for isotropy using the subsample test with 5x5 subsamples
#'  test.sub <- isotropy_test(dat, window.dims = c(5,5), method = "subsampling")
#'
#'  ## Test for isotropy using the blockpermutation test with 6x6 blocks
#'  test.sub <- isotropy_test(dat, window.dims = c(6,6), method = "blockpermutation")
#' }
#'
#' @importFrom robustbase covMcd
#' @importFrom robustbase Qn
#' @export isotropy_test



isotropy_test <- function(data,
                          lagmat = rbind(c(1,0), c(0,1), c(1,1), c(1,-1)),
                          A = rbind(c(1,-1,0,0),c(0,0,1,-1)),
                          estimator = "all",
                          method,
                          window.dims,
                          var.robust = TRUE,
                          edge.sub = TRUE,
                          c.regularisation = 0.1,
                          B = 1000){

  # To avoid error messages caused by incorrect capitalization
  method <- tolower(method)
  estimator <- tolower(estimator)

  if(!is.data.frame(data)) stop("data needs to be a data.frame.")
  if(ncol(data) != 3) stop("data needs to have three columns.")

  if(is.element("all", estimator)) estimator <- c("matheron", "genton", "mcd")
  if(any(!is.element(estimator, c("matheron", "genton", "mcd")))) stop("Invalid estimator. Only the estimators Matheron, Genton, MCD are allowed.")

  if(!any(is.element(c("blockpermutation", "subsampling"), method)))stop("Invalid method. Only Subsampling or Blockpermutation are allowed.")

  if(method != "subsampling" & !is.null(edge.sub)) warning("edge.sub only relevant for the subsampling method. Argument is ignored.")

  if(method == "subsampling") res <- isotropy_subsampling(data, lagmat, A, estimator, window.dims, var.robust, edge.sub, c.regularisation)

  if(method == "blockpermutation") res <- isotropy_blockpermutation(data, lagmat, A, estimator, window.dims, var.robust, c.regularisation, B)

  return(res)
}

