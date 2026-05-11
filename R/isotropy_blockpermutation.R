#' Isotropytest based on the variogram using blockpermutation
#'
#' Returns robust and non-robust isotropy test based on the variogram using blockpermutation approach for test decision
#'
#' @param data A data frame with three column containing the data (first column) and the grid (second and third column).
#' @param lagmat  A matrix of spatial lags, where each row represent a lag vector specified as \eqn{(lag.x, lag.y)} for which the semivariogram value is to be estimated.
#' @param A A contrast matrix, where each row defines a contrast of the estimated semivariogram evaluated at the lags specified in \code{lagmat}.
#' @param estimator A character string or vector of character strings specifying the estimators to be used. Possible values are \code{"matheron"}, \code{"genton"}, \code{"mcd"}. If \code{"all"} all three estimators are used. See details. Default is \code{"all"}.
#' @param window.dims A length-two vector specifying to the width and height of the subblocks -given as the number of columns and rows, respectively- used in the blockpermutation approach. See details.
#' @param var.robust Logical. If \code{TRUE} a robust covariance estimator is used for estimation of the covaraince matrix else the classical sample covariance matrix estimator is used. See details. Default is \code{TRUE}.
#' @param c.regularisation Numeric. A value between 0 and 1. Needed for the regularisation of the covariance estimation, if this is singular. See Gierse and Fried (2026); default: 0.1.
#' @param B Numeric. Number of permutations in the blockpermutations approach. Default is 1000.
#'
#' @details
#' This functions implements the blockpermutation test for isotropy proposed by Gierse and Fried (2026).  The test is based on the test idea of Guan et al. (2004).
#'
#' The idea of the test is to compare variogram estimations for lag vectors with the same distance but different directions using a contrast test. A blockpermutation approach
#' with small non-overlapping block is used for estimation of the covariance matrix of the variogram estimation and for the p-value calculation. For more details see  Gierse and Fried (2026).
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
#' \item gamma.sub: a matrix with the semivariogram estimations for each permutated data set
#' }
#' \item regularization: information if the regularized version of the covariance estimation is used in the teststatistic}
#'
#'@references
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
#' @seealso \code{\link{isotropy_test}}, \code{\link{isotropy_subsampling}}
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
#'  ## Test for isotropy using 6x6 blocks
#'  test <- isotropy_blockpermutation(dat, window.dims = c(6,6))
#' }
#'
#' @importFrom robustbase covMcd
#' @importFrom robustbase Qn
#' @importFrom stats cov
#' @export isotropy_blockpermutation

isotropy_blockpermutation <- function(data,
                                      lagmat = rbind(c(1,0), c(0,1), c(1,1), c(1,-1)),
                                      A = rbind(c(1,-1,0,0),c(0,0,1,-1)),
                                      estimator = "all",
                                      window.dims = c(6,6),
                                      var.robust = TRUE,
                                      c.regularisation = 0.1,
                                      B = 1000){

  # To avoid error messages caused by incorrect capitalization
  estimator <- tolower(estimator)

  if(is.element("all", estimator)) estimator <- c("matheron", "genton", "mcd")
  if(any(!is.element(estimator, c("matheron", "genton", "mcd")))) stop("Invalid estimator. Only the estimators Matheron, Genton, MCD are allowed.")

  ny <- length(unique(data[,3]))
  nx <- length(unique(data[,2]))
  n <- nx * ny
  n.lags <- nrow(lagmat)
  n.block <- prod(window.dims)

  rawdata <- .differences.lags(data)
  rawdata <- .raw.data.sub(rawdata, lagmat)

  b.x <- lapply(seq(1,nx,window.dims[1]), function(l) l:(l+window.dims[1] - 1))
  b.y <- lapply(seq(1,ny,window.dims[2]), function(l) l:(l+window.dims[2] - 1))

  block.grid <- expand.grid(1:length(b.x), 1:length(b.y))

  blocks <- lapply(1:nrow(block.grid), function(l) data[is.element(data[,2], b.x[[block.grid[l,1]]]) &
                                                          is.element(data[,3], b.y[[block.grid[l,2]]]),])

  blocks.rawdata <- lapply(blocks, .differences.lags)
  blocks.rawdata <- lapply(blocks.rawdata, .raw.data.sub, lagmat = lagmat)
  blocks.rawdata <- lapply(blocks.rawdata, as.data.frame)
  blocks.rawdata <- lapply(blocks.rawdata, function(x) cbind(x, "coord" = paste0(x[,"x.coord"], ".", x[,"y.coord"])))
  blocks.rawdata <- lapply(blocks.rawdata, function(x) cbind(x, "lag" = paste0(x[,"x.lag"], ".", x[,"y.lag"])))

  # Blockpermutierte Daten erstellen, um Wiederholung für jeden Schätzer zu vermeiden!
  block.rawdata.list <- .blockpermutation(blocks.rawdata, B)

  res <- list()
  ind.res <- 1

  if(is.element("matheron", estimator)){

    gamma.matheron <- .matheron.test(rawdata, lagmat)

    cov.matheron <- .estimation.sigma.blockpermutation(block.rawdata.list, lagmat, B, estimator = "matheron", robust = var.robust)

    # regularization if A%*%sigma.hat%*% singular
    cov.A.matheron <- A%*% cov.matheron$cov%*%t(A)
    reg.matheron = FALSE
    if(.is.almost.zero(det(cov.A.matheron))){
      cov.A.matheron <- cov.A.matheron +  c.regularisation * (sum(diag(cov.A.matheron))/nrow(A)) * diag(nrow(A))
      reg.matheron <- TRUE
    }

    # Calculation Teststatistik
    T.matheron <- n * t(A%*%gamma.matheron[,"gamma.hat"])%*%solve(cov.A.matheron)%*%(A%*%gamma.matheron[,"gamma.hat"])
    T.matheron <- as.numeric(T.matheron)

    # pvalue
    T.block.matheron <- c()
    for(i in 1:B){
      gh.matheron <- matrix(cov.matheron$gamma.block[i,], nrow = n.lags, ncol = 1)
      t.matheron <- n*t(A%*%gh.matheron) %*% solve(cov.A.matheron) %*% (A%*%gh.matheron)
      T.block.matheron <-  c(T.block.matheron, t.matheron)
    }
    pvalue.matheron <- sum(T.block.matheron >= T.matheron)/B

    res[[ind.res]] <- list("teststatistic" = T.matheron, "p.value" = pvalue.matheron, "cov" = cov.matheron, "regularization" = reg.matheron, "variogram.est" = gamma.matheron[,"gamma.hat"])
    names(res)[ind.res] <- "matheron"
    ind.res <- ind.res + 1
  }

  if(is.element("genton", estimator)){

    gamma.genton <- .genton.test(rawdata, lagmat)

    cov.genton <- .estimation.sigma.blockpermutation(block.rawdata.list, lagmat, B, estimator = "genton", robust = var.robust)

    # regularization if A%*%sigma.hat%*% singular
    cov.A.genton <- A%*% cov.genton$cov%*%t(A)
    reg.genton = FALSE
    if(.is.almost.zero(det(cov.A.genton))){

      cov.A.genton <- cov.A.genton +  c.regularisation * (sum(diag(cov.A.genton))/nrow(A)) * diag(nrow(A))
      reg.genton <- TRUE
    }

    # Calculation Teststatistik
    T.genton <- n * t(A%*%gamma.genton[,"gamma.hat"])%*%solve(cov.A.genton)%*%(A%*%gamma.genton[,"gamma.hat"])
    T.genton <- as.numeric(T.genton)

    # pvalue
    T.block.genton <- c()
    for(i in 1:B){
      gh.genton <- matrix(cov.genton$gamma.block[i,], nrow = n.lags, ncol = 1)
      t.genton <- n*t(A%*%gh.genton) %*% solve(cov.A.genton) %*% (A%*%gh.genton)
      T.block.genton <-  c(T.block.genton, t.genton)
    }
    pvalue.genton <- sum(T.block.genton >= T.genton)/B

    res[[ind.res]] <- list("teststatistic" = T.genton, "p.value" = pvalue.genton, "cov" = cov.genton, "regularization" = reg.genton, "variogram.est" = gamma.genton[,"gamma.hat"])
    names(res)[ind.res] <- "genton"
    ind.res <- ind.res + 1
  }

  if(is.element("mcd", estimator)){

    rawdata <- as.data.frame(rawdata)
    rawdata <- cbind(rawdata, "coord" = paste0(rawdata[,"x.coord"], ".", rawdata[,"y.coord"]))
    rawdata <- cbind(rawdata, "lag" = paste0(rawdata[,"x.lag"], ".", rawdata[,"y.lag"]))


    normalize.lag <- apply(rawdata, 1, function(x){
      x <- c(as.numeric(x["x.lag"]), as.numeric(x["y.lag"]))
      .normalize.row(x)
    })
    rawdata <- as.data.frame(rawdata)
    rawdata$lag.norm <- paste0(normalize.lag[1,], ".", normalize.lag[2,])

    vec.data <- .rawdata.to.vec(rawdata)

    gamma.mcd <- .mcd.test(vec.data, lagmat)

    cov.mcd <- .estimation.sigma.blockpermutation(block.rawdata.list, lagmat, B, estimator = "mcd", robust = var.robust)

    # regularization if A%*%sigma.hat%*% singular
    cov.A.mcd <- A%*% cov.mcd$cov%*%t(A)
    reg.mcd <- FALSE
    if(.is.almost.zero(det(cov.A.mcd))){

      cov.A.mcd <- cov.A.mcd +  c.regularisation * (sum(diag(cov.A.mcd))/nrow(A)) * diag(nrow(A))
      reg.mcd <- TRUE
    }

    # Calculation Teststatistik
    T.mcd <- n * t(A%*%gamma.mcd[,"gamma.hat"])%*%solve(cov.A.mcd)%*%(A%*%gamma.mcd[,"gamma.hat"])
    T.mcd <- as.numeric(T.mcd)

    # pvalue
    T.block.mcd <- c()
    for(i in 1:nrow(cov.mcd$gamma.block)){
      gh.mcd <- matrix(cov.mcd$gamma.block[i,], nrow = n.lags, ncol = 1)
      t.mcd <- n*t(A%*%gh.mcd) %*% solve(cov.A.mcd) %*% (A%*%gh.mcd)
      T.block.mcd <-  c(T.block.mcd, t.mcd)
    }
    pvalue.mcd <- sum(T.block.mcd >= T.mcd)/nrow(cov.mcd$gamma.block)

    res[[ind.res]] <- list("teststatistic" = T.mcd, "p.value" = pvalue.mcd, "cov" = cov.mcd, "regularization" = reg.mcd, "variogram.est" = gamma.mcd[,"gamma.hat"])
    names(res)[ind.res] <- "MCD"
    ind.res <- ind.res + 1
  }

  class(res) <- "isotropyRob"
  return(res)
}


.estimation.sigma.blockpermutation <- function(block.rawdata.list,
                                               lagmat,
                                               B,
                                               estimator,
                                               robust){


  nlags <-  dim(lagmat)[1]

  if(estimator == "matheron"){
    gamma.block.list <- lapply(block.rawdata.list, .matheron.test, lagmat = lagmat)
  }
  if(estimator == "genton"){
    gamma.block.list <- lapply(block.rawdata.list, .genton.test, lagmat = lagmat)
  }

  if(estimator == "mcd"){
    vec.block <- lapply(block.rawdata.list, .rawdata.to.vec)

    gamma.block.list <- lapply(vec.block, .mcd.test, lagmat = lagmat)
  }

  gamma.block.list <- lapply(gamma.block.list, function(x) x[,"gamma.hat"])
  gamma.block <- do.call(rbind, gamma.block.list)
  colnames(gamma.block) <- paste(lagmat[,1], lagmat[,2])


  # Variance estimation
  if(robust == TRUE){
    cov.block <- suppressWarnings(covMcd(gamma.block)$cov)
  } else{
    cov.block <- cov(gamma.block)
  }

  res <- list("cov" = cov.block, "gamma.block" = gamma.block)

  return(res)

}

.blockpermutation <- function(block.rawdata,
                              B){

  blocks.perm <- list()
  for(b in 1:B){
    perm <- block.rawdata
    move <- sample(c(T,F), length(block.rawdata), replace = TRUE)

    for(m in 1:length(move)){
      if(move[m] == TRUE){
        old <- block.rawdata[[m]]
        old.lag <- as.matrix(old[,c("x.lag", "y.lag")])

        new.lag <- apply(old.lag, 1, function(y){
          y[c("x.lag", "y.lag")] <- y[c("y.lag", "x.lag")]
          if(all(y[c("x.lag", "y.lag")] != 0)){
            y["y.lag"] = -y["y.lag"]
            if(y["x.lag"] < 0 & y["y.lag"] < 0){
              y[c("x.lag", "y.lag")] <- -y[c("x.lag", "y.lag")]
            }
          }
          return(y)
        })
        new <- old
        new[, c("x.lag", "y.lag")] <- t(new.lag)
        new$lag <- paste0(new$x.lag, ".", new$y.lag)
        perm[[m]] <- new
      }
    }
    perm <- do.call(rbind, perm)
    normalize.lag <- apply(perm, 1, function(x){
      x <- c(as.numeric(x["x.lag"]), as.numeric(x["y.lag"]))
      .normalize.row(x)
    })
    perm$lag.norm <- paste0(normalize.lag[1,], ".", normalize.lag[2,])

    blocks.perm[[b]] <- perm
  }

  return(blocks.perm)
}
