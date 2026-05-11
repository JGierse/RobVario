#' Test of isotropy based on the variogram using subsampling techniques
#'
#' Returns the nonparametric test of istoropy from Guan et. al. (2004) based on the variogram using a subsampling technique for test decision for data on a regular grid. See Guan et. al. (2004) for more details.
#'
#' @param data A dataframe with three column containing the data (first column) and the grid (second and third column).
#' @param lagmat  A matrix of spatial lags, where each row represent a lag vector specified as \eqn{(lag.x, lag.y)} for which the semivariogram value is to be estimated.
#' @param A A contrast matrix, where each row defines a contrast of the estimated semivariogram evaluated at the lags specified in \code{lagmat}.
#' @param estimator A character string or vector of character strings specifying the estimators to be used. Possible values are \code{"matheron"}, \code{"genton"}, \code{"mcd"}. If \code{"all"} all three estimators are used. See details. Default is \code{"all"}.
#' @param window.dims A length-two vector specifying to the width and height of the moving windows -given as the number of columns and rows, respectively- used in the subsampling approach.
#' @param var.robust Logical. If \code{TRUE} a robust covariance estimator is used for estimation of the covariance matrix else the classical sample covariance matrix estimator is used. See details. Default is \code{TRUE}.
#' @param edge Logical. If \code{TRUE} a correction for edge effects is used in the supsampling approach as suggested in Guan et. al. (2004). Default is \code{TRUE}.
#' @param c.regularisation Numeric. A value between 0 and 1. Needed for the regularisation of the covariance estimation, if this is singular. See Gierse and Fried (2026); default: 0.1.
#'
#' @details
#' This functions implements the subsampling test for isotropy of Guan et al. (2004) as well as a version of the test, where the Matheron variogram estimator is
#' replaced by robust variogram estimators. For more information about this see Gierse and Fried (2026).
#'
#' The idea of the test is to compare variogram estimations for lag vectors with the same distance but different directions using a contrast test. A subsampling approach
#' with small overlapping subsamples is used for estimation of the covariance matrix of the variogram estimation and for the p-value calculation.
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
#' \item gamma.sub: a matrix with the semivariogram estimations for each subsample
#' }
#' \item regularisation: information if the regularised version of the covariance estimation is used in the teststatistic}
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
#' @seealso \code{\link{isotropy_test}}, \code{\link{isotropy_blockpermutation}}
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
#'  ## Test for isotropy using 5x5 subsamples
#'  test <- isotropy_subsampling(dat, window.dims = c(5,5))
#' }
#'
#' @importFrom robustbase covMcd
#' @importFrom robustbase Qn
#' @importFrom stats cov
#' @export isotropy_subsampling

isotropy_subsampling <- function(data,
                                 lagmat = rbind(c(1,0), c(0,1), c(1,1), c(1,-1)),
                                 A = rbind(c(1,-1,0,0),c(0,0,1,-1)),
                                 estimator = "all",
                                 window.dims,
                                 var.robust = TRUE,
                                 edge = TRUE,
                                 c.regularisation = 0.1){

  # To avoid error messages caused by incorrect capitalization
  estimator <- tolower(estimator)

  if(is.element("all", estimator)) estimator <- c("matheron", "genton", "mcd")
  if(any(!is.element(estimator, c("matheron", "genton", "mcd")))) stop("Invalid estimator. Only the estimators Matheron, Genton, MCD are allowed.")

  ny <- length(unique(data[,3]))
  nx <- length(unique(data[,2]))
  n <- nx * ny
  n.lags <- nrow(lagmat)
  n.sub <- prod(window.dims)

  rawdata <- .differences.lags(data)
  rawdata <- .raw.data.sub(rawdata, lagmat)

  res <- list()
  ind.res <- 1

  if(edge == TRUE) rawdata <- .raw.data.edge(rawdata, lagmat)

  rawdata <- as.data.frame(rawdata)
  rawdata <- cbind(rawdata, "coord" = paste0(rawdata[,"x.coord"], ".", rawdata[,"y.coord"]))
  rawdata <- cbind(rawdata, "lag" = paste0(rawdata[,"x.lag"], ".", rawdata[,"y.lag"]))

  normalize.lag <- apply(rawdata, 1, function(x){
    x <- c(as.numeric(x["x.lag"]), as.numeric(x["y.lag"]))
    .normalize.row(x)
  })
  rawdata$lag.norm <- paste0(normalize.lag[1,], ".", normalize.lag[2,])

  # Subsampling einmal vorab durchführen, um duplikation bei den verschiedenen SChätzmethoden zu vermeiden
  sub.rawdata.list <- .subsampling(rawdata, window.dims)

  if(is.element("matheron", estimator)){

    gamma.matheron <- .matheron.test(rawdata, lagmat)

    cov.matheron <- .estimation.sigma.subsampling(sub.rawdata.list, lagmat, window.dims, estimator = "matheron", robust = var.robust)

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
    T.sub.matheron <- c()
    for(i in 1:nrow(cov.matheron$gamma.sub)){
      gh.matheron <- matrix(cov.matheron$gamma.sub[i,], nrow = n.lags, ncol = 1)
      t.matheron <- n.sub*t(A%*%gh.matheron) %*% solve(cov.A.matheron) %*% (A%*%gh.matheron)
      T.sub.matheron <-  c(T.sub.matheron, t.matheron)
    }
    pvalue.matheron <- sum(T.sub.matheron >= T.matheron)/nrow(cov.matheron$gamma.sub)

    res[[ind.res]] <- list("teststatistic" = T.matheron, "p.value" = pvalue.matheron, "cov" = cov.matheron, "regularization" = reg.matheron, "variogram.est" = gamma.matheron[,"gamma.hat"])
    names(res)[ind.res] <- "matheron"
    ind.res <- ind.res + 1
  }

  if(is.element("genton", estimator)){

    gamma.genton <- .genton.test(rawdata, lagmat)

    cov.genton <- .estimation.sigma.subsampling(sub.rawdata.list, lagmat, window.dims, estimator = "genton", robust = var.robust)

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
    T.sub.genton <- c()
    for(i in 1:nrow(cov.genton$gamma.sub)){
      gh.genton <- matrix(cov.genton$gamma.sub[i,], nrow = n.lags, ncol = 1)
      t.genton <- n.sub*t(A%*%gh.genton) %*% solve(cov.A.genton) %*% (A%*%gh.genton)
      T.sub.genton <-  c(T.sub.genton, t.genton)
    }
    pvalue.genton <- sum(T.sub.genton >= T.genton)/nrow(cov.genton$gamma.sub)

    res[[ind.res]] <- list("teststatistic" = T.genton, "p.value" = pvalue.genton, "cov" = cov.genton, "regularization" = reg.genton, "variogram.est" = gamma.genton[,"gamma.hat"])
    names(res)[ind.res] <- "genton"
    ind.res <- ind.res + 1
  }

  if(is.element("mcd", estimator)){

    vec.data <- .rawdata.to.vec(rawdata)

    gamma.mcd <- .mcd.test(vec.data, lagmat)

    cov.mcd <- .estimation.sigma.subsampling(sub.rawdata.list, lagmat, window.dims, estimator = "mcd", robust = var.robust)

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
    T.sub.mcd <- c()
    for(i in 1:nrow(cov.mcd$gamma.sub)){
      gh.mcd <- matrix(cov.mcd$gamma.sub[i,], nrow = n.lags, ncol = 1)
      t.mcd <- n.sub*t(A%*%gh.mcd) %*% solve(cov.A.mcd) %*% (A%*%gh.mcd)
      T.sub.mcd <-  c(T.sub.mcd, t.mcd)
    }
    pvalue.mcd <- sum(T.sub.mcd >= T.mcd)/nrow(cov.mcd$gamma.sub)

    res[[ind.res]] <- list("teststatistic" = T.mcd, "p.value" = pvalue.mcd, "cov" = cov.mcd, "regularization" = reg.mcd, "variogram.est" = gamma.mcd[,"gamma.hat"])
    names(res)[ind.res] <- "MCD"
    ind.res <- ind.res + 1
  }

  class(res) <- "isotropyRob"
  return(res)
}


.raw.data.edge <- function(rawdata,
                           lagmat){

  nlags <-  dim(lagmat)[1]

  index.set <-  c()
  clag <-  lagmat[1,]
  good <-  which((rawdata[,"x.lag"] == clag[1] & rawdata[,"y.lag"] == clag[2])|(rawdata[,"x.lag"] == -clag[1] & rawdata[,"y.lag"] == -clag[2]))
  index.set <-  rawdata[good,"index"]

  for(i in 2:nlags){
    clag <-  lagmat[i,]
    good <-  which((rawdata[,"x.lag"] == clag[1] & rawdata[,"y.lag"] == clag[2])|(rawdata[,"x.lag"] == -clag[1] & rawdata[,"y.lag"] == -clag[2]))
    index.tmp <-  rawdata[good,"index"]
    index.set <-  intersect(index.set, index.tmp)
  }

  good <-  which(rawdata[,"index"] %in% index.set)
  rawdata <-  rawdata[good,]

  return(rawdata)

}

.estimation.sigma.subsampling <- function(sub.rawdata.list,
                                          lagmat,
                                          window.dims,
                                          estimator,
                                          robust){

  n.subsample <- prod(window.dims)

  nlags <-  dim(lagmat)[1]

  if(estimator == "matheron"){
    gamma.sub.list <- lapply(sub.rawdata.list, .matheron.test, lagmat = lagmat)
  }
  if(estimator == "genton"){
    gamma.sub.list <- lapply(sub.rawdata.list, .genton.test, lagmat = lagmat)
  }

  if(estimator == "mcd"){

    vec.sub <- lapply(sub.rawdata.list, .rawdata.to.vec)

    gamma.sub.list <- lapply(vec.sub, .mcd.test, lagmat = lagmat)

  }

  gamma.sub.list <- lapply(gamma.sub.list, function(x) x[,"gamma.hat"])
  gamma.sub <- do.call(rbind, gamma.sub.list)
  colnames(gamma.sub) <- paste(lagmat[,1], lagmat[,2])

  # Variance estimation
  if(robust == TRUE){
    cov.sub <- suppressWarnings(covMcd(gamma.sub)$cov)
  } else{
    cov.sub <- cov(gamma.sub)
  }

  res <- list("cov" = cov.sub, "gamma.sub" = gamma.sub)

  return(res)

}

.subsampling <- function(rawdata,
                         window.dims){

  min.x <- min(rawdata[,"x.coord"])
  max.x <- max(rawdata[,"x.coord"])
  min.y <- min(rawdata[,"y.coord"])
  max.y <- max(rawdata[,"y.coord"])

  blk.x <- min.x:(max.x-window.dims[1]+1)
  blk.y <- min.y:(max.y-window.dims[1]+1)

  blk.coords <- expand.grid(blk.x, blk.y)
  blk.coords <- cbind(blk.coords[,1], blk.coords[,2])

  bw <- window.dims[1]
  bh <- window.dims[2]

  block.data <-  list()
  n.blks <- dim(blk.coords)[1]
  for(i in 1:n.blks){
    cblock <- blk.coords[i,]
    inb <- which(rawdata[,"x.coord"] >= cblock[1] & rawdata[,"x.coord"] < cblock[1]+bw & rawdata[,"y.coord"] >= cblock[2] & rawdata[,"y.coord"] < cblock[2]+bh)
    block.data[[i]] <- rawdata[inb,]
  }

  return(block.data)
}
