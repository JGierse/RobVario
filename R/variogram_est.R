#' Nonparametric directional variogram estimation
#'
#' Returns robust and non-robust directional variogram estimations in the four main directions south-north, east-west, southeast-northwest and southwest-northeast for data on a two-dimensional regular grid such as image data.
#'
#' @param data A dataframe with three columns containing the data (first column) and the grid (second and third column).
#' @param hmax A vector with the maximal lag to be estimated for each entry in \code{direction}. If it has length 1 the same hmax is used for all directions.
#' @param direction A vector containing the directions to be estimated. Possible values are \code{"E-W"} (east-west), \code{"S-N"} (south-north), \code{"SW-NE"} (southwest-northeast) and \code{"SE-NW"} (southeast-northwest). Default is all four directions.
#' @param estimator A character string or vector of character strings specifying the estimators to be used. Possible estimators are \code{"matheron"}, \code{"genton"}, \code{"mcd.org"} or \code{"mcd.diff"}. If \code{"all"} all four estimators are used. See details; Default is "all".
#' @param cp A vector containing possible finite sample correction factors for each direction for the MCD based variogram estimators. If it has length 1 the same correction factor is used for all directions.
#' @param reweighting Logical indicating if the reweighted version of the MCD based variogram estimators should be use. Default is \code{TRUE}.
#' @param ... Optional argument for the function \code{covMcd} of the package \pkg{robustbase} used for the MCD based variogram estimators
#'
#' @details
#' The variogram is a measure for the spatial dependency. For a given spatial lag vector \eqn{\boldsymbol{h}} and an intrinsic stationary and mean stationary random field \eqn{\{Z(\boldsymbol{s}), \boldsymbol{s}\in I \}} it is defined as
#' \deqn{2 \gamma(\boldsymbol{h}) = \text{Var}(Z(\boldsymbol{s}) - Z(\boldsymbol{s}+\boldsymbol{h})),~\boldsymbol{s},\boldsymbol{s}+\boldsymbol{h}\in I.}
#' Three different estimators can be used here to
#' estimate the variogram for a given two-dimensional lag \eqn{\boldsymbol{h}} vector in one of the four directions south-north, east-west, southeast-northwest and southwest-northeast for data on a two-dimensional regular grid.
#' For the different directions the lag vectors has the following forms:
#' \itemize{
#' \item south-north: \eqn{\boldsymbol{h} = \left(\begin{array}{c} 0 \\ x \end{array}\right)}
#' \item east-west: \eqn{\boldsymbol{h} = \left(\begin{array}{c} x \\ 0 \end{array}\right)}
#' \item southwest-northeast: \eqn{\boldsymbol{h} = \left(\begin{array}{c} x \\ x \end{array}\right)}
#' \item southeast-northwest: \eqn{\boldsymbol{h} = \left(\begin{array}{c} x \\ -x \end{array}\right)}.
#' }
#' If \code{estimator = "matheron"} the popular Matheron variogram estimator is used (Matheron, 1962):
#' \deqn{2\widehat{\gamma}(\boldsymbol{h}) = \frac{1}{|N(\boldsymbol{h})|} \sum_{(\boldsymbol{s}_i, \boldsymbol{s}_j)\in N(\boldsymbol{h})} (Z(\boldsymbol{s}_i) - Z(\boldsymbol{s}_j))^2 }
#' with \eqn{N(\boldsymbol{h}) = \{(\boldsymbol{s}_i, \boldsymbol{s}_j)\in I^2: \boldsymbol{s}_j - \boldsymbol{s}_i = \boldsymbol{h}\}} being the set of all pairs in distance \eqn{h}.
#' This estimator is non robust (Genton, 1998, Gierse & Fried 2025).
#'
#' A robust alternative especially for random fields with isolated outliers is the Genton variogram estimator proposed by Genton (1998) (\code{estimator = "genton"}).
#' \eqn{V_i(\boldsymbol{h})} is the \eqn{i}-th element of the set \eqn{V(\boldsymbol{h}) = \{Z(\boldsymbol{s}_l) - Z(\boldsymbol{s}_m): (\boldsymbol{s}_l, \boldsymbol{s}_m)\in N(\boldsymbol{h})\}} containing all differences for locations with distance \eqn{\boldsymbol{h}}.
#' The estimator is given by:
#' \deqn{2\widehat{\gamma}(\boldsymbol{h}) = \left(c \cdot \left(|V_i(\boldsymbol{h}) - V_j(\boldsymbol{h})|: i<j\right)_{(k)}\right)^2}
#' with \eqn{k = \begin{pmatrix}\left[\frac{|N(\boldsymbol{h})|}{2}\right] + 1 \\ 2\end{pmatrix}},  \eqn{[a]} denoting the integer part of \eqn{a}, \eqn{c} being a consistency factor,
#' which is approximately equal to \eqn{2.22} in large Gaussian samples, and \eqn{(\cdot)_{(k)}} denoting the \eqn{k}-th order statistic.
#'
#' Gierse & Fried (2025) proposed two estimators for random fields with an outlier block. Both are multivariate variogram estimators, which estimate the variogram for a set of lags.
#'  The option \code{estimator = "mcd.org"} uses the proposed MCD.org estimator. For a set of lag vectors \eqn{\boldsymbol{h}_1, \ldots, \boldsymbol{h}_{h_\text{max}}} vectors of the following
#'  structure are built
#'  \deqn{\boldsymbol{V} = \begin{pmatrix} Z(\boldsymbol{s}) \\ Z(\boldsymbol{s}+\boldsymbol{h}_1) \\ \vdots \\ Z(\boldsymbol{s}+\boldsymbol{h}_{h_{\max}})
#'  \end{pmatrix}}
#'  The (reweighted) multivariate covariance determinant estimator, a multivariate and highly robust covariance estimator, is used to robust estimate the variance-covariance matrix of the vectors. The function \code{covMcd} of the package \pkg{robustbase}
#'  is used for estimation. The covariance matrix has the following structure
#'  \deqn{	\boldsymbol{\Sigma}_{\boldsymbol{V}} = \begin{pmatrix}
#'  a_0 & a_1 & a_2 & \ldots & a_{h_{\max}} \\
#'  a_1 & a_0 & \ddots & \ddots & \vdots \\
#'  a_2 & \ddots & \ddots & \ddots & a_2    \\
#'  \vdots & \ddots & \ddots & \ddots & a_1 \\
#'  a_{h_{\max}} & \ldots & a_2 & a_1 & a_0
#'  \end{pmatrix}
#'  \in \mathbb{R}^{(h_{\max} + 1) \times (h_{\max}+1)}}
#'  with \eqn{a_0 = \text{Var}(Z(\boldsymbol{s})), a_1 = \text{Cov}(Z(\boldsymbol{s}), Z(\boldsymbol{s} +\boldsymbol{h}_1)), a_2 = \text{Cov}(Z(\boldsymbol{s}), Z(\boldsymbol{s}+\boldsymbol{h}_2))} and so on.
#'  The different estimates for the variance and covariances are averaged. The MCD.org variogram estimator is calculated by
#'  \eqn{2\widehat{\gamma}(\boldsymbol{h}) = 2[\widehat{\text{Var}}(Z(\boldsymbol{s})) - \widehat{\text{Cov}}(Z(\boldsymbol{s}), Z(\boldsymbol{s}+\boldsymbol{h}))]}.
#'
#'  The option \code{estimator = "mcd.diff"} performs the proposed MCD.diff estimator. For this estimator the following vectors of differences of the process are built
#'  \deqn{\boldsymbol{W} = \begin{pmatrix}
#'  Z(\boldsymbol{s}) - Z(\boldsymbol{s}+\boldsymbol{h}_1) \\
#'  Z(\boldsymbol{s}) - Z(\boldsymbol{s}+\boldsymbol{h}_2) \\
#'  \vdots \\
#'  Z(\boldsymbol{s}) - Z(\boldsymbol{s}+\boldsymbol{h}_{h_{\max}})
#'  \end{pmatrix}}
#'  Again the function \code{covMcd} of the package \pkg{robustbase} is used to estimate the variance-covariance matrix of these vectors. The estimated
#'  variogram values are then obtained from the diagonal of the estimated covariance matrix
#'  \deqn{2\widehat{\gamma}(\boldsymbol{h}) = \text{diag}\left(\widehat{\boldsymbol{\Sigma}}\right)}
#'
#' MCD.diff as well as MCD.org uses the MCD estimator for the robust estimation of the covariance matrix of the vectors.
#' These estimator requires at least twice as many vectors as the dimension of the vectors.
#' Therefore \deqn{h_{\max}} must  not be too large relative to the grid size, so that enough vectors can be generated.
#'
#' @return An object of class "varioRob" which is a list with an dataframe for each direction with columns
#' \itemize{
#' \item lag.x: lag in the E-W direction (along the x-axis)
#' \item lag.y: lag in the S-N direction (along the y-axis)
#' \item n: number of differences/vectors the estimation is based on
#' \item variogram: estimated variogram value
#' \item estimator: used estimator}
#'
#' @references
#' Genton, M. (1998). Highly robust variogram estimation. \emph{Mathematical Geology}, 30, 213-221. \doi{https://doi.org/10.1023/A:1021728614555}
#'
#' Gierse, J., & Fried, R. (2025). Nonparametric directional variogram estimation in the presence of outlier blocks. \emph{Statistical Papers}, 66(134). \doi{https://doi.org/10.1007/s00362-025-01754-2}
#'
#' Matheron, G. (1962). Traité de géostatistique appliquée, Tome I. \emph{Mémoires du Bureau de Recherches Géologiques et Miniéres}, no. 14, Editions Technip, 333
#'
#' @seealso  \code{\link{variogram_est_general}}
#'
#' @examples
#' \donttest{
#'
#' ## Simulate an isotropic gaussian random field without outliers
#' dat <- simulate_grf(gridsize = c(20, 20), param.variogram = c(1, 6), n.it = 1)
#'
#' dat <- cbind(dat$data, dat$grid)
#'
#' ## calculate the variogram using MCD.diff estimator for all four directions
#' varog <- variogram_est(data = dat, hmax = c(7,7,5,5),
#'                        direction = c("S-N", "E-W", "SW-NE", "SE-NW"),
#'                        estimator = "mcd.diff")
#' }
#'
#'
#' @importFrom robustbase covMcd
#' @importFrom robustbase Qn
#' @import stats
#' @export variogram_est


variogram_est <- function(data,
                          hmax,
                          direction = c("S-N", "E-W", "SW-NE", "SE-NW"),
                          estimator = "all",
                          cp = NULL,
                          reweighting = TRUE,
                          ...){

  # To avoid error messages caused by incorrect capitalization
  estimator <- tolower(estimator)

  if(!is.data.frame(data)) stop("data needs to be a data.frame.")
  if(ncol(data) != 3) stop("data needs to have three columns.")

  # numbers of directions
  n.direc <- length(direction)

  if(any(!is.element(direction, c("S-N", "E-W", "SW-NE", "SE-NW")))) stop("Only the directions S-N, E-W, SW-NE and SE-NW are allowed.")
  if(!is.logical(reweighting)) stop("reweighting needs to be TRUE or FALSE.")

  if(length(hmax) == 1){ # same hmax for all directions
    hmax <- rep(hmax, n.direc)
  }

  if(is.element("all", estimator)) estimator <- c("matheron", "genton", "mcd.diff", "mcd.org")
  if(any(!is.element(estimator, c("matheron", "genton", "mcd.diff", "mcd.org")))) stop("Invalid estimator. Only the estimators Matheron, Genton, MCD.diff and MCD.org are allowed.")

  if(length(cp) == 1 & !is.null(cp)){ # same correction factor for all direction
    cp <- rep(cp, n.direc)
  }
  ## order the cps in the order S-N, E-W, SW-NE, SE-NW
  cp.ord <- c(cp[direction == "S-N"], cp[direction == "E-W"], cp[direction == "SW-NE"], cp[direction == "SE-NW"])

  grid <- data[,c(2,3)]
  nx <- length(unique(round(grid[,1], 10))) # number of x-coordinates, i. e. number of coordinates in east-west direction
  ny <- length(unique(round(grid[,2], 10))) # number of y-coordinates, i. e. number of coordinates in north-south direction


  # Build a data grid.
  # Save the data as matrix with ny rows and nx columns
  # the element in the i-th row and j-th column is than the value of the process
  # with s = (j, i), i.e. with x-coordinate equal j and y-coordinate equal i
  data.mat <- matrix(NA, nrow = ny, ncol = nx)
  for(i in 1:nrow(grid)){
    data.mat[ny - round(grid[i,2], 10) + 1 , round(grid[i,1], 10)] <- data[i,1]
  } # 1. column in grid is the x-coordinate, i.e. along the x-axis
  # 2. column in grid is the y-coordinate, i.e. along the y-axis


  # Calculate all estimators contained in estimator
  if(is.element("mcd.diff", estimator)) est.mcd.diff <- .MCD.diff(data.mat, hmax, direction, cp.ord, reweighting, ...)
  if(is.element("mcd.org", estimator)) est.mcd.org <- .MCD.org(data.mat, hmax, direction, cp.ord, reweighting, ...)
  if(is.element("matheron", estimator)) est.matheron <- .Matheron(data.mat, hmax, direction)
  if(is.element("genton", estimator)) est.genton <- .Genton(data.mat, hmax, direction)

  # Save the estimators in a List (one element per direction)
  res <- list()
  for(d in 1:n.direc){
    res.list.element <- c()
    for(e in 1:length(estimator)){

      res.list.e <- get(paste0("est.", estimator[e]))
      dic.d  <- names(res.list.e)[d]
      res.list.e <- res.list.e[[d]]

      if(dic.d == "S-N"){
        res.list.e <- cbind("lag.x" = rep(0, hmax[d]), "lag.y" = 1:hmax[d], res.list.e[,-1])
      }
      if(dic.d == "E-W"){
        res.list.e <- cbind("lag.x" = 1:hmax[d], "lag.y" = rep(0, hmax[d]), res.list.e[,-1])
      }
      if(dic.d == "SW-NE"){
        res.list.e <- cbind("lag.x" = 1:hmax[d], "lag.y" = 1:hmax[d], res.list.e[,-1])
      }
      if(dic.d == "SE-NW"){
        res.list.e <- cbind("lag.x" = 1:hmax[d], "lag.y" = -1:-hmax[d], res.list.e[,-1])
      }
      res.list.element <- rbind(res.list.element, res.list.e)
    }
    res[[d]] <- res.list.element
  }
  names(res) <- names( get(paste0("est.", estimator[1])))
  class(res) <- "varioRob"
  return(res)
}


##########
## Function for determining the MCD.diff variogram estimator

.MCD.diff <- function(data.mat, hmax, direction, cp, reweighting, ...){


  ## Build the vectors
  vec.diff <- .build.vectors.MCD.diff(data.mat, hmax, direction)

  ## estimation of the variogram
  # default: use FastMCD Algorithm of Rousseeuw and van Driessen (1999)
  # alternative algorithm are possible through ... (see documentation of covMcd in robustbase)
  est.diff <- lapply(vec.diff, function(l){
    est <- covMcd(t(l), raw.only = !reweighting, use.correction = FALSE, ...)
    return(list("est" = diag(est$cov), "n" = est$n.obs))
  })

  ns <- lapply(est.diff, function(l) l$n)
  est.diff <- lapply(est.diff, function(l) l$est)

  if(!is.null(cp)){
    est.diff <- lapply(1:length(direction), function(l) cp[l] * est.diff[[l]])
  }

  ## Save as dataframe
  res <- lapply(1:length(direction), function(l){
    r <- cbind("lags" = 1:length(est.diff[[l]]), "n" = ns[[l]] ,"variogram" = est.diff[[l]])
    r <- as.data.frame(r)
    r$estimator <- "mcd.diff"
    return(r)
  })
  names(res) <- names(vec.diff)

  return(res)
}

###########
## Function for building the vectors for the MCD.diff variogram estimators

.build.vectors.MCD.diff <- function(data.mat, hmax, direction){

  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in north-south direction

  dic.all <- is.element(c("S-N", "E-W", "SW-NE", "SE-NW"), direction)

  res.diff <- list()

  count.dic <- 0

  ## S-N: build the vectors for the estimation in north-south direction
  if(dic.all[1] == TRUE){

    h <- hmax[direction == "S-N"]

    if(nx*(ny - h) <= h + 1) stop("hmax in S-N direction is to high for this gridsize. An estimation with mcd.diff is not possible for this combination.")

    # build a matrix for storing the vectors: in each column one vector
    vecS.N <- matrix(NA, nrow = h, ncol = nx*(ny - h))
    ind.vec <- 1 # counts the vectors already build


    for(i in 1:nx){ # for each column, i.e. for each x-coordinate,
      # build the following vectors (see also explanations in
      # section 3.2)
      for(j in ny:(h + 1)){
        vecS.N[,ind.vec] <- data.mat[j,i] - data.mat[(j-1):(j-h),i]
        ind.vec <- ind.vec + 1
      }
    }
    count.dic <- count.dic + 1
    res.diff[[count.dic]] <- vecS.N
  }

  ## E-W: build the vectors for the estimation in east-west direction
  if(dic.all[2] == TRUE){

    h <- hmax[direction == "E-W"]

    if(ny*(nx - h) <= h + 1) stop("hmax in E-W direction is to high for this gridsize. An estimation with mcd.diff is not possible for this combination.")

    # build a matrix for storing the vectors: in each column one vector
    vecE.W <- matrix(NA, nrow = h, ncol = ny*(nx - h))
    ind.vec <- 1 # counts the vectors already build

    for(j in 1:ny){ # for each row, i.e. for each y-coordinate,
      # build the following vectors (see also explanations in
      # section 3.2)
      for(i in 1:(nx - h)){
        vecE.W[,ind.vec] <- data.mat[j,i] - data.mat[j,(i+1):(i+h)]
        ind.vec <- ind.vec + 1
      }
    }
    count.dic <- count.dic + 1
    res.diff[[count.dic]] <- vecE.W
  }

  ## SW-NE: build the vectors for the estimation in southwest-northeast direction
  if(dic.all[3] == TRUE){

    h <- hmax[direction == "SW-NE"]

    if((nx - h)*(ny-h) <= h + 1) stop("hmax in SW-NE direction is to high for this gridsize. An estimation with mcd.diff is not possible for this combination.")

    # build a matrix for storing the vectors: in each column one vector
    vecSW.NE <- matrix(NA, nrow = h, ncol = (nx - h)*(ny-h))
    ind.vec <- 1 # counts the vectors already build

    for(i in 1:(nx - h)){ # i: x-coordinate (east-west)
      for(j in ny:(h+1)){ # j: y-coordinate (north-south)
        # (see also explanations in section 3.2)
        vecSW.NE[,ind.vec] <- data.mat[j, i ] - diag(data.mat[(j-1):(j-h), (i+1):(i+h), drop = FALSE])
        ind.vec <- ind.vec + 1
      }
    }
    count.dic <- count.dic + 1
    res.diff[[count.dic]] <- vecSW.NE
  }

  ## SE-NW: build the vectors for the estimation in southeast-northwest direction
  if(dic.all[4] == TRUE){

    h <- hmax[direction == "SE-NW"]

    if((nx - h)*(ny-h) <= h + 1) stop("hmax in SE-NW direction is to high for this gridsize. An estimation with mcd.diff is not possible for this combination.")

    # build a matrix for storing the vectors: in each column one vector
    vecSE.NW <- matrix(NA, nrow = h, ncol = (nx - h)*(ny-h))
    ind.vec <- 1 # counts the vectors already build

    for(i in 1:(nx - h)){    # i: x-coordinate (east-west)
      for(j in (ny -  h):1){ # j: y-coordinate (north-south)
        # (see also explanations in section 3.2)
        vecSE.NW[,ind.vec] <- data.mat[j, i ] - diag(data.mat[(j+1):(j+h), (i+1):(i+h), drop = FALSE])
        ind.vec <- ind.vec + 1
      }
    }
    count.dic <- count.dic + 1
    res.diff[[count.dic]] <- vecSE.NW
  }

  names.dic <- c("S-N", "E-W", "SW-NE", "SE-NW")
  names(res.diff) <- names.dic[dic.all]
  return(res.diff)
}

###########
## Function for determining the MCD.org variogram estimator

.MCD.org <- function(data.mat, hmax, direction, cp, reweighting, ...){

  # Build the vectors
  vec.org <- .build.vectors.MCD.org(data.mat, hmax, direction)

  ## estimation of the variogram
  # default: use FastMCD Algorithm of Rousseeuw and van Driessen (1999)
  # alternative algorithm are possible through ... (see documentation of covMcd in robustbase)
  est.org.cov <- lapply(vec.org, function(l){
    est <- covMcd(t(l), raw.only = !reweighting, use.correction = FALSE, ...)
    return(list("est" = est$cov, "n" = est$n.obs))
  })

  ns <- lapply(est.org.cov, function(l) l$n)
  est.org.cov <- lapply(est.org.cov, function(l) l$est)

  # Calculate the variogram estimation
  est.org <- lapply(est.org.cov, function(l){

    h <- nrow(l) - 1

    estimate <- rep(0, h)
    X <- l
    gamma.0 <- mean(diag(X)) # average the different variance estimations (C(0))
    # contained on the main diagonal
    Xr <- Xl <- X
    for(j in 1:(h - 1)){ # average the different estimations of c(1), ... , C(hmax-1)
      # see explanations in section 3.1
      Xr <- Xr[-1,-(h + 2 -j)]
      Xl <- Xl[-(h + 2 -j),-1]
      cov <- mean(c(diag(Xr), diag(Xl)))
      estimate[j] <- 2 * gamma.0 - 2 * cov # claculate Variogram for lags 1, ..., hmax
    }
    estimate[h] <- (2 * gamma.0 - 2 * mean(Xr[1,2], Xr[2,1])) # calculate variogram for lag hmax

    return(estimate)
  })

  if(!is.null(cp)){
    est.org <- lapply(1:length(direction), function(l) cp[l] * est.org[[l]])
  }

  ## Save as dataframe
  res <- lapply(1:length(direction), function(l){
    r <- cbind("lags" = 1:length(est.org[[l]]), "n" = ns[[l]],"variogram" = est.org[[l]])
    r <- as.data.frame(r)
    r$estimator <- "mcd.org"
    return(r)
  })
  names(res) <- names(vec.org)

  return(res)
}

###########
## Function for building the vectors for the MCD.org variogram estimators

.build.vectors.MCD.org <- function(data.mat, hmax, direction){

  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in north-south direction

  dic.all <- is.element(c("S-N", "E-W", "SW-NE", "SE-NW"), direction)

  res.org <- list()

  count.dic <- 0

  ## S-N: build the vectors for the estimation in north-south direction
  if(dic.all[1] == TRUE){

    h <- hmax[direction == "S-N"]

    if(nx*(ny - h) <= h + 2) stop("hmax in S-N direction is to high for this gridsize. An estimation with mcd.org is not possible for this combination.")

    # build a matrix for storing the vectors: in each column one vector
    vecS.N <- matrix(NA, nrow = h + 1, ncol = nx*(ny - h))
    ind.vec <- 1  # counts the vectors already build

    for(i in 1:nx){ # for each column, i.e. for each x-coordinate,
      # build the following vectors (see also explanations in
      # section 3.1)
      for(j in ny:(h + 1)){
        vecS.N[,ind.vec] <- data.mat[j:(j-h), i]
        ind.vec <- ind.vec + 1
      }
    }
    count.dic <- count.dic + 1
    res.org[[count.dic]] <- vecS.N
  }

  ## E-W: build the vectors for the estimation in east-west direction
  if(dic.all[2] == TRUE){

    h <- hmax[direction == "E-W"]

    if(ny*(nx - h) <= h + 2) stop("hmax in E-W direction is to high for this gridsize. An estimation with mcd.org is not possible for this combination.")

    # build a matrix for storing the vectors: in each column one vector
    vecE.W <- matrix(NA, nrow = h + 1, ncol = ny*(nx - h))
    ind.vec <- 1 # counts the vectors already build

    for(j in 1:ny){  # for each row, i.e. for each y-coordinate,
      # build the following vectors (see also explanations in
      # section 3.1)
      for(i in 1:(nx - h)){
        vecE.W[,ind.vec] <- data.mat[j,i:(i+h)]
        ind.vec <- ind.vec + 1
      }
    }
    count.dic <- count.dic + 1
    res.org[[count.dic]] <- vecE.W
  }

  ## SW-NE: build the vectors for the estimation in southwest-northeast direction
  if(dic.all[3] == TRUE){

    h <- hmax[direction == "SW-NE"]

    if((nx - h)*(ny-h) <= h + 2) stop("hmax in SW-NE direction is to high for this gridsize. An estimation with mcd.org is not possible for this combination.")

    # build a matrix for storing the vectors: in each column one vector
    vecSW.NE <- matrix(NA, nrow = h + 1, ncol = (nx - h)*(ny-h))
    ind.vec <- 1  # counts the vectors already build

    for(i in 1:(nx - h)){ # i: x-coordinate (east-west)
      for(j in ny:(h+1)){ # j: y-coordinate (north-south)
        # (see also explanations in section 3.1)
        vecSW.NE[,ind.vec] <- diag(data.mat[j:(j-h), i:(i+h)])
        ind.vec <- ind.vec + 1
      }
    }
    count.dic <- count.dic + 1
    res.org[[count.dic]] <- vecSW.NE
  }

  ## SE-NW: build the vectors for the estimation in southeast-northwest direction
  if(dic.all[4] == TRUE){

    h <- hmax[direction == "SE-NW"]

    if((nx - h)*(ny-h) <= h + 2) stop("hmax in SE-NW direction is to high for this gridsize. An estimation with mcd.org is not possible for this combination.")

    # build a matrix for storing the vectors: in each column one vector
    vecSE.NW <- matrix(NA, nrow = h + 1, ncol = (nx - h)*(ny-h))
    ind.vec <- 1   # counts the vectors already build

    for(i in 1:(nx - h)){    # i: x-coordinate (east-west)
      for(j in (ny -  h):1){ # j: y-coordinate (north-south)
        # (see also explanations in section 3.1)
        vecSE.NW[,ind.vec] <- diag(data.mat[j:(j+h), i:(i+h)])
        ind.vec <- ind.vec + 1
      }
    }
    count.dic <- count.dic + 1
    res.org[[count.dic]] <- vecSE.NW
  }

  names.dic <- c("S-N", "E-W", "SW-NE", "SE-NW")
  names(res.org) <- names.dic[dic.all]
  return(res.org)
}


###########
## Function for determining the Matheron variogram estimator

.Matheron <- function(data.mat, hmax, direction){

  ## Build differences for all lags
  diffs <- .diff(data.mat, hmax, direction)

  ## estimate the variogram using the Materon estimator
  est.mat <- lapply(diffs, function(d) sapply(1:length(d), function(l) mean(d[[l]]^2, na.rm = TRUE)))

  ## Determine the number of differences
  ns <- lapply(diffs, function(d) sapply(1:length(d), function(l) length(na.omit(d[[l]]))))

  ## Save as dataframe
  res <- lapply(1:length(direction), function(l){
    r <- cbind("lags" = 1:length(est.mat[[l]]), "n" = ns[[l]],"variogram" = est.mat[[l]])
    r <- as.data.frame(r)
    r$estimator <- "matheron"
    return(r)
  })
  names(res) <- names(diffs)

  return(res)
}


###########
## Function for determining the Genton variogram estimator

.Genton <- function(data.mat, hmax, direction){

  ## Build differences for all lags
  diffs <- .diff(data.mat, hmax, direction)

  ## estimate the variogram using the Genton estimator
  ## The Funktion uses the Qn function of the package robustbase
  est.gen <- lapply(diffs, function(d) sapply(1:length(d), function(l) Qn(d[[l]], na.rm = TRUE)^2))

  ## Determine the number of differences
  ns <- lapply(diffs, function(d) sapply(1:length(d), function(l) length(na.omit(d[[l]]))))

  ## Save as dataframe
  res <- lapply(1:length(direction), function(l){
    r <- cbind("lags" = 1:length(est.gen[[l]]), "n" = ns[[l]],"variogram" = est.gen[[l]])
    r <- as.data.frame(r)
    r$estimator <- "genton"
    return(r)
  })
  names(res) <- names(diffs)

  return(res)
}


###########
## Function for determining the differences for all lags and all directions for the Matheron and Genton variogram estimator

.diff <- function(data.mat, hmax, direction){

  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in north-south direction

  # estimation in which direction?
  dic.all <- is.element(c("S-N", "E-W", "SW-NE", "SE-NW"), direction)

  diffs <- list()

  count.dic <- 0

  ## determine lag vectors out of the information in hmax for each direction
  if(is.element("S-N", direction)){

    h <- hmax[direction == "S-N"]

    count.dic <- count.dic + 1

    # determine the differences for all lags of interest
    diffs[[count.dic]] <- lapply(1:h, function(l){
      as.vector(apply(data.mat, 2, function(x) x[1:(ny-l)]-x[(l+1):ny]))
    })
  }

  if(is.element("E-W", direction)){

    h <- hmax[direction == "E-W"]

    count.dic <- count.dic + 1

    # determine the differences for all lags of interest
    diffs[[count.dic]] <- lapply(1:h, function(l){
      as.vector(apply(data.mat, 1, function(x) x[1:(ny-l)]-x[(l+1):ny]))
    })
  }

  if(is.element("SW-NE", direction)){

    h <- hmax[direction == "SW-NE"]

    count.dic <- count.dic + 1

    # determine the differences for all lags of interest
    diffs[[count.dic]] <- lapply(1:h, function(l){
      as.vector(sapply(1:(nx-l), function(x) data.mat[(l+1):ny, x] - data.mat[1:(ny-l), x+1]))
    })
  }

  if(is.element("SE-NW", direction)){

    h <- hmax[direction == "SE-NW"]

    count.dic <- count.dic + 1

    # determine the differences for all lags of interest
    diffs[[count.dic]] <- lapply(1:h, function(l){
      as.vector(sapply(1:(nx-l), function(x) data.mat[1:(ny-l), x] - data.mat[(l+1):ny, x+1]))
    })

    sapply(1:(nx-1), function(x) data.mat[1:(ny-1), x] - data.mat[2:ny, x+1])
  }

  names.dic <- c("S-N", "E-W", "SW-NE", "SE-NW")
  names(diffs) <- names.dic[dic.all]

  return(diffs)

}

