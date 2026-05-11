#' Nonparametric directional variogram estimation
#'
#' Returns robust and non-robust directional variogram estimations in arbitrary directions for data on a two-dimensional regular grid.
#'
#' @param data A dataframe with three columns containing the data (first column) and the grid (second and third column).
#' @param h A matrix with two columns containing the lag vectors for which the variogram should be estimated. Each row contain one lag vector.
#' @param estimator A character string or vector of character strings specifying the estimators to be used. Possible estimators are \code{"matheron"}, \code{"genton"}, \code{"mcd.org"} or \code{"mcd.diff"}. If \code{"all"} all four estimators are used. See details; Default is "all".
#' @param cp A vector containing possible finite sample correction factors for each direction for the MCD based variogram estimators. If it has length 1 the same correction factor is used for all directions.
#' @param reweighting Logical indicating if the reweighted version of the MCD based variogram estimators should be used; default is \code{TRUE}.
#' @param ... Optional argument for the function \code{covMcd} of the package \pkg{robustbase} used for the MCD based variogram estimators
#'
#' @details
#' The variogram characterizes the spatial dependence of a random field. For a given lag vector \eqn{\boldsymbol{h}} and an intrinsically and mean stationary random field \eqn{\{Z(\boldsymbol{s}), \boldsymbol{s}\in I \}} it is defined as
#' \deqn{2 \gamma(\boldsymbol{h}) = \text{Var}(Z(\boldsymbol{s}) - Z(\boldsymbol{s}+\boldsymbol{h})),~\boldsymbol{s},\boldsymbol{s}+\boldsymbol{h}\in I.}
#' Three different estimators are implemented here to
#' estimate the variogram for arbitrary given two-dimensional lag vectors \eqn{h} for data on a two-dimensional regular grid.
#'
#' If \code{estimator = "matheron"} the popular Matheron variogram estimator is used (Matheron, 1962):
#' \deqn{2\widehat{\gamma}(\boldsymbol{h}) = \frac{1}{|N(\boldsymbol{h})|} \sum_{(\boldsymbol{s}_i, \boldsymbol{s}_j)\in N(\boldsymbol{h})} (Z(\boldsymbol{s}_i) - Z(\boldsymbol{s}_j))^2 }
#' with \eqn{N(\boldsymbol{h}) = \{(\boldsymbol{s}_i, \boldsymbol{s}_j)\in I^2: \boldsymbol{s}_j - \boldsymbol{s}_i = \boldsymbol{h}\}} being the set of all pairs in distance \eqn{h}.
#' This estimator is non robust (Genton, 1998, Gierse & Fried 2025).
#'
#' A robust alternative, particularly suitable for random fields containing isolated outliers, is the Genton variogram estimator proposed by Genton (1998) (\code{estimator = "genton"}).
#' Let \eqn{V_i(\boldsymbol{h})} denote the \eqn{i}-th element of the set \eqn{V(\boldsymbol{h}) = \{Z(\boldsymbol{s}_l) - Z(\boldsymbol{s}_m): (\boldsymbol{s}_l, \boldsymbol{s}_m)\in N(\boldsymbol{h})\}}
#' which contains all pairwise differences of observations separated by the lag vector \eqn{\boldsymbol{h}}.
#' The estimator is defined as:
#' \deqn{2\widehat{\gamma}(\boldsymbol{h}) = \left(c \cdot \left(|V_i(\boldsymbol{h}) - V_j(\boldsymbol{h})|: i<j\right)_{(k)}\right)^2}
#' where \eqn{k = \begin{pmatrix}\left[\frac{|N(\boldsymbol{h})|}{2}\right] + 1 \\ 2\end{pmatrix}},  \eqn{[a]} denoting the integer part of \eqn{a}, \eqn{c} is a consistency factor
#' (approximately equal to \eqn{2.22} in large Gaussian samples) and \eqn{(\cdot)_{(k)}} denoting the \eqn{k}-th order statistic.
#'
#' Gierse & Fried (2025) proposed two estimators designed for random fields affected by blocks of outliers. Both approaches are multivariate variogram estimators that simultaneously estimate the variogram over a set of lags.
#'  The option \code{estimator = "mcd.org"} uses the proposed MCD.org estimator. For a collection of lag vectors \eqn{\boldsymbol{h}_1, \ldots, \boldsymbol{h}_{h_\text{max}}}, vectors of the following
#'  \deqn{\boldsymbol{V} = \begin{pmatrix} Z(\boldsymbol{s}) \\ Z(\boldsymbol{s}+\boldsymbol{h}_1) \\ \vdots \\ Z(\boldsymbol{s}+\boldsymbol{h}_{h_{\max}})
#'  \end{pmatrix}}
#'  are constructed.
#'  To robustly estimate the associated variance–covariance matrix, the (reweighted) minimum covariance determinant (MCD) estimator is employed.
#'  This estimator is a highly robust multivariate method for covariance estimation. The function \code{covMcd} from the robustbase package is used for this.
#'  The resulting covariance matrix has the structure
#'  \deqn{	\boldsymbol{\Sigma}_{\boldsymbol{V}} = \begin{pmatrix}
#'  a_0 & a_1 & a_2 & \ldots & a_{h_{\max}} \\
#'  a_1 & a_0 & \ddots & \ddots & \vdots \\
#'  a_2 & \ddots & \ddots & \ddots & a_2    \\
#'  \vdots & \ddots & \ddots & \ddots & a_1 \\
#'  a_{h_{\max}} & \ldots & a_2 & a_1 & a_0
#'  \end{pmatrix}
#'  \in \mathbb{R}^{(h_{\max} + 1) \times (h_{\max}+1)}}
#'  where \eqn{a_0 = \text{Var}(Z(\boldsymbol{s})), a_1 = \text{Cov}(Z(\boldsymbol{s}), Z(\boldsymbol{s} +\boldsymbol{h}_1)), a_2 = \text{Cov}(Z(\boldsymbol{s}), Z(\boldsymbol{s}+\boldsymbol{h}_2))} and so forth.
#'  The variance and covariance estimates obtained from this matrix are averaged accordingly. The MCD.org variogram estimator is then given by
#'  \eqn{2\widehat{\gamma}(\boldsymbol{h}) = 2[\widehat{\text{Var}}(Z(\boldsymbol{s})) - \widehat{\text{Cov}}(Z(\boldsymbol{s}), Z(\boldsymbol{s}+\boldsymbol{h}))]}.
#'
#'  The option \code{estimator = "mcd.diff"} implements the proposed MCD.diff estimator. For this approach, vectors of pairwise differences of the process are constructed as
#'  \deqn{\boldsymbol{W} = \begin{pmatrix}
#'  Z(\boldsymbol{s}) - Z(\boldsymbol{s}+\boldsymbol{h}_1) \\
#'  Z(\boldsymbol{s}) - Z(\boldsymbol{s}+\boldsymbol{h}_2) \\
#'  \vdots \\
#'  Z(\boldsymbol{s}) - Z(\boldsymbol{s}+\boldsymbol{h}_{h_{\max}})
#'  \end{pmatrix}}
#'  As before, the function \code{covMcd} from the package \pkg{robustbase} is used to obtain a robust estimate of the variance–covariance matrix of these vectors.
#'  The variogram estimates are then derived from the diagonal elements of the estimated covariance matrix, i.e.,
#'  \deqn{2\widehat{\gamma}(\boldsymbol{h}) = \text{diag}\left(\widehat{\boldsymbol{\Sigma}}\right)}
#'
#' MCD.diff as well as MCD.org uses the MCD estimator for the robust estimation of the covariance matrix of the vectors.
#' These estimator requires at least twice as many vectors as the dimension of the vectors.
#' Therefore \eqn{h_{\max}} must not be too large relative to the grid size, so that enough vectors can be generated.
#'
#' Both the MCD.diff and MCD.org estimators rely on the MCD approach for the robust estimation of the covariance matrix of the constructed vectors.
#' These estimators require that the number of available vectors is at least twice the dimension of the vectors. Consequently, the choice of
#' \eqn{h_{\max}} should not be too large relative to the grid size, ensuring that a sufficient number of vectors can be formed.
#'
#' If the set of lag vectors of interest contains only vectors  in the following directions south-north, east-west, southeast-northwest and southwest-northeast,
#' then the function \code{variogram_est} should be used, as it has better runtime performance.
#'
#' @return An object of class "varioRobgen" which is a list with an dataframe for each direction with columns
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
#' @seealso \code{\link{variogram_est}}
#'
#' @examples
#' \donttest{
#'
#' ## Simulate an anisotropic gaussian random field without outliers
#' dat <- simulate_grf(gridsize = c(20, 20), param.variogram = c(1, 6),
#'                     aniso.param = c(pi/4, 2), n.it = 1)
#'
#' dat <- cbind(dat$data, dat$grid)
#'
#' ## calculate the variogram using all implemented estimator for different lags
#' varog <- variogram_est_general(data = dat,
#'                                h = rbind(c(2,1), c(4,2), c(1,2), c(2,4)),
#'                                estimator = "all")
#' }
#'
#' @importFrom robustbase covMcd
#' @importFrom robustbase Qn
#' @import stats
#' @export variogram_est_general

variogram_est_general <- function(data,
                                  h,
                                  estimator = "all",
                                  cp = NULL,
                                  reweighting = TRUE,
                                  ...){

  # To avoid error messages caused by incorrect capitalization
  estimator <- tolower(estimator)

  if(!is.data.frame(data)) stop("data needs to be a data.frame.")
  if(ncol(data) != 3) stop("data needs to have three columns.")

  # numbers of lag vectors
  n.lag <- nrow(h)

  if(is.element("all", estimator)) estimator <- c("matheron", "genton", "mcd.diff", "mcd.org")
  if(any(!is.element(estimator, c("matheron", "genton", "mcd.diff", "mcd.org")))) stop("Invalid estimator. Only the estimators Matheron, Genton, MCD.diff and MCD.org are allowed.")

  if(!is.logical(reweighting)) stop("reweighting needs to be TRUE or FALSE.")

  if(length(cp) == 1 & !is.null(cp)){ # same correction factor for each lag vector in h
    cp <- rep(cp, n.lag)
  }

  if(!is.matrix(h)) stop("h needs to be a matrix.")
  if(ncol(h) != 2) stop("h needs to have two colums. One for the east-west direction and one for the north-south direction.")

  # standadize the coordinates
  grid <- data[,c(2,3)]
  nx <- length(unique(round(grid[,1], 10))) # number of x-coordinates, i. e. number of coordinates in east-west direction
  ny <- length(unique(round(grid[,2], 10))) # number of y-coordinates, i. e. number of coordinates in north-south direction


  # Build a data grid.
  # Save the data as matrix with ny rows and nx columns
  # the element in the i-th row and j-th column is than the value of the process
  # with s = (j, i), i.e. with x-coordinate equal j and y-coordinate equal i
  data.mat <- matrix(NA, nrow = ny, ncol = nx)
  for(i in 1:nrow(grid)){
    data.mat[ny - round(grid[i,2], 10) + 1, round(grid[i,1], 10)] <- data[i,1]
  } # 1. column in grid is the x-coordinate, i.e. along the x-axis
  # 2. column in grid is the y-coordinate, i.e. along the y-axis


  # Calculate all estimators contained in estimator
  if(is.element("mcd.diff", estimator)) est.mcd.diff <- .MCD.diff.general(data.mat, h, cp, reweighting, ...)
  if(is.element("mcd.org", estimator)) est.mcd.org <- .MCD.org.general(data.mat, h, cp, reweighting, ...)
  if(is.element("matheron", estimator)) est.matheron <- .Matheron.general(data.mat, h)
  if(is.element("genton", estimator)) est.genton <- .Genton.general(data.mat, h)

  # Save the estimators in a List (one element per direction)
  normalize.h <- t(apply(h, 1, .normalize.row)) # normalize the lag vector to determine the direction of the lag vector
  h <- as.data.frame(h)
  h$factor <- paste0(normalize.h[,1],",", normalize.h[,2])
  split.h <- split(h, h$factor) # a list for each lag vector

  res.long <- c()
  for(e in 1:length(estimator)){
    res.long <- rbind(res.long, get(paste0("est.", estimator[e])))
  }
  res <- list()
  for(l in 1:length(split.h)){
    h.l <- split.h[[l]]
    h.l$lag <- paste(h.l[,1], h.l[,2])
    res.l <- lapply(1:nrow(res.long), function(x){
      if(is.element(res.long$lags[x], h.l$lag)) res.long[x,]
    })
    res.l <- do.call(rbind, res.l)
    h.new <- lapply(1:nrow(res.l), function(x){
      lag <- strsplit(res.l$lags[x], " ")
      return(as.numeric(c(lag[[1]][1], lag[[1]][2])))
    })
    h.new <- do.call(rbind, h.new)
    res.new <- data.frame("lag.x" = h.new[,1], "lag.y" = h.new[,2], "n" = res.l$n, "variogram" = res.l$variogram, "estimator" = res.l$estimator)
    res[[l]] <- res.new
  }

  names(res) <- names(split.h)
  class(res) <- "varioRobgen"
  return(res)
}


##########
## Function for determining the MCD.diff variogram estimator

.MCD.diff.general <- function(data.mat, h, cp, reweighting, ...){


  ## Build the vectors
  vec.diff <- .build.vectors.MCD.diff.general(data.mat, h)

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
    est.diff <- lapply(1:length(est.diff), function(l) cp[l] * est.diff[[l]])
  }

  res <- lapply(1:length(est.diff), function(l){
    data.frame("variogram" = est.diff[[l]], "n" = ns[[l]], "lags" = names(est.diff[[l]]), "estimator" = "mcd.diff")
  })
  res <- do.call(rbind, res)
  rownames(res) <- NULL

  ## Multiply with correction factor
  if(!is.null(cp)){

    # order of h
    h.ord <- paste(h[,1], h[,2])

    # multiply the correction factor
    res.cp <- lapply(1:nrow(res), function(x){
      res[x, "variogram"] <- as.numeric(res[x, "variogram"]) * cp[which(res[x, "lags"] == h.ord)]
      return(res[x,])
    })
    res <- do.call(rbind, res.cp)
  }

  return(res)
}

###########
## Function for building the vectors for the MCD.diff variogram estimators

.build.vectors.MCD.diff.general <- function(data.mat, h){

  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in north-south direction
  grid <- expand.grid(1:nx, 1:ny)

  # which lags are in the same direction
  normalize.h <- t(apply(h, 1, .normalize.row)) # normalize the lag vector to determine lag vectors in the same direction
  h <- as.data.frame(h)
  h$factor <- paste0(normalize.h[,1],",", normalize.h[,2])
  split.h <- split(h, h$factor) # for each direction one list

  # for each direction build the vectors needed for MCD.diff
  res.diff <- lapply(split.h, function(lag){
    grid.lag <- cbind(grid[,1] + lag[nrow(lag), 1], grid[,2] + lag[nrow(lag), 2])
    grid.true <- grid.lag[,1] > 0 & grid.lag[,1] <= nx & grid.lag[,2] > 0 & grid.lag[,2] <= ny
    grid.lag <- grid.lag[grid.true,]

    grid.org <- grid[grid.true, ]

    diff <- apply(grid.org, 1, function(l){
      l <- unlist(l)
      d <- sapply(1:nrow(lag), function(x){
        # notation in R: rownames column names are exactly the opposite
        data.mat[nx - l[2] + 1, l[1]] - data.mat[nx - (l[2] + lag[x,2]) + 1, l[1] + lag[x,1]]
      })
      d <- as.matrix(d, nrow = 1)
      return(d)
    })
    if(!is.matrix(diff)) diff <- matrix(diff, nrow = 1)
    rownames(diff) <- paste(lag[,1], lag[,2])

    if(ncol(diff) <= nrow(diff) + 1) stop("Too many lags for estimation for this gridsize. An estimation with mcd.diff is not possible in all directions for this combination.")

    return(diff)
  })

  return(res.diff)
}

###########
## Function for determining the MCD.org variogram estimator

.MCD.org.general <- function(data.mat, h, cp, reweighting, ...){

  # Build the vectors
  vec.org <- .build.vectors.MCD.org.general(data.mat, h)

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
    if(h > 1){
      for(j in 1:(h - 1)){ # average the different estimations of c(1), ... , C(hmax-1)
        # see explanations in section 3.1
        Xr <- Xr[-1,-(h + 2 -j)]
        Xl <- Xl[-(h + 2 -j),-1]
        cov <- mean(c(diag(Xr), diag(Xl)))
        estimate[j] <- 2 * gamma.0 - 2 * cov # claculate Variogram for lags 1, ..., hmax
      }
    }
    estimate[h] <- (2 * gamma.0 - 2 * mean(Xr[1,2], Xr[2,1])) # calculate variogram for lag hmax
    names(estimate) <- rownames(l)[-1]
    return(estimate)
  })


  res <- lapply(1:length(est.org), function(l){
    data.frame("variogram" = est.org[[l]], "n" = ns[[l]], "lags" = names(est.org[[l]]), "estimator" = "mcd.org")
  })
  res <- do.call(rbind, res)
  rownames(res) <- NULL

  ## Multiply with correction factor
  if(!is.null(cp)){

    # order of h
    h.ord <- paste(h[,1], h[,2])

    # multiply the correction factor
    res.cp <- lapply(1:nrow(res), function(x){
      res[x, "variogram"] <- as.numeric(res[x, "variogram"]) * cp[which(res[x, "lags"] == h.ord)]
      return(res[x,])
    })
    res <- do.call(rbind, res.cp)
  }

  return(res)
}

###########
## Function for building the vectors for the MCD.org variogram estimators

.build.vectors.MCD.org.general <- function(data.mat, h, direction){

  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in north-south direction
  grid <- expand.grid(1:nx, 1:ny)

  # which lags are in the same direction
  normalize.h <- t(apply(h, 1, .normalize.row)) # normalize the lag vector to determine lag vectors in the same direction
  h <- as.data.frame(h)
  h$factor <- paste0(normalize.h[,1],",", normalize.h[,2])
  split.h <- split(h, h$factor)

  # for each direction build the vectors needed for MCD.org
  res.org <- lapply(split.h, function(lag){
    grid.lag <- cbind(grid[,1] + lag[nrow(lag), 1], grid[,2] + lag[nrow(lag), 2])
    grid.true <- grid.lag[,1] > 0 & grid.lag[,1] <= nx & grid.lag[,2] > 0 & grid.lag[,2] <= ny
    grid.lag <- grid.lag[grid.true,]

    grid.org <- grid[grid.true, ]

    org <- apply(grid.org, 1, function(l){
      l <- unlist(l)
      # notation in R: rownames column names are exactly the opposite
      d <- c(data.mat[nx - l[2] + 1, l[1]], sapply(1:nrow(lag), function(x){
        data.mat[nx - (l[2] + lag[x,2]) + 1, l[1] + lag[x,1]]
      }))
      d <- as.matrix(d, nrow = 1)
      return(d)
    })
    if(!is.matrix(org)) org <- matrix(org, nrow = 1)
    rownames(org) <- c("", paste(lag[,1], lag[,2]))

    if(ncol(org) <= nrow(org) + 1) stop("Too many lags for estimation for this gridsize. An estimation with mcd.org is not possible in all directions for this combination.")

    return(org)
  })

  return(res.org)
}


###########
## Function for determining the Matheron variogram estimator

.Matheron.general <- function(data.mat, h){

  ## Build differences for all lags
  diffs <- .diff.general(data.mat, h)

  ## estimate the variogram using the Materon estimator
  est.mat <- lapply(diffs, function(d) mean(d^2, na.rm = TRUE))

  ## Determine the number of differences
  ns <- lapply(diffs, function(d) length(na.omit(d)))

  ## Save as dataframe
  res <- lapply(1:nrow(h), function(l){
    data.frame("variogram" = est.mat[[l]], "n" = ns[[l]], "lags" = paste(h[l,1], h[l,2]) , "estimator" = "matheron")
  })
  res <- do.call(rbind, res)

  return(res)
}


###########
## Function for determining the Genton variogram estimator

.Genton.general <- function(data.mat, h){

  ## Build differences for all lags
  diffs <- .diff.general(data.mat, h)

  ## estimate the variogram using the Genton estimator
  ## The Funktion uses the Qn function of the package robustbase
  est.gen <- lapply(diffs, function(d) Qn(d, na.rm = TRUE)^2)

  ## Determine the number of differences
  ns <- lapply(diffs, function(d) length(na.omit(d)))

  ## Save as dataframe
  res <- lapply(1:nrow(h), function(l){
    data.frame("variogram" = est.gen[[l]], "n" = ns[[l]], "lags" = paste(h[l,1], h[l,2]) , "estimator" = "genton")
  })
  res <- do.call(rbind, res)

  return(res)
}



###########
## Function for determining the differences for all lags and all directions for the Matheron and Genton variogram estimator

.diff.general <- function(data.mat, h){

  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in north-south direction

  grid <- expand.grid(1:nx, 1:ny)

  # calculate the differences for each lag vector of interest
  diffs <- lapply(1:nrow(h), function(x){
    lag <- h[x,]
    grid.lag <- cbind(grid[,1] + lag[1], grid[,2] + lag[2])
    grid.true <- grid.lag[,1] > 0 & grid.lag[,1] <= nx & grid.lag[,2] > 0 & grid.lag[,2] <= ny
    grid.lag <- grid.lag[grid.true,]

    grid.org <- grid[grid.true, ]

    diff <- sapply(1:nrow(grid.lag), function(l){
      data.mat[nx - grid.org[l,2] + 1, grid.org[l,1]] - data.mat[nx - grid.lag[l,2] + 1, grid.lag[l,1]]
    })
  })

  return(diffs)

}

##########
# Function for normalization of the lag vectors

.normalize.row <- function(x, tol = 1e-12){
  x/max(abs(x[x!=0]))
}
