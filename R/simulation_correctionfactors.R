#' Simulate of correction factors for the variogram estimators for small grids
#'
#' Returns simulated finite sample correction factors for (non-)robust directional variogram estimators using a gaussian random field.
#'
#' @param gridsize A two-dimensional vector containing the gridsize in x-direction (nx) and the gridsize in y-direction (ny).
#' @param variogram.est An object of class \code{varioRob} or \code{varioRob}, i.e. an output of the functions \code{\link{variogram_est}} or \code{\link{variogram_est_general}}.
#' @param n.it Numeric. Number of iterations used for the simulation of the correctionfactors; default: 1000.
#' @param estimator A character string or vector of character strings specifying the estimators to be used. Possible estimators are \code{"matheron"}, \code{"genton"}, \code{"mcd.org"} or \code{"mcd.diff"}. If \code{"all"} all four estimators are used. See details; Default is \code{"all"}.
#' @param reweighting  Logical indicating if the reweighted version of the MCD based variogram estimators should be used; default is \code{TRUE}.
#' @param variogram The variogram model. See \code{vgm} of the package \pkg{gstat} for more information; default \code{"Sph"}.
#' @param range Numeric. A start value for the range estimation. If \code{NULL} the start value is half of the largest lag; default is \code{NULL}.
#' @param sill Numeric. A start value for the sill estimation. If \code{NULL} the start value is the mean of the variogram estimation for the three highest lags; default is \code{NULL}.
#' @param nugget Numeric. A start value for the nugget estimation. If \code{NULL} the start value is the 0; default is \code{NULL}.
#' @param seed Numeric. A possible seed for randomnumber generation.
#'
#' @details
#' For small grid sizes and if many lags in one direction should be estimated, the robust estimator needs multiplicative finite sample correction factors (see Gierse & Fried, 2025).
#' The correction factors are simulated based on the grid size as well as the lags of interest using gaussian random fields. The procedure is the following:
#'
#' 1. estimate the range of the variogram given the variogram estimation for the variogram model specified in \code{variogram} using a least squares estimator
#'
#' 2. simulate data (gaussian random field) for the variogram model using the estimated range, a sill of 1 and a nugget of 0
#'
#' 3. estimate the variogram for all simulated data sets
#'
#' 4. calculate the correction factors (see Gierse & Fried, 2025)
#'
#' Only the estimated range is used for simulation of the data, since the sill as well as the nugget effect cannot be estimated from the uncorrected data.
#'
#' For the simulation of the data the function \code{predict} from the package \pkg{gstat} is used.
#'
#' If no start values for the variogram parameters are defined, they are specified as follows
#' \itemize{
#'     \item range: half of the largest estimated lag
#'     \item sill: mean of the last three variogram estimates
#'     \item nugget: 0
#'  }
#'
#' More information about the implementation can be found in documentations of the functions
#' \code{\link{variogram_est}} and \code{\link{variogram_est_general}}.
#'
#' @return A list with two elements
#' \itemize{
#' \item correction factors: a matrix with the simulated correction factors with a column for each estimator and a row for each direction
#' \item ranges: a vector with the estimated range for each direction
#' }
#'
#' @references
#' Gierse, J., & Fried, R. (2025). Nonparametric directional variogram estimation in the presence of outlier blocks. \emph{Statistical Papers}, 66(134). \doi{https://doi.org/10.1007/s00362-025-01754-2}
#'
#'
#' @examples
#' \donttest{
#' ## Simulate an isotropic gaussian random field without outliers
#' dat <- simulate_grf(gridsize = c(20, 20), param.variogram = c(1, 6), n.it = 1)
#'
#' dat <- cbind(dat$data, dat$grid)
#'
#' ## calculate the variogram using different robust estimators for all four directions
#' varog <- variogram_est(data = dat, hmax = c(7,7,5,5),
#'                        direction = c("S-N", "E-W", "SW-NE", "SE-NW"),
#'                        estimator = c("genton", "mcd.diff"))
#'
#' ## simulate correction factors
#' corr <- simulate_correctionfactors(gridsize = c(20, 20),
#'                                    variogram.est = varog,
#'                                    n.it = 1000,
#'                                    estimator = c("genton", "mcd.diff"))
#' }
#'
#' @importFrom robustbase covMcd
#' @importFrom robustbase Qn
#' @importFrom stats optim
#' @importFrom stats median
#' @importFrom stats predict
#' @import gstat
#' @import sp
#' @export simulate_correctionfactors

simulate_correctionfactors <- function(gridsize,
                                       variogram.est,
                                       n.it = 1000,
                                       estimator,
                                       reweighting = TRUE,
                                       variogram = "Sph",
                                       range = NULL,
                                       sill = NULL,
                                       nugget = NULL,
                                       seed = NULL){

  # To avoid error messages caused by incorrect capitalization
  estimator <- tolower(estimator)
  if(is.element("all", estimator)) estimator <- c("matheron", "genton", "mcd.diff", "mcd.org")

  if(!is.vector(gridsize) | length(gridsize) != 2) stop("gridsize needs to be an two-dimensional vector.")
  if(!inherits(variogram.est, "varioRob") & !inherits(variogram.est, "varioRobgen")) stop("variogram.est needs to be of class varioRob or varioRobgen.")
  if(!is.numeric(n.it)) stop("n.it needs to be numeric.")
  if(n.it %% 1 != 0) stop("n.it needs to be an integer number.")
  if(any(!is.element(estimator, c("matheron", "genton", "mcd.diff", "mcd.org")))) stop("Invalid estimator. Only the estimators Matheron, Genton, MCD.diff and MCD.org are allowed.")
  if(!is.logical(reweighting)) stop("reweighting needs to be TRUE or FALSE.")

  if(!is.numeric(range) & !is.null(range)) stop("range needs to be numeric or NULL.")
  if(!is.null(range)){if(range <= 0) stop("range needs to be a positive value.")}

  if(!is.numeric(sill) & !is.null(sill)) stop("sill needs to be numeric or NULL.")
  if(!is.null(sill)){if(sill <= 0) stop("sill needs to be a positive value.")}

  if(!is.numeric(nugget) & !is.null(nugget)) stop("nugget needs to be numeric or NULL.")
  if(!is.null(nugget)){if(nugget <= 0) stop("nugget needs to be a positive value.")}

  if(!is.numeric(seed) & !is.null(seed)) stop("seed needs to be numeric or NULL.")

  # amount of directions
  n.direc <- length(variogram.est)

  # amount of estimators, for which the correction factor should be calculated
  n.est <- length(estimator)

  # Create the grid
  grid <- expand.grid(1:gridsize[1], 1:gridsize[2])
  colnames(grid) <- c("x", "y")
  # needed for gstat
  coordinates(grid) <- ~x+y
  gridded(grid) <- TRUE

  ## simulate correctionfactors based on the variogram estimation for each direciton
  corrections <- c()
  ranges <- c()
  for(d in 1:n.direc){

    est.d <- variogram.est[[d]]

    ## 1. Estimate the Variogramparameter based on the variogram esimation for this direction
    ##    Of interest is the range of the variogram
    range.est <- rep(NA, length(estimator))
    for(e in 1:n.est){ # for each estimator
      est <- est.d[which(est.d$estimator == estimator[e]),]
      dists <- sqrt(est$lag.x^2 + est$lag.y^2)

      # If range, sill, or nugget is NULL, set default values (range and sill: based on the data)
      if(is.null(range)) range <- max(dists)/2  #max(dist)/2
      if(is.null(sill)) sill <- mean(est[(length(est) - 2):length(est),"variogram"])
      if(is.null(nugget)) nugget <- 0

      # estimation of the range
      par.est <- .est.variogram.param(model = variogram, est$variogram, dists, start = c(sill, range, nugget))
      range.est[e] <- par.est[2]
    }


    ## 2. Simulate datasets based on the (median) range estimation
    # variogram model
    model.sim <- vgm(psill = 1, model = variogram, range = median(range.est), nugget = 0) # ggfs. nugget noch anpassen

    # gstat object (no data set → unconditional simulation)
    g <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE, beta = 0, model = model.sim)

    # simulation
    if(!is.null(seed)) set.seed(seed)
    suppressMessages(sim.data <- predict(g, grid, nsim = n.it)@data)


    ## 3. Estimate the variogram with all estimators
    # transform data in the needed format
    sim.data <- lapply(1:n.it, function(it){
      data.frame(sim.data[,it], grid)
    })
    if(inherits(variogram.est, "varioRob")){
      est.vario <- lapply(sim.data, function(x){
        variogram_est(x, hmax = nrow(est), direction = names(variogram.est)[d], estimator = estimator, reweighting = reweighting)
      })
    }
    if(inherits(variogram.est, "varioRobgen")){
      est.vario <- lapply(sim.data, function(x){
        variogram_est_general(x, h = as.matrix(est[,1:2]), estimator = estimator, reweighting = reweighting)
      })
    }


    ## 4. Calculate the correction factor
    # calculate true value
    true <- 2*variogramLine(model.sim, dist_vector = dists)$gamma
    corr <- lapply(1:n.est, function(e){
      norm.est <-lapply(1:n.it, function(it){
        est.e <- est.vario[[it]][[1]][which(est.vario[[it]][[1]]$estimator == estimator[e]),]
        est.norm <- est.e$variogram/true
        return(est.norm)
      })
      norm.est <- do.call(rbind, norm.est)
      corr.est <- mean(rowMeans(norm.est)) # average over the different lags and than over the n.it
      return(1/corr.est)
    })
    corr <- unlist(corr)
    names(corr) <- estimator

    corrections <- rbind(corrections, corr)
    ranges <- c(ranges, mean(range.est))
  }
  rownames(corrections) <- names(variogram.est)
  names(ranges) <- names(variogram.est)

  res <- list("correctionfactors" = corrections, "ranges" = ranges)

  return(res)
}


##########
##  Function for KQ estimation of the of the variogram parameter

.est.variogram.param <- function(model, estimation, dist, start){

  # minimizing function
  min_func <- function(param){
    variogram.model <- vgm(psill = param[1], model = model, range = param[2], nugget = param[3])
    pred <- variogramLine(variogram.model, dist_vector = dist)$gamma
    return(sum((estimation - 2*pred)^2))
  }

  # optimization
  fit <- optim(start, min_func)$par

  return(fit)
}

