#' Simulation of gaussian random fields on a regular grid with different types of outliers
#'
#' Returns simulated (an)isotropy gaussian random fields on a regular grid without outliers, with isolated outliers, or with an outlier block.
#'
#' @param gridsize Integer vector of length 2 containing the gridsize in x-direction and the gridsize in y-direction
#' @param dist.outlier A function to draw random numbers from the distribution of the outliers; default: \code{rnorm}
#' @param variogram The variogram model. See function \code{cov.spatial} of the package \pkg{geoR} for more information; default: \code{spherical}.
#' @param param.variogram A two-dimensional vector containing the parameters of the variogram. The first element specifies the (partial) sill and the second one the range of the variogram.
#' @param nugget Numeric. Specifies the nugget effect; default: \code{NULL}
#' @param aniso.param A two-dimensional vector containing the parameters of the geometric anisotropy. The first element is the rotation angle and the second one the anisotropy ratio. If \code{NULL} an isotropic random field is simulated; default: \code{NULL}.
#' @param out.type A character string specifying the type of outlier (\code{"block"} or \code{"isolated"}. If \code{NULL} a random field without outliers is simulated. See also details; default: \code{NULL}.
#' @param block.type A character string only needed for generation of block outliers. It specifies wheather the outlier block should be nearly quadratic (\code{"square"}), as rectangular as possible (\code{"rectangle"}) or if it is random (\code{"random"}).
#' @param amount  Numeric. A value between 0 and 1 specifying the amount of outlier; default: \code{NULL}.
#' @param param.outlier  A vector specifying the parameters of the outlier distribution; default: \code{NULL}.
#' @param mixture logical. Specifying whether the outlier should be generated as mixture model (\code{TRUE}) or as additive outliers (\code{FALSE}); default \code{TRUE}.
#' @param n.it Numeric. Number of iterations; default: \code{1000}.
#' @param seed Numeric. A possible seed for randomnumber generation.
#'
#' @details Gaussian random fields (GRF) without outliers, with isolated outliers or with block outliers can be simulated. The GRF without outliers is simulated
#'          using the function \code{grf} of package \pkg{geoR}.
#'          GRFs with outliers are generated in the following manner:
#'
#'          Isolated outliers (occur randomly of the grid):
#'          \eqn{Y(\boldsymbol{s})|N_{n_0} = \begin{cases} Z(\boldsymbol{s}),~\text{ if } \boldsymbol{s} \notin N_{n_0} \\ W,~ \text{ if } \boldsymbol{s} \in N_{n_0}\end{cases}}
#'
#'          Block outliers (occur spatially aggregated in a block):
#'          \eqn{Y(\boldsymbol{s})|(U=\boldsymbol{s}_0) = \begin{cases} Z(\boldsymbol{s}),~\text{ if } \boldsymbol{s} \notin N(\boldsymbol{s}_0) \\W,~ \text{ if } \boldsymbol{s} \in N(\boldsymbol{s}_0)\end{cases}}
#'
#'          Thereby \eqn{Z(\boldsymbol{s})} is the GRF without outliers, \eqn{N_{n_0}} is a set of \eqn{n_0} randomly selected locations of the grid, \eqn{W} is the distribution of the outliers, \eqn{N(\boldsymbol{s}_0)} is a block neighbourhood around the gridpoint \eqn{\boldsymbol{s}_0} and \eqn{U} is a random variable which select one gridpoint.
#'          For more information see Gierse and Fried (2025).
#'
#' @return A list with two elements
#' \itemize{
#' \item data: a matrix with n.it columns; in each column is one simulated data set
#' \item grid: a data.frame with two columns; the first column contains the x-coordinate and the second column the y-coordinate
#' }
#'
#' @references
#' Gierse, J., & Fried, R. (2025). Nonparametric directional variogram estimation in the presence of outlier blocks. \emph{Statistical Papers}, 66(134). \doi{https://doi.org/10.1007/s00362-025-01754-2}
#'
#'
#' @examples
#' \donttest{
#'  ## Simulation of an anisotrop GRF without outliers
#'  dat1 <-  simulate_grf(gridsize = c(20, 20), param.variogram = c(1, 6),
#'                        aniso.param = c(pi/4, 2), n.it = 1)
#'
#'  ## Simulation of an anisotrop GRF with isolated outliers
#'  dat2 <-  simulate_grf(gridsize = c(20, 20), param.variogram = c(1, 6),
#'                        aniso.param = c(pi/4, 2), out.type = "isolated",
#'                        amount = 0.1, param.outlier = c(0,5), n.it = 1)
#'
#'  ## Simulation of an anisotrop GRF with an quadratic outlier block
#'  dat2 <-  simulate_grf(gridsize = c(20, 20), param.variogram = c(1, 6),
#'                        aniso.param = c(pi/4, 2), out.type = "block",
#'                        block.type = "square", amount = 0.1,
#'                        param.outlier = c(0,5), n.it = 1)
#'
#' }
#'
#' @import geoR
#' @importFrom stats rnorm
#' @export simulate_grf

simulate_grf <- function(gridsize,
                         dist.outlier = rnorm,
                         variogram = "spherical",
                         param.variogram,
                         nugget = NULL,
                         aniso.param = NULL,
                         out.type = NULL,
                         block.type = NULL,
                         amount = NULL,
                         param.outlier = NULL,
                         mixture = TRUE,
                         n.it = 1000,
                         seed = NULL){

  if(!is.vector(gridsize) | length(gridsize) != 2) stop("gridsize needs to be an two-dimensional vector.")
  if(!is.vector(param.variogram) | length(param.variogram) != 2) stop("param.variogram needs to be an two-dimensional vector.")
  if(!is.null(aniso.param)){if(!is.vector(aniso.param) | length(aniso.param) != 2) stop("aniso.param needs to be an two-dimensional vector.")}

  if(!is.function(dist.outlier)) stop("dist.outlier needs to be a function.")

  if(!is.numeric(nugget) & !is.null(nugget)) stop("nugget needs to be numeric or NULL.")
  if(!is.numeric(n.it)) stop("n.it needs to be numeric.")
  if(n.it %% 1 != 0) stop("n.it needs to be an integer number.")
  if(!is.numeric(amount) & !is.null(amount)) stop("amount needs to be numeric or NULL.")
  if(!is.null(amount)){if(amount < 0 | amount > 1) stop("amount needs to be number between 0 and 1")}

  if(!is.null(out.type)){
    out.type <- tolower(out.type)
    if(!is.element(out.type, c("block", "isolated"))) stop("Invalid outlier type. Only the outlier types block or isolated are allowed.")
    if(out.type == "isolated" & !is.null(block.type)) warning("block.type is only relevant if out.type = \"block\" and is therefore here ignored.")
  }

  if(is.null(out.type) & !is.null(block.type)) warning("block.type is only relevant if out.type = \"block\" and is therefore here ignored.")

  if(!is.null(block.type)){
    block.type <- tolower(block.type)
    if(out.type == "block" & !is.element(block.type, c("random", "square", "rectangle"))) stop("Invalid block type. Only the block types random, square or rectangle are allowed.")
  }

  if(!is.logical(mixture)) stop("mixture needs to be logical.")

  variogram <- tolower(variogram)


  # simulate data without outliers
  set.seed(seed)
  if(is.null(aniso.param)){ # isotropy
    if(is.null(nugget)){ # without nugget effect
      data <- grf(grid = expand.grid(1:gridsize[1], 1:gridsize[2]),
                  cov.model = variogram, cov.pars = param.variogram,
                  messages = FALSE, nsim = n.it)
    } else{ # nugget effect
      data <- grf(grid = expand.grid(1:gridsize[1], 1:gridsize[2]),
                  cov.model = variogram, cov.pars = param.variogram, nugget = nugget,
                  messages = FALSE, nsim = n.it)
    }
  } else{ # anisotropy
    if(is.null(nugget)){ # without nugget effect
      data <- grf(grid = expand.grid(1:gridsize[1], 1:gridsize[2]),
                  cov.model = variogram, cov.pars = param.variogram, aniso.pars = aniso.param,
                  messages = FALSE, nsim = n.it)
    } else{ # nugget effect
      data <- grf(grid = expand.grid(1:gridsize[1], 1:gridsize[2]),
                  cov.model = variogram, cov.pars = param.variogram, nugget = nugget, aniso.pars = aniso.param,
                  messages = FALSE, nsim = n.it)
    }
  }

  grid.coord <- expand.grid(1:gridsize[1], 1:gridsize[2]) # save the grid of the data
  data <- data$data # save the simulated data (for each iteration one column)
  if(n.it == 1){
    data <- matrix(data, ncol = 1)
  }

  # if desired include outliers
  if(!is.null(out.type)){
    for(n in 1:n.it){ # for each iteration

      # generate the coordinates of the outliers
      coords.outlier <- .gen.coords(type = out.type, amount = amount, gridsize = gridsize, block.type = block.type)
      row.out <- apply(coords.outlier, 1, function(i) which(as.numeric(i[1]) == grid.coord[,1] & as.numeric(i[2]) == grid.coord[,2]))
      names(row.out) <- NULL

      # simulate the values of the outliers
      if(length(param.outlier) == 2){
        cont <- dist.outlier(nrow(coords.outlier), param.outlier[1], param.outlier[2])
      } else{
        cont <- dist.outlier(nrow(coords.outlier), param.outlier[1])
      }

      # add the outliers to the simulated data
      if(mixture){ # replaced the original data by the outliers (mixture model)
        data[row.out, n] <- cont
      } else{ # add the outliers to the original data (additive model)
        data[row.out, n] <- data[row.out, n] + cont
      }

    }
  }

  return(res <- list(data = data, grid = grid.coord))

}

##########
## Function to determine the coordinates of the outlier

.gen.coords <- function(type,
                        amount,
                        gridsize,
                        block.type){

  # determine the grid
  grid.coord <- expand.grid("x" = 1:gridsize[1], "y" = 1:gridsize[2])
  # determine the gridsize
  n <- nrow(grid.coord)

  ## isolated outliers: draw the coordinates at random
  if(type == "isolated"){
    row.ind <- sample(1:n, ceiling(n*amount))
    coords <- grid.coord[row.ind, ]
  }

  ## block outliers: build a block of outliers (spatially aggregated outliers)
  if(type == "block"){

    # amount of outliers
    n.outlier <- ceiling(n * amount)

    if(block.type == "square"){ # nearly quadratic outlier blocks
      # determine the block size of the largest possible square block
      dim.block <- c(floor(sqrt(n.outlier)), floor(sqrt(n.outlier)))

      # change block dimension, if the grid in direction is smaller than the
      # largest possible square block
      if(gridsize[1] < dim.block[1]){
        dim.block[2] <- dim.block[2] + (dim.block[1] - gridsize[1])
        dim.block[1] <- gridsize[1]
      }
      if(gridsize[2] < dim.block[2]){
        dim.block[1] <- dim.block[1] + (dim.block[2] - gridsize[2])
        dim.block[2] <- gridsize[2]
      }

      # determine how much smaller the largest possible square block is compared to
      # the required number of outliers
      rest <- n.outlier - dim.block[1] * dim.block[2]

      # Determine a random grid point (s0) around which the block is to be constructed
      start.point <- grid.coord[sample(1:n, 1),]

      # 1. Determine the x-coordinates of the largest possible square block
      # with s0 as centred as possible
      dims <- sample(c(floor((dim.block[1]-1)/2), ceiling((dim.block[1]-1)/2)), 2)
      x.coords <- as.numeric((start.point[1] - dims[1])):as.numeric((start.point[1] + dims[2]))

      x.true.u <- sum(x.coords < 1)
      x.true.o <- sum(x.coords > gridsize[1])

      if(x.true.u > 0){
        x.coords <- x.coords + x.true.u
      }
      if(x.true.o > 0){
        x.coords <- x.coords - x.true.o
      }


      # 2. Determine the y-coordinates of the largest possible square block
      # with s0 as centred as possible
      dims <- sample(c(floor((dim.block[2]-1)/2), ceiling((dim.block[2]-1)/2)), 2)
      y.coords <- as.numeric((start.point[2] - dims[1])):as.numeric((start.point[2] + dims[2]))

      y.true.u <- sum(y.coords < 1)
      y.true.o <- sum(y.coords > gridsize[2])

      if(y.true.u > 0){
        y.coords <- y.coords + y.true.u
      }
      if(y.true.o > 0){
        y.coords <- y.coords - y.true.o
      }

      # Assemble the coordinates of the largest possible square block around s0
      # from the two components (1. & 2.)
      coords <- expand.grid(x.coords, y.coords)


      # Add the remaining grid points so that the outlier block contains exactly
      # n.outlier points; the block should be as square as possible
      if(rest > 0){
        dims.rest <- sample(c(floor(rest/2), ceiling(rest/2)), 2)

        if(min(x.coords) == 1 & max(x.coords) == gridsize[1]){
          dims.rest[2] <- rest
          dims.rest[1] <- 0
        }

        if(min(y.coords) == 1 & max(y.coords) == gridsize[2]){
          dims.rest[1] <- rest
          dims.rest[2] <- 0
        }


        if(!any(dims.rest == 0)){
          # x-direction
          if(min(x.coords) == 1){
            x.rest <- max(x.coords) + 1
          }
          if(max(x.coords) == gridsize[1]){
            x.rest <- min(x.coords) - 1
          }
          if(min(x.coords) != 1 & max(x.coords) != gridsize[1]){
            x.rest <- sample(c((max(x.coords) + 1), (min(x.coords) - 1)), 1)
          }
          if(dims.rest[1] == dim.block[2]){
            y.rest <- y.coords
          } else {
            r <- dim.block[2] - dims.rest[1]
            help <- sample(c(0,1), 1)

            if(help == 1){
              y.rest <- y.coords[-(1:r)]
            } else{
              y.rest <- y.coords[dim.block[2]:1]
              y.rest <- y.rest[-(1:r)]
            }
          }
          coords <- rbind(coords, cbind(Var1  = x.rest, Var2 = y.rest))

          # y-direction
          if(min(y.coords) == 1){
            y.rest <- max(y.coords) + 1
          }
          if(max(y.coords) == gridsize[2]){
            y.rest <- min(y.coords) - 1
          }
          if(min(y.coords) != 1 & max(y.coords) != gridsize[2]){
            y.rest <- sample(c((max(y.coords) + 1), (min(y.coords) - 1)), 1)
          }
          if(dims.rest[2] == dim.block[1]){
            x.rest <- x.coords
          } else {
            r <- dim.block[1] - dims.rest[2]
            help <- sample(c(0,1), 1)

            if(help == 1){
              x.rest <- x.coords[-(1:r)]
            } else{
              x.rest <- x.coords[dim.block[1]:1]
              x.rest <- x.rest[-(1:r)]
            }
          }
          coords <- rbind(coords, cbind(Var1  = x.rest, Var2 = y.rest))
        } else{
          dims.rest2 <- c(ceiling(rest/2), floor(rest/2))

          if(dims.rest[1] == 0){
            if(min(y.coords) == 1){
              y.rest <- c(max(y.coords) + 1, max(y.coords) + 2)
            }
            if(max(y.coords) == gridsize[2]){
              y.rest <- c(min(y.coords) - 1, min(y.coords) - 2)
            }
            if(min(y.coords) != 1 & max(y.coords) != gridsize[2]){
              y.rest <- c((max(y.coords) + 1), (min(y.coords) - 1))
            }

            if(dims.rest2[1] == dim.block[1]){
              x.rest <- x.coords
            } else {
              r1 <- dim.block[1] - dims.rest2[1]
              r2 <- dim.block[1] - dims.rest2[2]
              help <- sample(c(0,1), 1)

              if(help == 1){
                x.rest <- list(x.coords[-(1:r1)], x.coords[-(1:r2)])
              } else{
                x.rev <- x.coords[dim.block[1]:1]
                x.rest <- list(x.rev[-(1:r1)], x.rev[-(1:r2)])
              }
            }

            coords <- rbind(coords, cbind(Var1  = unlist(x.rest[1]), Var2 = y.rest[1]), cbind(Var1 = unlist(x.rest[2]), Var2 = y.rest[2]))

          }

          if(dims.rest[2] == 0){
            if(min(x.coords) == 1){
              x.rest <- c(max(x.coords) + 1, max(x.coords) + 2)
            }
            if(max(x.coords) == gridsize[2]){
              x.rest <- c(min(x.coords) - 1, min(x.coords) - 2)
            }
            if(min(x.coords) != 1 & max(x.coords) != gridsize[2]){
              x.rest <- c((max(x.coords) + 1), (min(x.coords) - 1))
            }

            if(dims.rest2[1] == dim.block[2]){
              y.rest <- y.coords
            } else {
              r1 <- dim.block[2] - dims.rest2[1]
              r2 <- dim.block[2] - dims.rest2[2]
              help <- sample(c(0,1), 1)

              if(help == 1){
                y.rest <- list(y.coords[-(1:r1)], y.coords[-(1:r2)])
              } else{
                y.rev <- y.coords[dim.block[1]:1]
                y.rest <- list(y.rev[-(1:r1)], y.rev[-(1:r2)])
              }
            }

            coords <- rbind(coords, cbind(Var1  = x.rest[1], Var2 = unlist(y.rest[1])), cbind(Var1 = x.rest[2], Var2 = unlist(y.rest[2])))

          }
        }
      }
    }

    if(block.type == "random"){ # a random outlier block

      # Determine a random grid point (s0) around which the block is to be constructed
      center_idx <- sample(1:nrow(grid.coord), 1)
      center <- grid.coord[center_idx, ]

      # Find possible dimensions for the block
      possible_dims <- list()
      for (w in 1:n.outlier) {
        h <- floor(n.outlier / w)
        possible_dims[[length(possible_dims) + 1]] <- c(w, h)
      }
      # get rid of the ones that don't fit!
      possible_dims <- lapply(possible_dims, function(d){
        if(d[1] <= gridsize[1] & d[2] <= gridsize[2]) return(d)})
      possible_dims <- do.call(rbind, possible_dims)
      possible_dims <- rbind(possible_dims, cbind(possible_dims[,2], possible_dims[,1]))
      possible_dims <- unique(possible_dims)

      # choose randomly on dimnesion
      dims <- possible_dims[sample(nrow(possible_dims), 1),]
      width <- dims[1]
      height <- dims[2]

      # Half the width/height for a “centered” position of the random starting point
      half_w <- floor(width / 2)
      half_h <- floor(height / 2)

      # Boundaries of the grid
      min_x <- min(grid.coord$x)
      max_x <- max(grid.coord$x)
      min_y <- min(grid.coord$y)
      max_y <- max(grid.coord$y)

      # Calculate block boundaries
      x_min <- center$x - half_w
      x_max <- center$x + (width - half_w - 1)

      y_min <- center$y - half_h
      y_max <- center$y + (height - half_h - 1)


      # If the block is outside the grid → move it
      if(x_min < min_x){
        shift <- min_x - x_min
        x_min <- x_min + shift
        x_max <- x_max + shift
      }
      if(x_max > max_x) {
        shift <- x_max - max_x
        x_min <- x_min - shift
        x_max <- x_max - shift
      }

      if(y_min < min_y) {
        shift <- min_y - y_min
        y_min <- y_min + shift
        y_max <- y_max + shift
      }
      if(y_max > max_y) {
        shift <- y_max - max_y
        y_min <- y_min - shift
        y_max <- y_max - shift
      }
      # Select points in the block
      block <- subset(grid.coord, grid.coord$x >= x_min & grid.coord$x <= x_max &
                        grid.coord$y >= y_min & grid.coord$y <= y_max)

      # Add the missing observations so that the block sizes correspond to the number of outliers
      if(nrow(block) < n.outlier){
        diff.n <- n.outlier - nrow(block)

        if(x_min == 1 & x_max == gridsize[1]){
          if(y_min == 1){
            y.rest <- y_max + 1
          }
          if(y_max == gridsize[2]){
            y.rest <- y_min - 1
          }
          if(y_min != 1 & y_max != gridsize[2]){
            y.rest <- sample(c(y_max+1, y_min -1), 1)
          }
          x.rest.start <- sample(block$x, 1)
          x.rest <- x.rest.start:(x.rest.start+diff.n-1)
          if(max(x.rest) > gridsize[1]){
            x.rest <- x.rest - (max(x.rest) - gridsize[1])
          }
          block <- rbind(block, cbind("x" = x.rest, "y" = rep(y.rest, diff.n)))
        }
        if(y_min == 1 & y_max == gridsize[2]){
          if(x_min == 1){
            x.rest <- x_max + 1
          }
          if(x_max == gridsize[1]){
            x.rest <- x_min - 1
          }
          if(x_min != 1 & x_max != gridsize[1]){
            x.rest <- sample(c(x_max+1, x_min -1), 1)
          }
          y.rest.start <- sample(block$y, 1)
          y.rest <- y.rest.start:(y.rest.start+diff.n-1)
          if(max(y.rest) > gridsize[2]){
            y.rest <- y.rest - (max(y.rest) - gridsize[2])
          }
          block <- rbind(block, cbind("x" = rep(x.rest, diff.n), "y" = y.rest))
        }
        if(!(x_min == 1 & x_max == gridsize[1]) & !(y_min == 1 & y_max == gridsize[2])){

          direc <- sample(c("x", "y"), 1)

          if(direc == "y"){
            if(y_min == 1){
              y.rest <- y_max + 1
            }
            if(y_max == gridsize[2]){
              y.rest <- y_min - 1
            }
            if(y_min != 1 & y_max != gridsize[2]){
              y.rest <- sample(c(y_max+1, y_min -1), 1)
            }
            x.rest.start <- sample(block$x, 1)
            x.rest <- x.rest.start:(x.rest.start+diff.n-1)
            if(max(x.rest) > gridsize[1]){
              x.rest <- x.rest - (max(x.rest) - gridsize[1])
            }
            block <- rbind(block, cbind("x" = x.rest, "y" = rep(y.rest, diff.n)))
          }
          if(direc == "x"){
            if(x_min == 1){
              x.rest <- x_max + 1
            }
            if(x_max == gridsize[1]){
              x.rest <- x_min - 1
            }
            if(x_min != 1 & x_max != gridsize[1]){
              x.rest <- sample(c(x_max+1, x_min -1), 1)
            }
            y.rest.start <- sample(block$y, 1)
            y.rest <- y.rest.start:(y.rest.start+diff.n-1)
            if(max(y.rest) > gridsize[2]){
              y.rest <- y.rest - (max(y.rest) - gridsize[2])
            }
            block <- rbind(block, cbind("x" = rep(x.rest, diff.n), "y" = y.rest))
          }
        }
      }
      coords <- block
      if(nrow(coords) != n.outlier) stop("Block size does not match")
    }

    if(block.type == "rectangle"){# a block that is as long as possible

      # Determine a random grid point (s0) around which the block is to be constructed
      center_idx <- sample(1:nrow(grid.coord), 1)
      center <- grid.coord[center_idx, ]

      # Find possible dimensions for the block
      possible_dims <- list()
      for (w in 1:n.outlier) {
        h <- floor(n.outlier / w)
        possible_dims[[length(possible_dims) + 1]] <- c(w, h)
      }
      # get rid of the ones that don't fit!
      possible_dims <- lapply(possible_dims, function(d){
        if(d[1] <= gridsize[1] & d[2] <= gridsize[2]) return(d)})
      possible_dims <- do.call(rbind, possible_dims)
      possible_dims <- rbind(possible_dims, cbind(possible_dims[,2], possible_dims[,1]))
      possible_dims <- unique(possible_dims)

      # Choose rectangles that are as long as possible
      diff.dims <- abs(possible_dims[,1]-possible_dims[,2])
      possible_dims <- possible_dims[which(diff.dims == max(diff.dims)) ,]

      # Select a random dimension
      dims <- possible_dims[sample(nrow(possible_dims), 1),]
      width <- dims[1]
      height <- dims[2]

      # Half the width/height for a “centered” position of the random starting point
      half_w <- floor(width / 2)
      half_h <- floor(height / 2)

      # Boundaries of the grid
      min_x <- min(grid.coord$x)
      max_x <- max(grid.coord$x)
      min_y <- min(grid.coord$y)
      max_y <- max(grid.coord$y)

      # Calculate block boundaries
      x_min <- center$x - half_w
      x_max <- center$x + (width - half_w - 1)

      y_min <- center$y - half_h
      y_max <- center$y + (height - half_h - 1)


      # If the block is outside the grid → move it
      if(x_min < min_x) {
        shift <- min_x - x_min
        x_min <- x_min + shift
        x_max <- x_max + shift
      }
      if(x_max > max_x) {
        shift <- x_max - max_x
        x_min <- x_min - shift
        x_max <- x_max - shift
      }

      if(y_min < min_y) {
        shift <- min_y - y_min
        y_min <- y_min + shift
        y_max <- y_max + shift
      }
      if(y_max > max_y) {
        shift <- y_max - max_y
        y_min <- y_min - shift
        y_max <- y_max - shift
      }
      # Select points in the block
      block <- subset(grid.coord, grid.coord$x >= x_min & grid.coord$x <= x_max &
                        grid.coord$y >= y_min & grid.coord$y <= y_max)

      # Add the missing observations so that the block sizes correspond to the number of outliers
      if(nrow(block) < n.outlier){
        diff.n <- n.outlier - nrow(block)

        if(x_min == 1 & x_max == gridsize[1]){
          if(y_min == 1){
            y.rest <- y_max + 1
          }
          if(y_max == gridsize[2]){
            y.rest <- y_min - 1
          }
          if(y_min != 1 & y_max != gridsize[2]){
            y.rest <- sample(c(y_max+1, y_min -1), 1)
          }
          x.rest.start <- sample(block$x, 1)
          x.rest <- x.rest.start:(x.rest.start+diff.n-1)
          if(max(x.rest) > gridsize[1]){
            x.rest <- x.rest - (max(x.rest) - gridsize[1])
          }
          block <- rbind(block, cbind("x" = x.rest, "y" = rep(y.rest, diff.n)))
        }
        if(y_min == 1 & y_max == gridsize[2]){
          if(x_min == 1){
            x.rest <- x_max + 1
          }
          if(x_max == gridsize[1]){
            x.rest <- x_min - 1
          }
          if(x_min != 1 & x_max != gridsize[1]){
            x.rest <- sample(c(x_max+1, x_min -1), 1)
          }
          y.rest.start <- sample(block$y, 1)
          y.rest <- y.rest.start:(y.rest.start+diff.n-1)
          if(max(y.rest) > gridsize[2]){
            y.rest <- y.rest - (max(y.rest) - gridsize[2])
          }
          block <- rbind(block, cbind("x" = rep(x.rest, diff.n), "y" = y.rest))
        }
        if(!(x_min == 1 & x_max == gridsize[1]) & !(y_min == 1 & y_max == gridsize[2])){

          direc <- sample(c("x", "y"), 1)

          if(direc == "y"){
            if(y_min == 1){
              y.rest <- y_max + 1
            }
            if(y_max == gridsize[2]){
              y.rest <- y_min - 1
            }
            if(y_min != 1 & y_max != gridsize[2]){
              y.rest <- sample(c(y_max+1, y_min -1), 1)
            }
            x.rest.start <- sample(block$x, 1)
            x.rest <- x.rest.start:(x.rest.start+diff.n-1)
            if(max(x.rest) > gridsize[1]){
              x.rest <- x.rest - (max(x.rest) - gridsize[1])
            }
            block <- rbind(block, cbind("x" = x.rest, "y" = rep(y.rest, diff.n)))
          }
          if(direc == "x"){
            if(x_min == 1){
              x.rest <- x_max + 1
            }
            if(x_max == gridsize[1]){
              x.rest <- x_min - 1
            }
            if(x_min != 1 & x_max != gridsize[1]){
              x.rest <- sample(c(x_max+1, x_min -1), 1)
            }
            y.rest.start <- sample(block$y, 1)
            y.rest <- y.rest.start:(y.rest.start+diff.n-1)
            if(max(y.rest) > gridsize[2]){
              y.rest <- y.rest - (max(y.rest) - gridsize[2])
            }
            block <- rbind(block, cbind("x" = rep(x.rest, diff.n), "y" = y.rest))
          }
        }
      }
      coords <- block
      if(nrow(coords) != n.outlier) stop("Block size does not match")
    }
  }
  return(coords)
}

