###########
## Function for building differences/lag vectors between all data points

.differences.lags <- function(data){

  locations <- data[,2:3]
  z <-  data[,1]
  z <-  z-mean(z, na.rm = TRUE)
  n <-  nrow(data)

  index <-  sort(rep(1:n,n))
  splags <-  c()
  differences <-  c()
  origin.pts <-  matrix(data = NA, nrow = 0, ncol = 2)

  for(i in 1:n){
    cpt.mat <-  cbind(rep(locations[i,1], n), rep(locations[i,2], n))
    origin.pts <-  rbind(origin.pts, cpt.mat)
    lags.d <-  -1*cbind(locations[i,1] - locations[,1], locations[i,2] - locations[,2])
    splags <-  rbind(splags, lags.d)
    differences <-  c(differences, (z[i] - z))
  }

  x.coord <- origin.pts[,1]
  y.coord <- origin.pts[,2]
  x.lag <- splags[,1]
  y.lag <- splags[,2]
  diff <- differences
  res <- cbind(x.coord, y.coord, index, x.lag, y.lag, diff)
  row.names(res) <- NULL
  return(res)
}

###########
## Function for filtering the data of interest, i.e. data which lag vector is an element of lagmat

.raw.data.sub <- function(rawdata,
                          lagmat){

  nlags <-  dim(lagmat)[1]
  good <-  c()
  for(i in 1:nlags){
    gd <- which(rawdata[,"x.lag"] == lagmat[i,1] & rawdata[,"y.lag"] == lagmat[i,2])
    good <- c(good, gd)
  }

  return(rawdata[good,])
}

###########
## Function for calculation of the Matheron semivariogram estimator for all lags in lagmat

.matheron.test <- function(rawdata,
                           lagmat){

  nlags <-  dim(lagmat)[1]
  n.bin <-  c()
  gamma.hat <-  c()


  for(i in 1:nlags){
    clag <-  lagmat[i,]
    good <-  which((rawdata[,"x.lag"] == clag[1] & rawdata[,"y.lag"] == clag[2])|(rawdata[,"x.lag"] == -clag[1] & rawdata[,"y.lag"] == -clag[2]))
    n.clag <-  length(good)
    n.bin <-  c(n.bin, n.clag)
    gh <-  sum(rawdata[good,6]^2, na.rm = TRUE)/(2*n.clag)
    gamma.hat <-  c(gamma.hat, gh)
  }

  gamma.est <-  cbind(lagmat, gamma.hat, n.bin)
  colnames(gamma.est) <-  c("lag.x","lag.y","gamma.hat", "n.bin")
  return(gamma.est)
}

###########
## Function for calculation of the Genton semivariogram estimator for all lags in lagmat

.genton.test <- function(rawdata,
                         lagmat){

  nlags <-  dim(lagmat)[1]
  n.bin <-  c()
  gamma.hat <-  c()


  for(i in 1:nlags){
    clag <-  lagmat[i,]
    good <-  which((rawdata[,"x.lag"] == clag[1] & rawdata[,"y.lag"] == clag[2])|(rawdata[,"x.lag"] == -clag[1] & rawdata[,"y.lag"] == -clag[2]))
    n.clag <-  length(good)
    n.bin <-  c(n.bin, n.clag)
    gh <-  0.5*Qn(rawdata[good,6], na.rm = TRUE)^2
    gamma.hat <-  c(gamma.hat, gh)
  }

  gamma.est <-  cbind(lagmat, gamma.hat, n.bin)
  colnames(gamma.est) <-  c("lag.x","lag.y","gamma.hat", "n.bin")
  return(gamma.est)
}

###########
## Function for calculation of the MCD.diff semivariogram estimator for al lags in lagmat

.mcd.test <- function(vec.data,
                      lagmat,
                      ...){

  gamma.list <- lapply(vec.data, function(X) 0.5*diag(covMcd(t(X), use.correction = FALSE)$cov))#, ...)$cov))

  n.vec <- unlist(lapply(vec.data, ncol))
  n.vec <- unlist(lapply(1:length(vec.data), function(x) rep(n.vec[x], nrow(vec.data[[x]]))))

  lags <- unlist(lapply(vec.data, rownames))
  lags <- sapply(lags, function(x) strsplit(x, split = "\\."))
  lags <- do.call(rbind, lags)

  gamma.est <- data.frame(lags, unlist(gamma.list), n.vec)
  rownames(gamma.est) <- NULL
  colnames(gamma.est) <-  c("lag.x","lag.y","gamma.hat", "n.bin")

  # rearrange according to lagmat
  gamma.est$lag <- paste0(gamma.est$lag.x, ".", gamma.est$lag.y)
  lagmat <- as.data.frame(lagmat)
  lagmat$lag <- paste0(lagmat[,1], ".", lagmat[,2])
  gamma.est <- apply(lagmat, 1, function(x){
    gamma.est[which(gamma.est$lag == x["lag"]),]
  })
  gamma.est <- do.call(rbind, gamma.est)
  gamma.est <- gamma.est[,-5]

  return(gamma.est)
}

###########
## Function for checking if two values are equal up to a given tolerance

.is.almost.zero <- function(x, tol = .Machine$double.eps) {
  abs(x) < tol
}


###########
## Function to bulid the vectors for the MCD.diff semivariogram estimator

.rawdata.to.vec <- function(rawdata){

  vecs <- list()

  # for each direction one element
  split.rawdata <- split(rawdata, rawdata$lag.norm)

  n.l <- length(split.rawdata) # amount of directions

  for(l in 1:n.l){ # for each direction
    rawdata.l <- split.rawdata[[l]]

    # for each lag in this direction one element
    rawdata.l.split <- split(rawdata.l, rawdata.l$lag)

    # build the vectors
    if(length(rawdata.l.split) > 1){
      vec.l <- c()
      for(i in 1:nrow(rawdata.l.split[[1]])){
        ind.vec <- sapply(2:(length(rawdata.l.split)), function(ind){
          is.element(rawdata.l.split[[1]]$coord[i], rawdata.l.split[[ind]]$coord)
        })
        if(all(ind.vec)){
          vec.l <- cbind(vec.l, sapply(1:length(rawdata.l.split), function(ind){
            rawdata.l.split[[ind]][which(rawdata.l.split[[ind]]$coord == rawdata.l.split[[1]]$coord[i]), "diff"]
          }))
        }
      }
    } else{
      vec.l <- matrix(rawdata.l.split[[1]]$diff, nrow = 1)
    }
    rownames(vec.l) <- names(rawdata.l.split)
    vecs[[l]] <- vec.l
  }

  return(vecs)
}


