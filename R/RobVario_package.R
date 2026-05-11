#' @docType package
#' @name RobVario-package
#' @details
#' The package provides functions for analyzing the spatial dependency structure of spatial data on a regular grid with isolated as well as spatially aggregated outliers (outlier block).
#'
#' The functions \code{\link{variogram_est}} and \code{\link{variogram_est_general}} implement robust and non-robust directional variogram estimators. The main difference is that \code{\link{variogram_est}}
#' can be only used for estimation in the four directions south-north, east-west, southwest-northeast and southeast-northwest, while \code{\link{variogram_est_general}} can be used for estimation in an arbitrary direction.
#'
#' The function \code{\link{isotropy_test}} can be used to test for isotropy using either a subsampling or a blockpermutation approach. Thereby robust and non-robust test are implemented.
#'
#' Gaussian random fields with isolated or an outlier block can be simulated using the function
#' \code{\link{simulate_grf}}.
#'
#' For small grids the robust variogram estimators need correction factors. This can be simulated using the function \code{\link{simulate_correctionfactors}}.
#'
#' For more details on the models see the documentation of the functions.
#'
#' @references
#' - Genton, M. (1998). Highly robust variogram estimation. \emph{Mathematical Geology}, 30, 213-221. \doi{https://doi.org/10.1023/A:1021728614555}
#'
#' - Gierse, J., & Fried, R. (2025). Nonparametric directional variogram estimation in the presence of outlier blocks. \emph{Statistical Papers}, 66(134). \doi{https://doi.org/10.1007/s00362-025-01754-2}
#'
#' - Gierse, J. & Fried. R. (2026). EINFÜGEN
#'
#' - Guan, Y., Sherman, M. & Calvin, J. A. (2004). A nonparametric test for spatial isotropy using subsampling. \emph{J. Americ. Statist. Assoc.}, 99, 810-821. \doi{https://doi.org/10.1198/016214504000001150}
#'
#' - Matheron, G. (1962). Traité de géostatistique appliquée, Tome I. \emph{Mémoires du Bureau de Recherches Géologiques et Miniéres}, no. 14, Editions Technip, 333
#'
"_PACKAGE"
