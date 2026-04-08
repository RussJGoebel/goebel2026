#' Build a row-normalised spatial adjacency matrix (W) from an sf grid
#'
#' @description
#' Constructs a sparse row-normalised adjacency matrix using queen contiguity
#' (i.e. shared edges or corners) via \code{spdep::poly2nb}. Suitable for
#' any polygon grid, regular or irregular.
#'
#' @param target_grid An sf POLYGON object representing the target grid.
#' @param zero.policy Logical passed to \code{spdep::nb2mat}. If TRUE, allows
#'   regions with no neighbours (returns zero row). Default TRUE.
#'
#' @return A sparse \code{Matrix::Matrix} of class \code{dgCMatrix} with
#'   nrow = ncol = nrow(target_grid). Row sums are 1 for all non-isolated cells.
#'
#' @examples
#' \dontrun{
#' W <- make_W_matrix(target_grid)
#' }
#'
#' @export
make_W_matrix <- function(target_grid, zero.policy = TRUE) {
  nb <- spdep::poly2nb(target_grid)
  W  <- Matrix::Matrix(
    spdep::nb2mat(nb, style = "W", zero.policy = zero.policy),
    sparse = TRUE
  )
  return(W)
}


#' Build a SAR precision matrix Q from an sf grid
#'
#' @description
#' Constructs the simultaneous autoregressive (SAR) precision matrix
#' \eqn{Q = (I - \rho W)^\top (I - \rho W)}, where \eqn{W} is the
#' row-normalised queen adjacency matrix of the grid. When \eqn{\rho = 1}
#' this yields an intrinsic SAR prior.
#'
#' @param target_grid An sf POLYGON object representing the target grid.
#' @param rho Spatial autoregression parameter in (0, 1]. Default 1 (intrinsic prior).
#' @param zero.policy Logical passed to \code{spdep::nb2mat}. Default TRUE.
#'
#' @return A sparse symmetric \code{Matrix::Matrix} of class \code{dsCMatrix}
#'   with nrow = ncol = nrow(target_grid).
#'
#' @details
#' The scaling parameter \eqn{\tau} is not included here — it is handled
#' via the \code{phi} argument in \code{fastblm::fit_fastblm} and
#' \code{fastblm::tune_cv}.
#'
#' @examples
#' \dontrun{
#' Q <- make_sar_precision(target_grid)
#' Q <- make_sar_precision(target_grid, rho = 0.9)
#' }
#'
#' @export
make_sar_precision <- function(target_grid, rho = 1, zero.policy = TRUE) {
  W       <- make_W_matrix(target_grid, zero.policy = zero.policy)
  IminusW <- Matrix::Diagonal(nrow(W)) - rho * W
  Q       <- Matrix::crossprod(IminusW)
  return(Q)
}
