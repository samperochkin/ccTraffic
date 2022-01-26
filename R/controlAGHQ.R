



#' Set aghq::quad parameters for model fitting.
#'
#' .
#'
#' So far, controlAGHQ() holds only one parameter, `k`, the number of quadrature points.
#'
#' @param ... See details.
#' @return A list of parameters for use by aghq::aghq.
#' @examples
#' controlAGHQ()
#' controlAGHQ(k=5)
#' @export
controlAGHQ <- function(...){

  params = list(...)

  # default cc parameters
  control <- list(k = 3)

  control[names(params)] <- params
  return(control)
}
