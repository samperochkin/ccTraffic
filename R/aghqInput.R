



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
aghqInput <- function(...){

  params = list(...)

  # default cc parameters
  aghq_input <- list(k = 3, control = aghq::default_control_tmb())

  aghq_input[names(params)] <- params
  return(aghq_input)
}
