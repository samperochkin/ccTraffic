#' Specify a fixed effect.
#'
#' @param prior Prior distribution on the (beta) coefficient. Only the gaussian prior is implemented.
#' @examples
#' fixedEffect()
#' @export
fixedEffect <- function(beta_prec = .01, is_factor = F, beta_mean = 0){
  list(beta_prec = beta_prec, is_factor = is_factor, beta_mean = beta_mean)
}
