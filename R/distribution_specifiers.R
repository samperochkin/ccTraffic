#' Specify an exponential (penalized complexity) prior distribution.
#'
#' PC prior distribution; see \insertCite{Simpson:2017;textual}{bayesEpi}. The parameters are such that \eqn{P(\sigma > u) = \alpha}, where \eqn{\sigma = -2log(\theta)}.
#'
#' @importFrom Rdpack reprompt
#' @references{
#'   \insertRef{Simpson:2017}{bayesEpi}
#' }
#' @param alpha A probability.
#' @param u A value greater than 0.
#' @return A list specifying a pc prep prior distribution.
#' @examples
#' pc_prec_prior(alpha = .5, u = .01)
#' @export
pc_prec_prior <- function(alpha = .5, u = .001){
  list(type="pc_prec", params=mget(names(formals()),sys.frame(sys.nframe())))
}
