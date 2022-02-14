#' Fit a model
#'
#' Function used for fitting the desired case crossover model.
#'
#' @param model A model (list of info, see template).
#' @param data The dataset on which the model is fitted. It should contain all the variables involved in the model.
#' @return A fitted model.
#' @export
#' @useDynLib ccTraffic
#' @import TMB
#' @import aghq
#' @importFrom magrittr %>%
#' @importFrom purrr map
fitModel <- function(model, data){

  data <- as.data.frame(data)
  checkups(model, data)

  # 1) case_day, control_days. 2) z_init, z_pos. 3) X, beta_prec
  list2env(getCaseControl(data, model), envir = environment())
  list2env(setOverdispersion(data, model, case_day, control_days), envir = environment())
  list2env(constructX(model, data), envir = environment())

  # hyperparameters
  theta_hypers <- model$overdispersion$theta_prior$params %>% unlist # *** apply if more than one theta

  # Model fit ---------------------------------------------------------------
  tmb_data <- list(count = data[case_day, model$count_variable],
                   case_day = case_day, control_days = control_days,
                   X = X,
                   beta_prec = beta_prec, theta_hypers = theta_hypers,
                   z_pos = z_pos)

  theta_init <- getPriorInit(model)
  parameters <- list(beta = rep(0,ncol(tmb_data$X)), z = z_init, theta = theta_init)

  obj <- TMB::MakeADFun(tmb_data, parameters, random = c("beta","z"), DLL="ccTraffic", hessian=T)

  quad <- aghq::marginal_laplace_tmb(obj, model$aghq_input$k, theta_init, control = model$aghq_input$control)
  list(quad = quad, obj = obj, model = model)
}
