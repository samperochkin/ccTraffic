#' Extract results from a fitted case-crossover model.
#' @param fit A fitted case-crossover model (obtained via `fitModel(model)`).
#' @param probs Probabilities at which to recover the quantiles of the posterior (cumulative) distribution of fixed and random effects.
#' @param M Number of replicates to draw for the posterior (to approximate it, the more the better).
#' @return Three summary data.frame's. One each for beta, z (overdispersion) and theta.
#' @export
getResults <- function(fit, probs = c(.01, .025, .05, .95, .975, .99), M = 1e4){

  quad_samples <- aghq::sample_marginal(fit$quad, M)
  model <- fit$model
  fixed_names <- names(model$fixed)

  #### beta output ####
  beta_output <- data.frame(parameter_type = "beta", variable_fname = names(fit$obj$env$data$beta_prec))
  beta_output$variable_name <- sapply(strsplit(as.character(beta_output$variable_fname), "_"), "[[", 1)
  beta_output$factor_label <- NA
  beta_output$factor_id <- as.integer(NA)

  fixed_factors <- c(fixed_names[sapply(model$fixed, "[[", "is_factor")])

  for(fixed_name in fixed_factors){
    ind <- which(beta_output$variable_name == fixed_name)
    beta_output$factor_id[ind] <- seq_along(ind) + 1
    beta_output$factor_label[ind] <- model$fixed[[fixed_name]]$factor_labels[beta_output$factor_id[ind]]
  }

  ind <- which(rownames(quad_samples$samps) == "beta")
  beta_output$mean <- rowMeans(quad_samples$samps[ind,,drop=F])
  beta_output$median <- apply(quad_samples$samps[ind,,drop=F],1,median)
  beta_output$sd <- apply(quad_samples$samps[ind,,drop=F],1,sd)
  beta_output <- cbind(beta_output,
                       t(apply(quad_samples$samps[ind,,drop=F],1,quantile, probs = probs)))
  colnames(beta_output)[ncol(beta_output) - (length(probs)-1):0] <- paste0("perc_", probs*100)
  rownames(beta_output) <- NULL


  #### z output ####
  z_output <- cbind(data.frame(parameter_type = "z"), model$overdispersion$cluster_mat)

  ind <- which(rownames(quad_samples$samps) == "z")
  z_output$mean <- rowMeans(quad_samples$samps[ind,,drop=F])
  z_output$median <- apply(quad_samples$samps[ind,,drop=F],1,median)
  z_output$sd <- apply(quad_samples$samps[ind,,drop=F],1,sd)
  z_output <- cbind(z_output,
                    t(apply(quad_samples$samps[ind,,drop=F],1,quantile, probs = probs)))
  colnames(z_output)[ncol(z_output) - (length(probs)-1):0] <- paste0("perc_", probs*100)
  rownames(z_output) <- NULL


  #### theta output ####
  theta_output <- data.frame(parameter_type = "theta",
                             variable_name = paste(model$overdispersion$cluster_variables, collapse = "_"),
                             mean = sapply(quad_samples$thetasamples, mean),
                             median = sapply(quad_samples$thetasamples, median),
                             sd = sapply(quad_samples$thetasamples, sd))

  theta_output <- cbind(theta_output,
                        do.call("rbind", lapply(quad_samples$thetasamples, function(samp) t(quantile(samp, probs = probs)))))
  colnames(theta_output)[ncol(theta_output) - (length(probs)-1):0] <- paste0("perc_", probs*100)

  return(list(beta = beta_output, z = z_output, theta = theta_output))
}

#' Compute posterior probabilities from a posterior distribution (fitted model).
#' @param fit A fitted case-crossover model (obtained via `fitModel(model)`).
#' @param beta_names Names of beta parameters concerned.
#' @param quantiles Quantiles at which to compute the probabilities.
#' @param M Number of replicates to draw for the posterior (to approximate the probability, the more the better).
#' @return A vector of (approximate) probabilities.
#' @export
getPostProb <- function(fit, beta_names, quantiles, M = 1e4){
  ind <- match(beta_names, names(fit$obj$env$data$beta_prec))
  probs <- rowMeans(aghq::sample_marginal(fit$quad, M)$samps[ind,,drop=F] < quantiles)
  names(probs) <- paste0("P(", beta_names, " < ", quantiles, ")")
  return(probs)
}
