# Dev script for fitting a case-crossover model to traffic data
# Original author: Samuel Perreault (samuel.perreault[at]utoronto[dot]ca)
# This document might not be the original version.



# packages ----------------------------------------------------------------
library(ccTraffic)
library(ggplot2)



# synthetic data ----------------------------------------------------------
n_int <- 200 # number of intersections
n_cen <- 8 # number of census tracks

# Int and Cen variables
data <- data.frame(Int = sample(1e10,n_int), Cen = sample(n_cen,n_int,replace=T))

# add time dimension
time_grid <- expand.grid(k = 1:nrow(data), time = 1:8)
data <- cbind(data[time_grid$k,], time = time_grid$time)

# remove some rows to amke sure algorithm is robust to that
data <- data[-sample(nrow(time_grid), n_int),]

# create fake exposure (random)
data$exposure <- sample(0:2, nrow(data), replace = T)

# create overdispersion terms (per Cen and time)
z_id <- data[,cbind("Cen", "time")]
z_frame <- unique(z_id) # for later
z_id <- apply(z_id, 1, paste, collapse = "_") # paste the two columns to create unique ids
u_z_id <- unique(z_id) # collect unique ids
data$z_id <- match(z_id, u_z_id) # switch from weird x_y id to a single number
theta <- 4.5
sigma <- sqrt(exp(-theta)); sigma
z_effect <- rnorm(length(u_z_id), 0, sigma)
z_frame$z_true <- z_effect

# create effect of time
time_effect <- runif(max(data$time), 0, .5)

# generate counts from true model (with intercept of 1 that we cannot estimate)
data$count <- rpois(nrow(data), exp(1 + time_effect[data$time] - .05 * data$exposure + z_effect[data$z_id]))



# specify model -----------------------------------------------------------

# main model
model <- list(
  count_variable = "count",
  strata_variable = "Int",
  fixed = list("exposure" = fixedEffect(beta_prec = 0.01),
               "time" = fixedEffect(beta_prec = 0.01, is_factor = T)),
  overdispersion = list(cluster_variables = c("Cen", "time"),
                        theta_prior = pc_prec_prior(u = .1, alpha = .5)),
  control_aghq = controlAGHQ()
)

# alternate model
# model <- list(
#   count_variable = "count",
#   strata_variable = "Int",
#   fixed = list("exposure" = fixedEffect(beta_prec = 0.01),
#                "time" = fixedEffect(beta_prec = 0.01, is_factor = T),
#                "exposure:time" = fixedEffect(beta_prec = 0.01)),
#   overdispersion = list(cluster_variables = c("Cen", "time"),
#                         theta_prior = pc_prec_prior(u = .1, alpha = .5)),
#   control_aghq = controlAGHQ()
# )

# https://inla.r-inla-download.org/r-inla.org/doc/prior/pc.prec.pdf
# Voir eq. (1) pour la densite de theta et le sens de alpha et u.



# Fit model ---------------------------------------------------------------
fit <- fitModel(model, data, dll = "cc")
results <- getResults(fit)

# beta estimates
results$beta

# overdispersion estimates (easier to investigate via plots)
results$z

# estimate of theta (and sigma), the variance parameter for the overdispersion
matrix(c(results$theta$median, theta, sqrt(exp(-results$theta$median)), sigma),
       2, 2, dimnames = list(c("estimate", "true"), c("theta", "sigma")))

# how to get posterior probabilities
getPostProb(fit, beta_names = c("exposure"), quantiles = c(-.05))



# plots -------------------------------------------------------------------

# plots of beta estimates with credible intervals
# Not that we use time_effect[-1]-time_effect[1] to create the green points (true effect)
# true effects make sense only when there is no interaction...
ggplot(results$beta, aes(x = as.factor(factor_label), y=median)) +
  theme_bw() +
  geom_errorbar(aes(ymin=perc_2.5, ymax=perc_97.5), width=.1) +
  geom_point() +
  geom_abline(data=data.frame(variable_name = "exposure", true_value = -.05),
              aes(slope=0, intercept = true_value), lty = 2, col="green") +
  geom_point(data = data.frame(variable_name = "time", factor_label = as.factor(2:8), y = time_effect[-1] - time_effect[1]),
             aes(x = factor_label, y = y), col = "green") +
  facet_wrap(~variable_name, scales = "free_x")


# plots of z estimates with credible intervals
# Here we remove the effect of the overdispersion term at time=1 (not fitted)
z_frame1 <- z_frame[z_frame$time == 1,]
names(z_frame1)[3] <- "z_ref"
z_frame <- merge(z_frame, z_frame1[,-2], by = c("Cen"))
z_frame$z_adj <- z_frame$z_true - z_frame$z_ref

# actual plots
# if you have lots of census tracks, the following chunk will crash your R
# take a subset of, say, 10 census tracks...
u_Cen <- unique(results$z$Cen)
Cen_subset <- u_Cen[1:min(length(u_Cen), 10)]
ggplot(results$z[results$z$Cen %in% Cen_subset,], aes(x = time, y=median)) +
  theme_bw() +
  geom_point() +
  geom_errorbar(aes(ymin=perc_2.5, ymax=perc_97.5), width=.1) +
  geom_line(alpha=.25) +
  geom_point(data = z_frame[z_frame$Cen %in% Cen_subset,], aes(y = z_adj), col="green") +
  # geom_point(data = z_frame[z_frame$Cen %in% Cen_subset,], aes(y = z_true), col="green") +
  facet_wrap(~Cen, scales = "free_x")


