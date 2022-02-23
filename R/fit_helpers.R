# Helpers -----------------------------------------------------------------

# Checkups
checkups <- function(model, data){
  findVariables(model, data)
  checkNA(model, data)
}

# check if all variables are actually in the data.
findVariables <- function(model, data){
  var_names <- c(model$count_variable, names(model$fixed), names(model$offset))
  var_names <- var_names[!grepl(":", var_names)]
  if(!all(var_names %in% names(data))){
    stop("One of the names provided for the count_variable or fixed effects could not be matched in the data provided.")
  }
}

# check for NA.
checkNA <- function(model, data){
  var_names <- c(model$count_variable, names(model$fixed), names(model$offset))
  var_names <- var_names[!grepl(":", var_names)]
  if(any(is.na(data[,var_names]))){
    stop("The dataset contains NA... please remove them before fitting the model")
  }
  NULL
}


# Create the case_day vector and the corresponding control_days matrix.
#' @import purrr
getCaseControl <- function(data, model){

  strata_id <- data[, model$strata_variable]
  u_strata_id <- unique(strata_id)
  strata_id <- match(strata_id, u_strata_id)
  case_day <- which(data[, model$count_variable] > 0)

  stratum <- split(1:nrow(data), strata_id)
  max_len <- max(sapply(stratum, length)) - 1
  control_days <- lapply(case_day, function(c_day_id){
    con <- setdiff(stratum[[strata_id[c_day_id]]], c_day_id)
    con <- c(con, rep(0, max_len-length(con)))
    con
  }) %>% Reduce(f="rbind")

  # filter out case day with no control days
  keep <- apply(control_days != 0, 1, any)
  case_day <- case_day[keep]
  control_days <- control_days[keep,,drop=F]

  if(any(rowSums(control_days) == 0)) stop("Error in selecting the control days")

  list(data = data, case_day = case_day, control_days = control_days)
}


# Identify where to remove (not to put) overdispersion terms
setOverdispersion <- function(data, model, case_day, control_days){

  clus_vars <- model$overdispersion$cluster_variables
  clus_values <- data[, clus_vars]
  u_clus_values <- lapply(1:2, function(k) unique(clus_values[,k]))
  clus_ids <- sapply(1:2, function(k) match(clus_values[,k], u_clus_values[[k]]))

  colnames(clus_ids) <- paste0(clus_vars, "_", "id")
  rownames(clus_values) <- NULL

  u_ids <- unique(clus_ids)
  u_ids <- u_ids[order(u_ids[,1],u_ids[,2]),]
  # first_of_each <- sapply(unique(u_ids[,1]), function(id) which(u_ids[,1] == id)[1])
  # u_ids <- u_ids[-first_of_each,]

  z_init <- rep(0,  nrow(u_ids))
  z_pos <- match(apply(clus_ids, 1, paste, collapse = "__"),
                 apply(u_ids, 1, paste, collapse = "__"))
  # z_pos[is.na(z_pos)] <- 0
  if(any(is.na(z_pos))) return("Something went wrong assigning overdispersion terms... Call me.")

  model$overdispersion$cluster_mat <- cbind(sapply(1:2, function(k) u_clus_values[[k]][u_ids[,k]]), u_ids)
  colnames(model$overdispersion$cluster_mat)[1:2] <- c(clus_vars)

  list(model = model, z_init = z_init, z_pos = z_pos)
}

constructX <- function(model, data){

  fixed_names <- names(model$fixed)
  is_factor <- sapply(model$fixed, "[[", "is_factor")
  is_interaction <- grepl(":", fixed_names)
  names(is_interaction) <- fixed_names

  # put interactions at the end
  fixed_names <- c(fixed_names[!is_interaction],fixed_names[is_interaction])
  is_factor <- c(is_factor[!is_interaction],is_factor[is_interaction])
  is_interaction <- c(is_interaction[!is_interaction],is_interaction[is_interaction])

  # some checkups
  if(any(is_factor & is_interaction)) stop("Some fixed effect is a both a factor (is_factor = T)
                                           and an interaction (specified with ':'). It's most likely
                                           that you did not mean to have is_factor = T.")

  interacting_variables <- unlist(strsplit(fixed_names[is_interaction], ":"))
  if(!all(interacting_variables %in% fixed_names)) stop("Variables involved in interactions must be
                                                        introduced as fixed effects on their own as well.")

  # design matrix
  X <- matrix(nrow = nrow(data), ncol = 0)
  beta_prec <- numeric(0)

  for(fixed_name in fixed_names){
    if(!is_factor[fixed_name] & !is_interaction[fixed_name]){
      x <- as.matrix(data[,fixed_name])
      colnames(x) <- fixed_name
      beta_len <- 1
    }else if(is_factor[fixed_name]){
      u <- unique(data[,fixed_name])
      model$fixed[[fixed_name]]$factor_labels <- u
      x <- matrix(0, nrow(data), length(u) - 1)
      x[cbind(1:nrow(x),match(data[,fixed_name], u[-1]))] <- 1
      colnames(x) <- paste0(fixed_name, "_", u[-1])
      model$fixed[[fixed_name]]$column_names <- colnames(x)
      beta_len <- length(u) - 1
    }else if(is_interaction[fixed_name]){
      interacting_variables <- unlist(strsplit(fixed_name, ":"))
      columns_list <- lapply(interacting_variables, function(int_var){
        if(is_factor[int_var]){
          return(model$fixed[[int_var]]$column_names)
        }
        return(int_var)
      })
      column_pairs <- expand.grid(columns_list)
      x <- apply(column_pairs, 1, function(var){
        apply(X[,var],1,prod)
      })
      colnames(x) <- apply(column_pairs,1,paste,collapse=":")
      beta_len <- nrow(column_pairs)
      if(any(is_factor[interacting_variables])){
        labs <- lapply(interacting_variables, function(var){
          lab <- model$fixed[[var]]$factor_labels
          if(is.null(lab)) return(NA)
          lab
        })
        model$fixed[[fixed_name]]$factor_labels <- apply(expand.grid(labs), 1, paste, collapse = ":")
        model$fixed[[fixed_name]]$is_factor <- T
      }
    }else{
      stop("Something went wrong. Write to me.")
    }
    X <- cbind(X, x)
    beta_prec <- c(beta_prec, rep(model$fixed[[fixed_name]]$beta_prec, beta_len))
  }
  names(beta_prec) <- colnames(X)

  list(model = model, X = X, beta_prec = beta_prec)
}

# Compute initial theta parameter to be passed to aghq::quad.
getPriorInit <- function(model){
  priors <- list(model$overdispersion$theta_prior) # ************** NEED APPLY IF MORE THAN ONE?
  sapply(priors, function(prior){
    return(-2*log(-prior$params$u/log(prior$params$alpha)))
  })
}
