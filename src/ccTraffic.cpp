// Overdispersed case crossover.
#include <TMB.hpp>

/* prior */
template <class Type>
Type log_prior(Type theta, vector<Type> hypers)
{
  Type phi = -log(hypers(0)) / hypers(1);
  return log(0.5 * phi) - phi * exp(-0.5*theta) - 0.5*theta;
}


template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(count);
  DATA_IVECTOR(case_day);
  DATA_IMATRIX(control_days);

  DATA_MATRIX(X);

  DATA_VECTOR(beta_prec);
  DATA_VECTOR(theta_hypers);

  DATA_IVECTOR(z_pos);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(z);
  PARAMETER_VECTOR(theta);

  int beta_dim = beta.size();
  int z_dim = z.size();
  int theta_dim = theta.size();

  int eta_dim = 0;
  if(beta_dim != 0){
    eta_dim += X.col(0).size();
  }

  vector<Type> eta(eta_dim);
  eta.setZero();

  if(beta_dim != 0) eta += X*beta;
  for(int i = 0; i < eta_dim; i++){
    // if(z_pos(i) == 0) continue;
    eta(i) += z(z_pos(i)-1);
  }


  /*--------------------------------------------------------------------------*/
  /* LOG-LIKELIHOOD --------------------------------------------------------- */
  /*--------------------------------------------------------------------------*/
  Type log_likelihood = 0.0;
  Type log_hazard_ratio_sum;

  int n_case_day = case_day.size();
  int n_control_days = control_days.row(0).size();

  for (int i = 0;i<n_case_day;i++) {
    log_hazard_ratio_sum = 0.0;
    for(int j = 0;j<n_control_days;j++) {
      if(control_days(i,j) == 0) continue;
      log_hazard_ratio_sum = logspace_add(log_hazard_ratio_sum, eta(control_days(i,j) - 1) - eta(case_day(i) - 1));
    }
    log_likelihood -= count(i) * log_hazard_ratio_sum;
  }
  REPORT(log_likelihood);
  // Rcout << "ll : " << log_likelihood << "\n";


  /*--------------------------------------------------------------------------*/
  /* LOG LIKELIHOOD BETA -----------------------------------------------------*/
  /*--------------------------------------------------------------------------*/
  Type log_pi_beta = 0;
  for(int i=0;i<beta_dim;i++) log_pi_beta += dnorm(beta(i), Type(0), 1/sqrt(beta_prec(i)), true);
  REPORT(log_pi_beta);
  // Rcout << "lbeta : " << log_pi_beta << "\n";


  /*--------------------------------------------------------------------------*/
  /* LOG LIKELIHOOD Z --------------------------------------------------------*/
  /*--------------------------------------------------------------------------*/
  /* WATCH OUT HERE IF MORE THAN ONE THETA PARAMETER */
  Type log_pi_z = 0;
  for(int i=0; i<theta_dim; i++) log_pi_z += dnorm(z, Type(0), exp(-theta(i)/2), true).sum();
  REPORT(log_pi_z);
  // Rcout << "lz : " << log_pi_z << "\n";


  /*--------------------------------------------------------------------------*/
  /* LOG PRIOR FOR THETA -----------------------------------------------------*/
  /*--------------------------------------------------------------------------*/
  Type log_prior_theta = 0;
  int k = 0;
  vector<Type> hypers_i(2);
  for (int i=0; i<theta_dim; i++){
    for(int j=0;j<2;j++) hypers_i(j) = theta_hypers(k+j);
    log_prior_theta += log_prior(theta(i), hypers_i);
    k += 2;
  }
  REPORT(log_prior_theta);
  // Rcout << "ltheta : " << log_prior_theta << "\n";


  Type nll = -log_likelihood - log_pi_beta - log_pi_z - log_prior_theta;
  return nll;
}
