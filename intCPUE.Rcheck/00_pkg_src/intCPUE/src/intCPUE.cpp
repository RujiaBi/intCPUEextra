#define TMB_LIB_INIT R_init_intCPUE
#define EIGEN_DONT_PARALLELIZE
#include <TMB.hpp>
#define _USE_MATH_DEFINES
#include <cmath>

// List Of Sparse Matrices (LOSM): read a List from R into vector<SparseMatrix<> >
// following sdmTMB, taken from kaskr, https://github.com/kaskr/adcomp/issues/96
using namespace Eigen;
using namespace tmbutils;
template<class Type>
struct LOSM_t : vector<SparseMatrix<Type> > {
  LOSM_t(SEXP x){
    this->resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asSparseMatrix<Type>(sm);
    }
  }
};

// Bias corrected log-normal
template<class Type>
Type dlnorm_bc(const Type& x, const Type& meanlog, const Type& sdlog, int give_log = 0) {
    Type sd2 = sdlog * sdlog;
    Type adjusted_meanlog = meanlog - sd2 / 2;
    Type logres = dnorm(log(x), adjusted_meanlog, sdlog, true) - log(x);
	
    if (give_log)
        return logres;
    else
        return exp(logres);
}

// PC prios on range and sigma
template <class Type>
Type pc_prior_matern(        
  Type range, 
  Type sigma,   
  Type matern_range,   // r₀       — the user‐supplied “Pr(range<r₀)=range_prob” cutoff  
  Type matern_sigma,   // σ₀       — the user‐supplied “Pr(σ>σ₀)=SD_prob” cutoff
  Type range_prob,     // Pr(range<r₀)  
  Type sigma_prob,     // Pr(σ>σ₀)  
  int  give_log=0,     // if 1, return log‐density; if 0, return density  
  int  share_range=0   // if 1, omit the range‐prior term (useful in hierarchical models) 
){
  // dimension and half‐dimension
  Type d     = 2.;
  Type dhalf = d / 2.;

  // compute λ parameters
  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(sigma_prob) / matern_sigma;

  // log‐density of PC‐prior for range (unless share_range==1)
  Type range_ll = log(dhalf) 
                + log(lam1) 
                + (-1.0 - dhalf) * log(range)
                - lam1 * pow(range, -dhalf);

  // log‐density of PC‐prior for σ
  Type sigma_ll = log(lam2) 
                - lam2 * sigma;

  // combine
  Type penalty = sigma_ll;
  if(!share_range) penalty += range_ll;

  if(give_log)   return penalty;  
  else           return exp(penalty);
}


// Main TMB objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // Dimensions
  DATA_INTEGER(n_i); // Number of observations (stacked across all years)
  DATA_INTEGER(n_t); // Number of time-indices
  DATA_INTEGER(n_v); // Number of vessels (i.e., levels for the factor explaining overdispersion)
  DATA_INTEGER(n_f); // Number of flags (e.g., 0 = the reference, 1 = other 1, 2 = other 2, etc.)
  DATA_INTEGER(n_g); // Number of extrapolation-grid cells
  DATA_INTEGER(use_vessel_effect);
  DATA_INTEGER(use_pop_spatial);
  DATA_INTEGER(use_pop_spatiotemporal);
  DATA_INTEGER(use_pop_spatiotemporal_rw);
  DATA_INTEGER(use_pop_spatiotemporal_ar1);
  DATA_INTEGER(use_q_diffs_system);
  DATA_INTEGER(use_q_diffs_time);
  DATA_INTEGER(use_q_diffs_spatial);
  DATA_INTEGER(use_flag_sd);
  
  // Data
  DATA_VECTOR(b_i); // Response (Positive CPUE, biomass / effort) for each observation
  DATA_IVECTOR(e_i); // Response (0 or 1) for each observation
  DATA_IVECTOR(t_i); // Time index for each observation
  DATA_IVECTOR(v_i); // Vessel for each observation
  DATA_IVECTOR(f_i); // Flag (gear) type (e.g., 0 = the reference, 1 = other 1, 2 = other 2, etc.) for each observation  
  
  DATA_IMATRIX(has_tf); // n_t x (n_f-1), entries 0/1
  
  DATA_VECTOR(area_g); // Area for each extrapolation-grid cell
  
  // Projection matrices from knots s to data i or extrapolation-grid cells x
  DATA_SPARSE_MATRIX(A_is); // Project vertices to samples 
  DATA_SPARSE_MATRIX(A_gs); // Project vertices to integration points
  DATA_IMATRIX(Ais_ij);
  DATA_VECTOR(Ais_x);
  
  // PC priors on SPDE hyperparamters
  DATA_SCALAR(matern_range);    // r₀: e.g. diff(range(mesh$loc[,1]))/5
  DATA_SCALAR(range_prob);      // P(range < r₀)
  DATA_SCALAR(matern_sigma_0);  // σ₀ for omega
  DATA_SCALAR(matern_sigma_t);  // σ₀ for epsilon
  DATA_SCALAR(matern_sigma_flag);  // σ₀ for epsilon
  DATA_SCALAR(sigma_prob);      // P(σ > σ₀)
  
  // Catchability smooth objects
  DATA_INTEGER(has_smooths_catch);
  DATA_MATRIX(Xs_catch);
  DATA_STRUCT(Zs_catch, LOSM_t);
  DATA_IVECTOR(b_smooth_start_catch);

  // Population smooth objects (observation rows + projection rows)
  DATA_INTEGER(has_smooths_pop);
  DATA_MATRIX(Xs_pop_i);
  DATA_STRUCT(Zs_pop_i, LOSM_t);
  DATA_MATRIX(Xs_pop_g);
  DATA_STRUCT(Zs_pop_g, LOSM_t);
  DATA_IVECTOR(b_smooth_start_pop);

  // Aniso SPDE objects
  DATA_STRUCT(spde, spde_aniso_t);

  // Parameters  
  PARAMETER(ln_sd); // SD for Log-normal distribution
  PARAMETER_VECTOR(ln_sd_flag); // Flag-specific SDs for Log-normal distribution
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  
  PARAMETER(ln_range_1); // SPDE hyper-parameter
  PARAMETER(ln_sigma_0_1); // Time-constant SPDE hyper-parameter
  PARAMETER(ln_sigma_t_1); // Time-varying SPDE hyper-parameter
  PARAMETER(transf_rho_1); // Spatiotemporal AR1 correlation (working scale)
  
  PARAMETER(ln_range_2); // SPDE hyper-parameter
  PARAMETER(ln_sigma_0_2); // Time-constant SPDE hyper-parameter
  PARAMETER(ln_sigma_t_2); // Time-varying SPDE hyper-parameter
  PARAMETER(transf_rho_2); // Spatiotemporal AR1 correlation (working scale)

  // Vessel random effects
  PARAMETER_VECTOR(ves_v_1);
  PARAMETER_VECTOR(ves_v_2);
  PARAMETER(ves_ln_std_dev_1);
  PARAMETER(ves_ln_std_dev_2);
  
  // Fixed temporal effects 
  PARAMETER_VECTOR(yq_t_1); 
  PARAMETER_VECTOR(yq_t_2);
  
  // Spatial fields
  PARAMETER_VECTOR(omega_s_1); // Time-constant SPDE spatial field (over mesh nodes)
  PARAMETER_VECTOR(omega_s_2); // Time-constant SPDE spatial field (over mesh nodes)
  PARAMETER_MATRIX(epsilon_st_1); // Time-varying SPDE spatial field (over mesh nodes) across time [s,c]
  PARAMETER_MATRIX(epsilon_st_2); // Time-varying SPDE spatial field (over mesh nodes) across time [s,c]
  
  // Flag systematic differences
  PARAMETER_VECTOR(flag_f_1);
  PARAMETER_VECTOR(flag_f_2);
  PARAMETER(flag_ln_std_dev_1); 
  PARAMETER(flag_ln_std_dev_2); 
  
  // Flag temporal differences
  PARAMETER_MATRIX(flag_t_1); 
  PARAMETER_MATRIX(flag_t_2);
  PARAMETER(flag_t_ln_std_dev_1); 
  PARAMETER(flag_t_ln_std_dev_2); 

  // Flag spatial differences
  PARAMETER_MATRIX(flag_s_1); 
  PARAMETER_MATRIX(flag_s_2);
  PARAMETER(ln_sigma_flag_1); 
  PARAMETER(ln_sigma_flag_2);
  
  // Catchability smoothers (two components: encounter=0, positive=1)
  PARAMETER_MATRIX(bs_catch);
  PARAMETER_MATRIX(b_smooth_catch);
  PARAMETER_MATRIX(ln_smooth_sigma_catch);

  // Population smoothers (two components: encounter=0, positive=1)
  PARAMETER_MATRIX(bs_pop);
  PARAMETER_MATRIX(b_smooth_pop);
  PARAMETER_MATRIX(ln_smooth_sigma_pop);
  
  // Global variables
  Type nll = 0;
  Type nll_prior = 0;
  Type nll_penalty = 0;
  Type sd = exp(ln_sd); 
  vector<Type> sd_flag(n_f);
  for (int f = 0; f < n_f; f++) {
    sd_flag(f) = exp(ln_sd_flag(f));
  }
  Type range_1      = exp(ln_range_1);
  Type sigma_0_1    = exp(ln_sigma_0_1);
  Type sigma_t_1    = exp(ln_sigma_t_1);
  Type sigma_flag_1 = exp(ln_sigma_flag_1);  
  Type ves_std_dev_1 = exp(ves_ln_std_dev_1); 
  Type flag_std_dev_1 = exp(flag_ln_std_dev_1); 
  Type flag_t_std_dev_1 = exp(flag_t_ln_std_dev_1); 
  Type rho_1 = Type(2.0) / (Type(1.0) + exp(-Type(2.0) * transf_rho_1)) - Type(1.0);
  Type range_2      = exp(ln_range_2);
  Type sigma_0_2    = exp(ln_sigma_0_2);
  Type sigma_t_2    = exp(ln_sigma_t_2);
  Type sigma_flag_2 = exp(ln_sigma_flag_2); 
  Type ves_std_dev_2 = exp(ves_ln_std_dev_2);
  Type flag_std_dev_2 = exp(flag_ln_std_dev_2);
  Type flag_t_std_dev_2 = exp(flag_t_ln_std_dev_2); 
  Type rho_2 = Type(2.0) / (Type(1.0) + exp(-Type(2.0) * transf_rho_2)) - Type(1.0);
  
  // SPDE hyper transforms
  Type kappa_1 = sqrt(Type(8.0)) / range_1;
  Type tau_0_1 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_1*sigma_0_1);
  Type tau_t_1 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_1*sigma_t_1);
  Type tau_flag_1 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_1*sigma_flag_1);
  Type kappa_2 = sqrt(Type(8.0)) / range_2;
  Type tau_0_2 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_2*sigma_0_2);
  Type tau_t_2 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_2*sigma_t_2);
  Type tau_flag_2 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_2*sigma_flag_2);
  
  // PC‐priors
  int share_range_1 = 0;
  if (use_pop_spatial == 1) {
    nll_prior -= pc_prior_matern(
      range_1, sigma_0_1,
      matern_range, matern_sigma_0,
      range_prob, sigma_prob,
      1,
      share_range_1
    );
    share_range_1 = 1;
  }
  if (use_pop_spatiotemporal == 1) {
    nll_prior -= pc_prior_matern(
      range_1, sigma_t_1,
      matern_range, matern_sigma_t,
      range_prob, sigma_prob,
      1,
      share_range_1
    );
    share_range_1 = 1;
  }
  if (use_q_diffs_spatial == 1) {
    nll_prior -= pc_prior_matern(
      range_1, sigma_flag_1,
      matern_range, matern_sigma_flag,
      range_prob, sigma_prob,
      1,
      share_range_1
    );
    share_range_1 = 1;
  }

  int share_range_2 = 0;
  if (use_pop_spatial == 1) {
    nll_prior -= pc_prior_matern(
      range_2, sigma_0_2,
      matern_range, matern_sigma_0,
      range_prob, sigma_prob,
      1,
      share_range_2
    );
    share_range_2 = 1;
  }
  if (use_pop_spatiotemporal == 1) {
    nll_prior -= pc_prior_matern(
      range_2, sigma_t_2,
      matern_range, matern_sigma_t,
      range_prob, sigma_prob,
      1,
      share_range_2
    );
    share_range_2 = 1;
  }
  if (use_q_diffs_spatial == 1) {
    nll_prior -= pc_prior_matern(
      range_2, sigma_flag_2,
      matern_range, matern_sigma_flag,
      range_prob, sigma_prob,
      1,
      share_range_2
    );
    share_range_2 = 1;
  }
  
  
  const bool use_any_spde_1 = (use_pop_spatial == 1 || use_pop_spatiotemporal == 1 || use_q_diffs_spatial == 1);
  const bool use_any_spde_2 = (use_pop_spatial == 1 || use_pop_spatiotemporal == 1 || use_q_diffs_spatial == 1);

  // Need to parameterize H matrix such that det(H)=1 (preserving volume)
  // Note that H appears in (20) in Lindgren et al 2011
  matrix<Type> H(2,2);
  SparseMatrix<Type> Q_1;
  SparseMatrix<Type> Q_2;
  if (use_any_spde_1 || use_any_spde_2) {
    H(0,0) = exp(ln_H_input(0));
    H(1,0) = ln_H_input(1);
    H(0,1) = ln_H_input(1);
    H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
    if (use_any_spde_1) Q_1 = Q_spde(spde, kappa_1, H); // Precision matrix
    if (use_any_spde_2) Q_2 = Q_spde(spde, kappa_2, H); // Precision matrix
  }
  
  
  // Encounter probability
  // Time-constant SPDE (Omega)
  if (use_pop_spatial == 1) {
    nll += SCALE(GMRF(Q_1), 1./tau_0_1)(omega_s_1);
  }

  // Time-varying SPDE (Epsilon)
  if (use_pop_spatiotemporal == 1) {
    if (use_pop_spatiotemporal_rw == 1) {
      nll += SCALE(GMRF(Q_1), 1./tau_t_1)(epsilon_st_1.col(0)); // t=0
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_1), 1./tau_t_1)(epsilon_st_1.col(t) - epsilon_st_1.col(t - 1));
      }
    } else if (use_pop_spatiotemporal_ar1 == 1) {
      nll += SCALE(GMRF(Q_1), 1./(tau_t_1 * sqrt(Type(1.0) - rho_1 * rho_1)))(epsilon_st_1.col(0));
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_1), 1./tau_t_1)(epsilon_st_1.col(t) - rho_1 * epsilon_st_1.col(t - 1));
      }
    } else {
      for (int t = 0; t < n_t; t++) {
        nll += SCALE(GMRF(Q_1), 1./tau_t_1)(epsilon_st_1.col(t));
      }
    }
  }
   
  // Random vessel effect
  if (use_vessel_effect == 1) {
    for(int i=0; i<n_v; i++){
      nll -= dnorm(ves_v_1(i), Type(0.0), ves_std_dev_1, true);
    }
  }
  
  // Differences in gear-speciafic catchability
  if (use_q_diffs_system == 1) {
    for(int i=0; i<n_f-1; i++){
      nll -= dnorm(flag_f_1(i), Type(0.0), flag_std_dev_1, true);
    }
  }

  // Differences in gear-specific catchability over time
  if (use_q_diffs_time == 1) {
    for(int j=0; j<n_f-1; j++){
      for(int t=0; t<n_t; t++){
        if (has_tf(t, j) == 1) {
          nll -= dnorm(flag_t_1(t, j), Type(0.0), flag_t_std_dev_1, true);
        }
      }
    }
  }
  
  // Differences in gear-specific catchability over space
  if (use_q_diffs_spatial == 1) {
    for(int i=0; i<n_f-1; i++){
      nll += SCALE(GMRF(Q_1), 1./tau_flag_1)(flag_s_1.col(i));
    }
  }
  
  
  // Positive catch
  // Time-constant SPDE (Omega)
  if (use_pop_spatial == 1) {
    nll += SCALE(GMRF(Q_2), 1./tau_0_2)(omega_s_2);
  }
  
  // Time-varying SPDE (Epsilon)
  if (use_pop_spatiotemporal == 1) {
    if (use_pop_spatiotemporal_rw == 1) {
      nll += SCALE(GMRF(Q_2), 1./tau_t_2)(epsilon_st_2.col(0)); // t=0
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_2), 1./tau_t_2)(epsilon_st_2.col(t) - epsilon_st_2.col(t - 1));
      }
    } else if (use_pop_spatiotemporal_ar1 == 1) {
      nll += SCALE(GMRF(Q_2), 1./(tau_t_2 * sqrt(Type(1.0) - rho_2 * rho_2)))(epsilon_st_2.col(0));
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_2), 1./tau_t_2)(epsilon_st_2.col(t) - rho_2 * epsilon_st_2.col(t - 1));
      }
    } else {
      for (int t = 0; t < n_t; t++) {
        nll += SCALE(GMRF(Q_2), 1./tau_t_2)(epsilon_st_2.col(t));
      }
    }
  }
   
  // Random vessel effect
  if (use_vessel_effect == 1) {
    for(int i=0; i<n_v; i++){
      nll -= dnorm(ves_v_2(i), Type(0.0), ves_std_dev_2, true);
    }
  }

  // Differences in gear-specific catchability
  if (use_q_diffs_system == 1) {
    for(int i=0; i<n_f-1; i++){
      nll -= dnorm(flag_f_2(i), Type(0.0), flag_std_dev_2, true);
    }
  }
  
  // Differences in gear-specific catchability over time
  if (use_q_diffs_time == 1) {
    for(int j=0; j<n_f-1; j++){
      for(int t=0; t<n_t; t++){
        if (has_tf(t, j) == 1) {
          nll -= dnorm(flag_t_2(t, j), Type(0.0), flag_t_std_dev_2, true);
        }
      }
    }
  }
  
  // Differences in gear-specific catchability over space
  if (use_q_diffs_spatial == 1) {
    for(int i=0; i<n_f-1; i++){
      nll += SCALE(GMRF(Q_2), 1./tau_flag_2)(flag_s_2.col(i));
    }
  }
  
  
  // ---- mean-centering for flag_t over observed times (per flag column) ----
  vector<Type> flag_t_mean_1(n_f - 1);
  vector<Type> flag_t_mean_2(n_f - 1);
  flag_t_mean_1.setZero();
  flag_t_mean_2.setZero();

  if (use_q_diffs_time == 1 && n_f > 1) {
    for (int j = 0; j < n_f - 1; j++) {
      Type sum1 = 0.0, sum2 = 0.0;
      Type cnt  = 0.0;
      for (int t = 0; t < n_t; t++) {
        if (has_tf(t, j) == 1) {
          sum1 += flag_t_1(t, j);
          sum2 += flag_t_2(t, j);
          cnt  += 1.0;
        }
      }
      if (cnt > 0) {
        flag_t_mean_1(j) = sum1 / cnt;
        flag_t_mean_2(j) = sum2 / cnt;
      } else {
        // no data at all for this flag column across time: mean stays 0
        flag_t_mean_1(j) = 0.0;
        flag_t_mean_2(j) = 0.0;
      }
    }
  }
    
	
  // ---- Smooth contributions + priors ----
  vector<Type> eta_smooth_catch_1(n_i);
  vector<Type> eta_smooth_catch_2(n_i);
  vector<Type> eta_smooth_pop_i_1(n_i);
  vector<Type> eta_smooth_pop_i_2(n_i);
  int n_gt = n_g * n_t;
  vector<Type> eta_smooth_pop_g_1(n_gt);
  vector<Type> eta_smooth_pop_g_2(n_gt);
  eta_smooth_catch_1.setZero();
  eta_smooth_catch_2.setZero();
  eta_smooth_pop_i_1.setZero();
  eta_smooth_pop_i_2.setZero();
  eta_smooth_pop_g_1.setZero();
  eta_smooth_pop_g_2.setZero();

  // Observation-level fitted values for diagnostics
  vector<Type> eta_hat_encounter_i(n_i);
  vector<Type> eta_hat_positive_i(n_i);
  vector<Type> encounter_prob_i(n_i);
  vector<Type> log_positive_mean_i(n_i);
  vector<Type> mu_hat_i(n_i);
  eta_hat_encounter_i.setZero();
  eta_hat_positive_i.setZero();
  encounter_prob_i.setZero();
  log_positive_mean_i.setZero();
  mu_hat_i.setZero();

  if (has_smooths_catch) {
    int n_smooth = b_smooth_start_catch.size();

    for (int s = 0; s < n_smooth; s++) {
      int k_s = Zs_catch(s).cols();
      int start = b_smooth_start_catch(s);
      Type smooth_sd0 = exp(ln_smooth_sigma_catch(s, 0));
      Type smooth_sd1 = exp(ln_smooth_sigma_catch(s, 1));

      vector<Type> beta0(k_s);
      for (int j = 0; j < k_s; j++) {
        beta0(j) = b_smooth_catch(start + j, 0);
        nll -= dnorm(beta0(j), Type(0.0), smooth_sd0, true);
      }
      eta_smooth_catch_1 += Zs_catch(s) * beta0;

      vector<Type> beta1(k_s);
      for (int j = 0; j < k_s; j++) {
        beta1(j) = b_smooth_catch(start + j, 1);
        nll -= dnorm(beta1(j), Type(0.0), smooth_sd1, true);
      }
      eta_smooth_catch_2 += Zs_catch(s) * beta1;
    }

    if (Xs_catch.cols() > 0) {
      eta_smooth_catch_1 += Xs_catch * vector<Type>(bs_catch.col(0));
      eta_smooth_catch_2 += Xs_catch * vector<Type>(bs_catch.col(1));
    }
  }

  if (has_smooths_pop) {
    int n_smooth = b_smooth_start_pop.size();

    for (int s = 0; s < n_smooth; s++) {
      int k_s = Zs_pop_i(s).cols();
      int start = b_smooth_start_pop(s);
      Type smooth_sd0 = exp(ln_smooth_sigma_pop(s, 0));
      Type smooth_sd1 = exp(ln_smooth_sigma_pop(s, 1));

      vector<Type> beta0(k_s);
      for (int j = 0; j < k_s; j++) {
        beta0(j) = b_smooth_pop(start + j, 0);
        nll -= dnorm(beta0(j), Type(0.0), smooth_sd0, true);
      }
      eta_smooth_pop_i_1 += Zs_pop_i(s) * beta0;
      eta_smooth_pop_g_1 += Zs_pop_g(s) * beta0;

      vector<Type> beta1(k_s);
      for (int j = 0; j < k_s; j++) {
        beta1(j) = b_smooth_pop(start + j, 1);
        nll -= dnorm(beta1(j), Type(0.0), smooth_sd1, true);
      }
      eta_smooth_pop_i_2 += Zs_pop_i(s) * beta1;
      eta_smooth_pop_g_2 += Zs_pop_g(s) * beta1;
    }

    if (Xs_pop_i.cols() > 0) {
      vector<Type> bs0 = vector<Type>(bs_pop.col(0));
      vector<Type> bs1 = vector<Type>(bs_pop.col(1));
      eta_smooth_pop_i_1 += Xs_pop_i * bs0;
      eta_smooth_pop_i_2 += Xs_pop_i * bs1;
      eta_smooth_pop_g_1 += Xs_pop_g * bs0;
      eta_smooth_pop_g_2 += Xs_pop_g * bs1;
    }
  }

	
  // Response
  vector<Type> s_effect_1(n_i);
  vector<Type> s_effect_2(n_i);
  vector<Type> st_effect_1(n_i);
  vector<Type> st_effect_2(n_i);
  vector<Type> flag_s_effect_1(n_i);
  vector<Type> flag_s_effect_2(n_i);

  s_effect_1.setZero();
  s_effect_2.setZero();
  st_effect_1.setZero();
  st_effect_2.setZero();
  flag_s_effect_1.setZero();
  flag_s_effect_2.setZero();

  if (use_pop_spatial == 1) {
    s_effect_1 = A_is * omega_s_1;
    s_effect_2 = A_is * omega_s_2;
  }
  
  if (use_pop_spatiotemporal == 1 && use_q_diffs_spatial == 1) {
    for(int r=0; r<Ais_ij.rows(); r++){
      int i = Ais_ij(r, 0);  // observation
      int s = Ais_ij(r, 1);  // knot
      int t_id = t_i(i);
      int f_id = f_i(i);
      Type a_is = Ais_x(r);

      st_effect_1(i) += a_is * epsilon_st_1(s, t_id);
      st_effect_2(i) += a_is * epsilon_st_2(s, t_id);

      if (f_id > 0) {
        flag_s_effect_1(i) += a_is * flag_s_1(s, f_id-1);
        flag_s_effect_2(i) += a_is * flag_s_2(s, f_id-1);
      }
    }
  } else if (use_pop_spatiotemporal == 1) {
    for(int r=0; r<Ais_ij.rows(); r++){
      int i = Ais_ij(r, 0);  // observation
      int s = Ais_ij(r, 1);  // knot
      int t_id = t_i(i);
      Type a_is = Ais_x(r);
      st_effect_1(i) += a_is * epsilon_st_1(s, t_id);
      st_effect_2(i) += a_is * epsilon_st_2(s, t_id);
    }
  } else if (use_q_diffs_spatial == 1) {
    for(int r=0; r<Ais_ij.rows(); r++){
      int i = Ais_ij(r, 0);  // observation
      int s = Ais_ij(r, 1);  // knot
      int f_id = f_i(i);
      if (f_id > 0) {
        Type a_is = Ais_x(r);
        flag_s_effect_1(i) += a_is * flag_s_1(s, f_id-1);
        flag_s_effect_2(i) += a_is * flag_s_2(s, f_id-1);
      }
    }
  }
  
  for(int i=0; i<n_i; i++){	
	int tid = t_i(i);
	int vid = v_i(i);
	int fid = f_i(i);
	
    Type ves_effect_1 = Type(0.0);
    if (use_vessel_effect == 1) ves_effect_1 = ves_v_1(vid);
	Type yq_effect_1 = yq_t_1(tid);
	Type flag_effect_1 = Type(0.0);
    if (use_q_diffs_system == 1 && fid > 0) flag_effect_1 = flag_f_1(fid-1);
	
    Type ves_effect_2 = Type(0.0);
    if (use_vessel_effect == 1) ves_effect_2 = ves_v_2(vid);
	Type yq_effect_2 = yq_t_2(tid);
	Type flag_effect_2 = Type(0.0);
    if (use_q_diffs_system == 1 && fid > 0) flag_effect_2 = flag_f_2(fid-1);
	
	Type flag_t_effect_1 = Type(0);
    Type flag_t_effect_2 = Type(0);

    if (use_q_diffs_time == 1 && fid > 0) {
      int j = fid - 1; // column index in flag_t
      // only apply if (t,flag) is observed; otherwise force 0
      if (has_tf(tid, j) == 1) {
        flag_t_effect_1 = flag_t_1(tid, j) - flag_t_mean_1(j);
        flag_t_effect_2 = flag_t_2(tid, j) - flag_t_mean_2(j);
      } else {
        flag_t_effect_1 = 0.0;
        flag_t_effect_2 = 0.0;
      }
    }

	Type eta1 = ves_effect_1 + yq_effect_1 + flag_effect_1 + flag_t_effect_1 + flag_s_effect_1(i) +
      s_effect_1(i) + st_effect_1(i) + eta_smooth_catch_1(i) + eta_smooth_pop_i_1(i);
    Type eta2 = ves_effect_2 + yq_effect_2 + flag_effect_2 + flag_t_effect_2 + flag_s_effect_2(i) +
      s_effect_2(i) + st_effect_2(i) + eta_smooth_catch_2(i) + eta_smooth_pop_i_2(i);
	  
    // Poisson-link; i.e., complementary log–log (cloglog)
    Type log_one_minus_p = -1.0*exp(eta1);
    Type logp = logspace_sub( Type(0.0), log_one_minus_p );

    // Positive
	Type logcat = eta1 + eta2 - logp;

    eta_hat_encounter_i(i) = eta1;
    eta_hat_positive_i(i) = eta2;
    encounter_prob_i(i) = exp(logp);
    log_positive_mean_i(i) = logcat;
    mu_hat_i(i) = exp(eta1 + eta2);
   
	//Likelihood
	// Probability
	nll -= e_i(i) * logp + (1 - e_i(i)) * log_one_minus_p;
				 
	//Positive catch
	if(e_i(i) == 1){
      Type sd_i = sd;
      if (use_flag_sd == 1) sd_i = sd_flag(fid);
	  nll -= dlnorm_bc( b_i(i), logcat, sd_i, true ); 
	}
  }
  
  
  //Projection
  vector<Type> s_effect_proj_1(n_g);
  vector<Type> s_effect_proj_2(n_g);
  matrix<Type> st_effect_proj_1(n_g, n_t);
  matrix<Type> st_effect_proj_2(n_g, n_t);
  s_effect_proj_1.setZero();
  s_effect_proj_2.setZero();
  st_effect_proj_1.setZero();
  st_effect_proj_2.setZero();

  if (use_pop_spatial == 1) {
    s_effect_proj_1 = A_gs * omega_s_1;
    s_effect_proj_2 = A_gs * omega_s_2;
  }
  if (use_pop_spatiotemporal == 1) {
    st_effect_proj_1 = A_gs * epsilon_st_1; // [g, t]
    st_effect_proj_2 = A_gs * epsilon_st_2;
  }
  
  matrix<Type> cpue_density(n_g, n_t);
  vector<Type> mu_total(n_t); // Area-weighted CPUE in weight
  vector<Type> link_total(n_t);
  
  mu_total.setZero();

  for(int t=0; t<n_t; t++){ 
	Type yq_effect_proj_1 = yq_t_1(t); 
	Type yq_effect_proj_2 = yq_t_2(t); 
	
    for(int g=0; g<n_g; g++){	
      int gt = g + n_g * t;
	  Type eta1_proj = yq_effect_proj_1 + s_effect_proj_1(g) + st_effect_proj_1(g, t) + eta_smooth_pop_g_1(gt);
	  Type eta2_proj = yq_effect_proj_2 + s_effect_proj_2(g) + st_effect_proj_2(g, t) + eta_smooth_pop_g_2(gt);
  	      
	  // CPUE in weight
	  Type cpue = exp(eta1_proj + eta2_proj);
	  cpue_density(g, t) = cpue;
	  
	  Type adds = cpue * area_g(g);
		
	  // Area-weighted CPUE
      mu_total(t) += adds;
    }
	
	link_total(t) = log(mu_total(t));
  }
  
  
  // --- Bias correction "epsilon trick" ---
  PARAMETER_VECTOR(eps_index); // length 0 normally; length n_t for bias correction
  
  if (eps_index.size() > 0) {
    Type S;
    for (int t=0; t < n_t; t++) {
      S = mu_total(t);
      nll_penalty += eps_index(t) * S;
    }
  }
  
  nll += nll_prior;   
  nll += nll_penalty; 
  
  // Reporting;
  REPORT(nll_prior);
  REPORT(nll_penalty);
  REPORT(use_vessel_effect);
  REPORT(use_pop_spatial);
  REPORT(use_pop_spatiotemporal);
  REPORT(use_pop_spatiotemporal_rw);
  REPORT(use_pop_spatiotemporal_ar1);
  REPORT(use_q_diffs_system);
  REPORT(use_q_diffs_time);
  REPORT(use_q_diffs_spatial);
  REPORT(use_flag_sd);
  REPORT(sd_flag);
  REPORT(eta_hat_encounter_i);
  REPORT(eta_hat_positive_i);
  REPORT(encounter_prob_i);
  REPORT(log_positive_mean_i);
  REPORT(mu_hat_i);
  REPORT(cpue_density);
  ADREPORT(link_total);
  return nll;
}
