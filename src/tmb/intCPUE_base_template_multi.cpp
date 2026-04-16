#define EIGEN_DONT_PARALLELIZE
#include <TMB.hpp>
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef INTCPUE_ENABLE_Q_TIME
#define INTCPUE_ENABLE_Q_TIME 1
#endif

#ifndef INTCPUE_ENABLE_Q_SPATIAL
#define INTCPUE_ENABLE_Q_SPATIAL 1
#endif

using namespace Eigen;
using namespace tmbutils;

template<class Type>
struct LOSM_t : vector<SparseMatrix<Type> > {
  LOSM_t(SEXP x) {
    this->resize(LENGTH(x));
    for (int i = 0; i < LENGTH(x); i++) {
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asSparseMatrix<Type>(sm);
    }
  }
};

template<class Type>
struct SPDEList_t {
  std::vector<R_inla::spde_aniso_t<Type> > x;
  SPDEList_t() {}
  SPDEList_t(SEXP y) {
    int n = LENGTH(y);
    x.reserve(n);
    for (int i = 0; i < n; i++) {
      x.push_back(R_inla::spde_aniso_t<Type>(VECTOR_ELT(y, i)));
    }
  }
  int size() const { return x.size(); }
  const R_inla::spde_aniso_t<Type>& operator()(int i) const { return x[i]; }
  const R_inla::spde_aniso_t<Type>& operator[](int i) const { return x[i]; }
};

template<class Type>
Type dlnorm_bc(const Type& x, const Type& meanlog, const Type& sdlog, int give_log = 0) {
  Type sd2 = sdlog * sdlog;
  Type adjusted_meanlog = meanlog - sd2 / Type(2.0);
  Type logres = dnorm(log(x), adjusted_meanlog, sdlog, true) - log(x);
  if (give_log) return logres;
  return exp(logres);
}

template <class Type>
Type pc_prior_matern(
  Type range,
  Type sigma,
  Type matern_range,
  Type matern_sigma,
  Type range_prob,
  Type sigma_prob,
  int give_log = 0,
  int share_range = 0
) {
  Type d = Type(2.0);
  Type dhalf = d / Type(2.0);
  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(sigma_prob) / matern_sigma;
  Type range_ll = log(dhalf) +
    log(lam1) +
    (-Type(1.0) - dhalf) * log(range) -
    lam1 * pow(range, -dhalf);
  Type sigma_ll = log(lam2) - lam2 * sigma;
  Type penalty = sigma_ll;
  if (!share_range) penalty += range_ll;
  if (give_log) return penalty;
  return exp(penalty);
}

template<class Type>
matrix<Type> make_H(Type h1, Type h2) {
  matrix<Type> H(2, 2);
  H(0, 0) = exp(h1);
  H(1, 0) = h2;
  H(0, 1) = h2;
  H(1, 1) = (Type(1.0) + h2 * h2) / exp(h1);
  return H;
}

template<class Type>
Type objective_function<Type>::operator() () {
  using namespace density;
  using namespace R_inla;

  DATA_INTEGER(n_a);
  DATA_INTEGER(n_i);
  DATA_INTEGER(n_t);
  DATA_INTEGER(n_v);
  DATA_INTEGER(n_f);
  DATA_INTEGER(n_g);
  DATA_IVECTOR(n_i_area);
  DATA_IVECTOR(n_g_area);
  DATA_IVECTOR(n_s_area);
  DATA_INTEGER(use_pop_spatiotemporal_rw);
  DATA_INTEGER(use_pop_spatiotemporal_ar1);
  DATA_INTEGER(use_q_diffs_time);
  DATA_INTEGER(use_q_diffs_spatial);
  DATA_INTEGER(use_flag_sd);

  DATA_VECTOR(b_i);
  DATA_IVECTOR(e_i);
  DATA_IVECTOR(t_i);
  DATA_IVECTOR(v_i);
  DATA_IVECTOR(f_i);
  DATA_IVECTOR(area_i);
  DATA_IMATRIX(flag_t_index);
  DATA_VECTOR(area_g);

  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);
#if INTCPUE_ENABLE_Q_SPATIAL
  DATA_SPARSE_MATRIX(A_flag_is);
#endif
  DATA_IMATRIX(Ais_ij);
  DATA_VECTOR(Ais_x);

  DATA_VECTOR(matern_range);
#if INTCPUE_ENABLE_Q_SPATIAL
  DATA_SCALAR(matern_range_flag);
#endif
  DATA_SCALAR(range_prob);
  DATA_VECTOR(matern_sigma_0);
  DATA_VECTOR(matern_sigma_t);
#if INTCPUE_ENABLE_Q_SPATIAL
  DATA_SCALAR(matern_sigma_flag);
#endif
  DATA_SCALAR(sigma_prob);

  DATA_INTEGER(has_smooths_catch);
  DATA_MATRIX(Xs_catch);
  DATA_STRUCT(Zs_catch, LOSM_t);
  DATA_IVECTOR(b_smooth_start_catch);

  DATA_INTEGER(has_smooths_pop);
  DATA_MATRIX(Xs_pop_i);
  DATA_STRUCT(Zs_pop_i, LOSM_t);
  DATA_MATRIX(Xs_pop_g);
  DATA_STRUCT(Zs_pop_g, LOSM_t);
  DATA_IVECTOR(b_smooth_start_pop);

  DATA_STRUCT(spdes, SPDEList_t);
#if INTCPUE_ENABLE_Q_SPATIAL
  DATA_STRUCT(flag_spde, SPDEList_t);
#endif

  PARAMETER_VECTOR(ln_sd);
  PARAMETER_MATRIX(ln_sd_flag);
  PARAMETER_MATRIX(ln_H_input);
  PARAMETER_VECTOR(ln_range_1);
  PARAMETER_VECTOR(ln_sigma_0_1);
  PARAMETER_VECTOR(ln_sigma_t_1);
  PARAMETER_VECTOR(transf_rho_1);
  PARAMETER_VECTOR(ln_range_2);
  PARAMETER_VECTOR(ln_sigma_0_2);
  PARAMETER_VECTOR(ln_sigma_t_2);
  PARAMETER_VECTOR(transf_rho_2);
  PARAMETER_VECTOR(ves_v_1);
  PARAMETER_VECTOR(ves_v_2);
  PARAMETER(ves_ln_std_dev_1);
  PARAMETER(ves_ln_std_dev_2);
  PARAMETER_MATRIX(yq_t_1);
  PARAMETER_MATRIX(yq_t_2);
  PARAMETER_VECTOR(omega_s_1);
  PARAMETER_VECTOR(omega_s_2);
  PARAMETER_MATRIX(epsilon_st_1);
  PARAMETER_MATRIX(epsilon_st_2);
  PARAMETER_VECTOR(flag_f_1);
  PARAMETER_VECTOR(flag_f_2);
  PARAMETER(flag_ln_std_dev_1);
  PARAMETER(flag_ln_std_dev_2);
#if INTCPUE_ENABLE_Q_TIME
  PARAMETER_VECTOR(flag_t_1);
  PARAMETER_VECTOR(flag_t_2);
  PARAMETER(flag_t_ln_std_dev_1);
  PARAMETER(flag_t_ln_std_dev_2);
#endif
#if INTCPUE_ENABLE_Q_SPATIAL
  PARAMETER_MATRIX(ln_H_flag_input);
  PARAMETER_MATRIX(flag_s_1);
  PARAMETER_MATRIX(flag_s_2);
  PARAMETER(ln_range_flag_1);
  PARAMETER(ln_range_flag_2);
  PARAMETER(ln_sigma_flag_1);
  PARAMETER(ln_sigma_flag_2);
#endif
  PARAMETER_MATRIX(bs_catch);
  PARAMETER_MATRIX(b_smooth_catch);
  PARAMETER_MATRIX(ln_smooth_sigma_catch);
  PARAMETER_MATRIX(bs_pop);
  PARAMETER_MATRIX(b_smooth_pop);
  PARAMETER_MATRIX(ln_smooth_sigma_pop);
  PARAMETER_VECTOR(eps_index);

  Type nll = 0;
  Type nll_prior = 0;
  Type nll_penalty = 0;
  Type ves_std_dev_1 = exp(ves_ln_std_dev_1);
  Type ves_std_dev_2 = exp(ves_ln_std_dev_2);
  Type flag_std_dev_1 = exp(flag_ln_std_dev_1);
  Type flag_std_dev_2 = exp(flag_ln_std_dev_2);
#if INTCPUE_ENABLE_Q_TIME
  Type flag_t_std_dev_1 = exp(flag_t_ln_std_dev_1);
  Type flag_t_std_dev_2 = exp(flag_t_ln_std_dev_2);
#endif

  matrix<Type> sd_flag(n_f, n_a);
  for (int a = 0; a < n_a; a++) {
    for (int f = 0; f < n_f; f++) {
      sd_flag(f, a) = exp(ln_sd_flag(f, a));
    }
  }

  vector<Type> eta_smooth_catch_1(n_i);
  vector<Type> eta_smooth_catch_2(n_i);
  vector<Type> eta_smooth_pop_i_1(n_i);
  vector<Type> eta_smooth_pop_i_2(n_i);
  vector<Type> eta_smooth_pop_g_1(n_g * n_t);
  vector<Type> eta_smooth_pop_g_2(n_g * n_t);
  eta_smooth_catch_1.setZero();
  eta_smooth_catch_2.setZero();
  eta_smooth_pop_i_1.setZero();
  eta_smooth_pop_i_2.setZero();
  eta_smooth_pop_g_1.setZero();
  eta_smooth_pop_g_2.setZero();

  if (has_smooths_catch) {
    int n_smooth = b_smooth_start_catch.size();
    for (int s = 0; s < n_smooth; s++) {
      int k_s = Zs_catch(s).cols();
      int start = b_smooth_start_catch(s);
      Type smooth_sd0 = exp(ln_smooth_sigma_catch(s, 0));
      Type smooth_sd1 = exp(ln_smooth_sigma_catch(s, 1));
      vector<Type> beta0(k_s);
      vector<Type> beta1(k_s);
      for (int j = 0; j < k_s; j++) {
        beta0(j) = b_smooth_catch(start + j, 0);
        beta1(j) = b_smooth_catch(start + j, 1);
        nll -= dnorm(beta0(j), Type(0.0), smooth_sd0, true);
        nll -= dnorm(beta1(j), Type(0.0), smooth_sd1, true);
      }
      eta_smooth_catch_1 += Zs_catch(s) * beta0;
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
      vector<Type> beta1(k_s);
      for (int j = 0; j < k_s; j++) {
        beta0(j) = b_smooth_pop(start + j, 0);
        beta1(j) = b_smooth_pop(start + j, 1);
        nll -= dnorm(beta0(j), Type(0.0), smooth_sd0, true);
        nll -= dnorm(beta1(j), Type(0.0), smooth_sd1, true);
      }
      eta_smooth_pop_i_1 += Zs_pop_i(s) * beta0;
      eta_smooth_pop_i_2 += Zs_pop_i(s) * beta1;
      eta_smooth_pop_g_1 += Zs_pop_g(s) * beta0;
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

  vector<Type> s_effect_proj_1(n_g);
  vector<Type> s_effect_proj_2(n_g);
  matrix<Type> st_effect_proj_1(n_g, n_t);
  matrix<Type> st_effect_proj_2(n_g, n_t);
  s_effect_proj_1.setZero();
  s_effect_proj_2.setZero();
  st_effect_proj_1.setZero();
  st_effect_proj_2.setZero();

#if INTCPUE_ENABLE_Q_TIME
  int n_flag_t_cols = n_f - 1;
  matrix<Type> flag_t_full_1(n_t, n_flag_t_cols);
  matrix<Type> flag_t_full_2(n_t, n_flag_t_cols);
  flag_t_full_1.setZero();
  flag_t_full_2.setZero();
#endif

  int s_offset = 0;
  for (int a = 0; a < n_a; a++) {
    int ns = n_s_area(a);

    Type range_1 = exp(ln_range_1(a));
    Type sigma_0_1 = exp(ln_sigma_0_1(a));
    Type sigma_t_1 = exp(ln_sigma_t_1(a));
    Type rho_1 = Type(2.0) / (Type(1.0) + exp(-Type(2.0) * transf_rho_1(a))) - Type(1.0);
    Type kappa_1 = sqrt(Type(8.0)) / range_1;

    Type range_2 = exp(ln_range_2(a));
    Type sigma_0_2 = exp(ln_sigma_0_2(a));
    Type sigma_t_2 = exp(ln_sigma_t_2(a));
    Type rho_2 = Type(2.0) / (Type(1.0) + exp(-Type(2.0) * transf_rho_2(a))) - Type(1.0);
    Type kappa_2 = sqrt(Type(8.0)) / range_2;

    matrix<Type> H = make_H(ln_H_input(a, 0), ln_H_input(a, 1));
    SparseMatrix<Type> Q_1 = Q_spde(spdes(a), kappa_1, H);
    SparseMatrix<Type> Q_2 = Q_spde(spdes(a), kappa_2, H);

    int share_range_1 = 0;
    int share_range_2 = 0;

    nll_prior -= pc_prior_matern(range_1, sigma_0_1, matern_range(a), matern_sigma_0(a), range_prob, sigma_prob, 1, share_range_1);
    share_range_1 = 1;
    Type tau_0_1 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_1 * sigma_0_1);
    nll += SCALE(GMRF(Q_1), Type(1.0) / tau_0_1)(omega_s_1.segment(s_offset, ns));

    nll_prior -= pc_prior_matern(range_2, sigma_0_2, matern_range(a), matern_sigma_0(a), range_prob, sigma_prob, 1, share_range_2);
    share_range_2 = 1;
    Type tau_0_2 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_2 * sigma_0_2);
    nll += SCALE(GMRF(Q_2), Type(1.0) / tau_0_2)(omega_s_2.segment(s_offset, ns));

    nll_prior -= pc_prior_matern(range_1, sigma_t_1, matern_range(a), matern_sigma_t(a), range_prob, sigma_prob, 1, share_range_1);
    share_range_1 = 1;
    Type tau_t_1 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_1 * sigma_t_1);
    if (use_pop_spatiotemporal_rw == 1) {
      nll += SCALE(GMRF(Q_1), Type(1.0) / tau_t_1)(epsilon_st_1.col(0).segment(s_offset, ns));
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_1), Type(1.0) / tau_t_1)(
          epsilon_st_1.col(t).segment(s_offset, ns) - epsilon_st_1.col(t - 1).segment(s_offset, ns)
        );
      }
    } else if (use_pop_spatiotemporal_ar1 == 1) {
      nll += SCALE(GMRF(Q_1), Type(1.0) / (tau_t_1 * sqrt(Type(1.0) - rho_1 * rho_1)))(epsilon_st_1.col(0).segment(s_offset, ns));
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_1), Type(1.0) / tau_t_1)(
          epsilon_st_1.col(t).segment(s_offset, ns) - rho_1 * epsilon_st_1.col(t - 1).segment(s_offset, ns)
        );
      }
    }

    nll_prior -= pc_prior_matern(range_2, sigma_t_2, matern_range(a), matern_sigma_t(a), range_prob, sigma_prob, 1, share_range_2);
    share_range_2 = 1;
    Type tau_t_2 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_2 * sigma_t_2);
    if (use_pop_spatiotemporal_rw == 1) {
      nll += SCALE(GMRF(Q_2), Type(1.0) / tau_t_2)(epsilon_st_2.col(0).segment(s_offset, ns));
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_2), Type(1.0) / tau_t_2)(
          epsilon_st_2.col(t).segment(s_offset, ns) - epsilon_st_2.col(t - 1).segment(s_offset, ns)
        );
      }
    } else if (use_pop_spatiotemporal_ar1 == 1) {
      nll += SCALE(GMRF(Q_2), Type(1.0) / (tau_t_2 * sqrt(Type(1.0) - rho_2 * rho_2)))(epsilon_st_2.col(0).segment(s_offset, ns));
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_2), Type(1.0) / tau_t_2)(
          epsilon_st_2.col(t).segment(s_offset, ns) - rho_2 * epsilon_st_2.col(t - 1).segment(s_offset, ns)
        );
      }
    }

    s_offset += ns;
  }

  for (int j = 0; j < n_f - 1; j++) {
    nll -= dnorm(flag_f_1(j), Type(0.0), flag_std_dev_1, true);
    nll -= dnorm(flag_f_2(j), Type(0.0), flag_std_dev_2, true);
  }

#if INTCPUE_ENABLE_Q_TIME
  if (use_q_diffs_time == 1) {
    for (int j = 0; j < n_f - 1; j++) {
      for (int t = 0; t < n_t; t++) {
        int idx = flag_t_index(t, j);
        if (idx >= 0) {
          flag_t_full_1(t, j) = flag_t_1(idx);
          flag_t_full_2(t, j) = flag_t_2(idx);
          nll -= dnorm(flag_t_1(idx), Type(0.0), flag_t_std_dev_1, true);
          nll -= dnorm(flag_t_2(idx), Type(0.0), flag_t_std_dev_2, true);
        }
      }
    }
  }
#endif

#if INTCPUE_ENABLE_Q_SPATIAL
  if (use_q_diffs_spatial == 1) {
    Type range_flag_1 = exp(ln_range_flag_1);
    Type range_flag_2 = exp(ln_range_flag_2);
    Type sigma_flag_1 = exp(ln_sigma_flag_1);
    Type sigma_flag_2 = exp(ln_sigma_flag_2);
    Type kappa_flag_1 = sqrt(Type(8.0)) / range_flag_1;
    Type kappa_flag_2 = sqrt(Type(8.0)) / range_flag_2;
    matrix<Type> H_flag = make_H(ln_H_flag_input(0, 0), ln_H_flag_input(0, 1));
    SparseMatrix<Type> Q_flag_1 = Q_spde(flag_spde(0), kappa_flag_1, H_flag);
    SparseMatrix<Type> Q_flag_2 = Q_spde(flag_spde(0), kappa_flag_2, H_flag);
    nll_prior -= pc_prior_matern(range_flag_1, sigma_flag_1, matern_range_flag, matern_sigma_flag, range_prob, sigma_prob, 1, 0);
    nll_prior -= pc_prior_matern(range_flag_2, sigma_flag_2, matern_range_flag, matern_sigma_flag, range_prob, sigma_prob, 1, 0);
    Type tau_flag_1 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_flag_1 * sigma_flag_1);
    Type tau_flag_2 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_flag_2 * sigma_flag_2);
    for (int j = 0; j < n_f - 1; j++) {
      nll += SCALE(GMRF(Q_flag_1), Type(1.0) / tau_flag_1)(flag_s_1.col(j));
      nll += SCALE(GMRF(Q_flag_2), Type(1.0) / tau_flag_2)(flag_s_2.col(j));
    }
  }
#endif

  for (int i = 0; i < n_v; i++) {
    nll -= dnorm(ves_v_1(i), Type(0.0), ves_std_dev_1, true);
    nll -= dnorm(ves_v_2(i), Type(0.0), ves_std_dev_2, true);
  }

  s_effect_1 = A_is * omega_s_1;
  s_effect_2 = A_is * omega_s_2;
  s_effect_proj_1 = A_gs * omega_s_1;
  s_effect_proj_2 = A_gs * omega_s_2;
  st_effect_proj_1 = A_gs * epsilon_st_1;
  st_effect_proj_2 = A_gs * epsilon_st_2;

#if INTCPUE_ENABLE_Q_SPATIAL
  matrix<Type> flag_s_linpred_1(n_i, n_f - 1);
  matrix<Type> flag_s_linpred_2(n_i, n_f - 1);
  flag_s_linpred_1.setZero();
  flag_s_linpred_2.setZero();
  if (use_q_diffs_spatial == 1 && n_f > 1) {
    flag_s_linpred_1 = A_flag_is * flag_s_1;
    flag_s_linpred_2 = A_flag_is * flag_s_2;
  }
#endif

  for (int r = 0; r < Ais_ij.rows(); r++) {
    int i = Ais_ij(r, 0);
    int s = Ais_ij(r, 1);
    int t_id = t_i(i);
    Type a_is = Ais_x(r);
    st_effect_1(i) += a_is * epsilon_st_1(s, t_id);
    st_effect_2(i) += a_is * epsilon_st_2(s, t_id);
  }

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

  for (int i = 0; i < n_i; i++) {
    int tid = t_i(i);
    int vid = v_i(i);
    int fid = f_i(i);
    int aid = area_i(i);

    Type ves_effect_1 = Type(0.0);
    Type ves_effect_2 = Type(0.0);
    if (n_v > 0) {
      ves_effect_1 = ves_v_1(vid);
      ves_effect_2 = ves_v_2(vid);
    }

    Type yq_effect_1 = yq_t_1(tid, aid);
    Type yq_effect_2 = yq_t_2(tid, aid);
    Type flag_effect_1 = Type(0.0);
    Type flag_effect_2 = Type(0.0);
    if (fid > 0 && n_f > 1) {
      flag_effect_1 = flag_f_1(fid - 1);
      flag_effect_2 = flag_f_2(fid - 1);
    }

    Type flag_t_effect_1 = Type(0.0);
    Type flag_t_effect_2 = Type(0.0);
#if INTCPUE_ENABLE_Q_TIME
    if (use_q_diffs_time == 1 && fid > 0 && n_f > 1) {
      int col = fid - 1;
      flag_t_effect_1 = flag_t_full_1(tid, col);
      flag_t_effect_2 = flag_t_full_2(tid, col);
    }
#endif
#if INTCPUE_ENABLE_Q_SPATIAL
    if (use_q_diffs_spatial == 1 && fid > 0 && n_f > 1) {
      int col = fid - 1;
      flag_s_effect_1(i) = flag_s_linpred_1(i, col);
      flag_s_effect_2(i) = flag_s_linpred_2(i, col);
    }
#endif

    Type eta1 = ves_effect_1 + yq_effect_1 + flag_effect_1 + flag_t_effect_1 + flag_s_effect_1(i) +
      s_effect_1(i) + st_effect_1(i) + eta_smooth_catch_1(i) + eta_smooth_pop_i_1(i);
    Type eta2 = ves_effect_2 + yq_effect_2 + flag_effect_2 + flag_t_effect_2 + flag_s_effect_2(i) +
      s_effect_2(i) + st_effect_2(i) + eta_smooth_catch_2(i) + eta_smooth_pop_i_2(i);

    Type log_one_minus_p = -exp(eta1);
    Type logp = logspace_sub(Type(0.0), log_one_minus_p);
    Type logcat = eta1 + eta2 - logp;

    eta_hat_encounter_i(i) = eta1;
    eta_hat_positive_i(i) = eta2;
    encounter_prob_i(i) = exp(logp);
    log_positive_mean_i(i) = logcat;
    mu_hat_i(i) = exp(eta1 + eta2);

    nll -= e_i(i) * logp + (Type(1.0) - e_i(i)) * log_one_minus_p;
    if (e_i(i) == 1) {
      Type sd_i = exp(ln_sd(aid));
      if (use_flag_sd == 1) {
        sd_i = sd_flag(fid, aid);
      }
      nll -= dlnorm_bc(b_i(i), logcat, sd_i, true);
    }
  }

  matrix<Type> cpue_density(n_g, n_t);
  matrix<Type> mu_total(n_t, n_a);
  vector<Type> link_total(n_t * n_a);
  cpue_density.setZero();
  mu_total.setZero();

  int g_start = 0;
  for (int a = 0; a < n_a; a++) {
    int ng = n_g_area(a);
    for (int t = 0; t < n_t; t++) {
      Type yq_effect_proj_1 = yq_t_1(t, a);
      Type yq_effect_proj_2 = yq_t_2(t, a);
      for (int g_local = 0; g_local < ng; g_local++) {
        int g = g_start + g_local;
        int gt = g + n_g * t;
        Type eta1_proj = yq_effect_proj_1 + s_effect_proj_1(g) + st_effect_proj_1(g, t) + eta_smooth_pop_g_1(gt);
        Type eta2_proj = yq_effect_proj_2 + s_effect_proj_2(g) + st_effect_proj_2(g, t) + eta_smooth_pop_g_2(gt);
        Type cpue = exp(eta1_proj + eta2_proj);
        cpue_density(g, t) = cpue;
        mu_total(t, a) += cpue * area_g(g);
      }
      link_total(t + n_t * a) = log(mu_total(t, a));
    }
    g_start += ng;
  }

  if (eps_index.size() > 0) {
    for (int a = 0; a < n_a; a++) {
      for (int t = 0; t < n_t; t++) {
        Type S = newton::Tag(mu_total(t, a));
        nll_penalty += eps_index(t + n_t * a) * S;
      }
    }
  }

  nll += nll_prior;
  nll += nll_penalty;

  REPORT(nll_prior);
  REPORT(nll_penalty);
  REPORT(use_pop_spatiotemporal_rw);
  REPORT(use_pop_spatiotemporal_ar1);
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
