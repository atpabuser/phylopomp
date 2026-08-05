// BDSS: Linear birth-death with superspreading (C++)
#include "master.h"
#include "popul_proc.h"
#include "generics.h"
#include "internal.h"

static const int normal = 1;
static const int superspreader = 2;

//! BDSS process state.
typedef struct {
  int N;
  int S;
} bdss_state_t;

//! BDSS process parameters.
typedef struct {
  double lambda_nn;
  double lambda_ns;
  double lambda_sn;
  double lambda_ss;
  double mu;
  double chi;
  double pop;
  double N0;
  double S0;
} bdss_parameters_t;

using bdss_proc_t = popul_proc_t<bdss_state_t,bdss_parameters_t,8>;
using bdss_genealogy_t = master_t<bdss_proc_t,2>;

template<>
std::string bdss_proc_t::yaml (std::string tab) const {
  std::string t = tab + "  ";
  std::string p = tab + "parameter:\n"
    + YAML_PARAM(lambda_nn)
    + YAML_PARAM(lambda_ns)
    + YAML_PARAM(lambda_sn)
    + YAML_PARAM(lambda_ss)
    + YAML_PARAM(mu)
    + YAML_PARAM(chi)
    + YAML_PARAM(pop)
    + YAML_PARAM(N0)
    + YAML_PARAM(S0);
  std::string s = tab + "state:\n"
    + YAML_STATE(N)
    + YAML_STATE(S);
  return p+s;
}

template<>
void bdss_proc_t::update_params (double *p, int n) {
  int m = 0;
  PARAM_SET(lambda_nn);
  PARAM_SET(lambda_ns);
  PARAM_SET(lambda_sn);
  PARAM_SET(lambda_ss);
  PARAM_SET(mu);
  PARAM_SET(chi);
  if (m != n) err("wrong number of parameters!");
}

template<>
void bdss_proc_t::update_IVPs (double *p, int n) {
  int m = 0;
  PARAM_SET(pop);
  PARAM_SET(N0);
  PARAM_SET(S0);
  if (m != n) err("wrong number of initial-value parameters!");
}

template<>
double bdss_proc_t::event_rates (double *rate, int n) const {
  int m = 0;
  double total = 0;
  RATE_CALC(params.lambda_nn * state.N);
  RATE_CALC(params.lambda_ns * state.N);
  RATE_CALC(params.lambda_sn * state.S);
  RATE_CALC(params.lambda_ss * state.S);
  RATE_CALC(params.mu * state.N);
  RATE_CALC(params.mu * state.S);
  RATE_CALC(params.chi * state.N);
  RATE_CALC(params.chi * state.S);
  if (m != n) err("wrong number of events!");
  return total;
}

template<>
void bdss_genealogy_t::rinit (void) {
  double m = params.pop/(params.N0 + params.S0);
state.N = nearbyint(m*params.N0);
state.S = nearbyint(m*params.S0);
graft(normal, state.N);
graft(superspreader, state.S);
}

template<>
void bdss_genealogy_t::jump (int event) {
  switch (event) {
  case 0:
      state.N += 1; birth(normal, normal);
      break;
    case 1:
      state.S += 1; birth(normal, superspreader);
      break;
    case 2:
      state.N += 1; birth(superspreader, normal);
      break;
    case 3:
      state.S += 1; birth(superspreader, superspreader);
      break;
    case 4:
      state.N -= 1; death(normal);
      break;
    case 5:
      state.S -= 1; death(superspreader);
      break;
    case 6:
      state.N -= 1; sample_death(normal);
      break;
    case 7:
      state.S -= 1; sample_death(superspreader);
      break;
  default:                      // #nocov
    assert(0);                  // #nocov
    break;                      // #nocov
  }
}

GENERICS(BDSS,bdss_genealogy_t)
