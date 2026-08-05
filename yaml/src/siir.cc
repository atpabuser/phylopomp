// SIIR: Two-strain SIR model (C++)
#include "master.h"
#include "popul_proc.h"
#include "generics.h"
#include "internal.h"

static const int strain1 = 1;
static const int strain2 = 2;

//! SIIR process state.
typedef struct {
  int S;
  int I1;
  int I2;
  int R;
} siir_state_t;

//! SIIR process parameters.
typedef struct {
  double Beta1;
  double Beta2;
  double gamma;
  double psi1;
  double psi2;
  double sigma12;
  double sigma21;
  double omega;
  double pop;
  double S_0;
  double I1_0;
  double I2_0;
  double R_0;
} siir_parameters_t;

using siir_proc_t = popul_proc_t<siir_state_t,siir_parameters_t,9>;
using siir_genealogy_t = master_t<siir_proc_t,2>;

template<>
std::string siir_proc_t::yaml (std::string tab) const {
  std::string t = tab + "  ";
  std::string p = tab + "parameter:\n"
    + YAML_PARAM(Beta1)
    + YAML_PARAM(Beta2)
    + YAML_PARAM(gamma)
    + YAML_PARAM(psi1)
    + YAML_PARAM(psi2)
    + YAML_PARAM(sigma12)
    + YAML_PARAM(sigma21)
    + YAML_PARAM(omega)
    + YAML_PARAM(pop)
    + YAML_PARAM(S_0)
    + YAML_PARAM(I1_0)
    + YAML_PARAM(I2_0)
    + YAML_PARAM(R_0);
  std::string s = tab + "state:\n"
    + YAML_STATE(S)
    + YAML_STATE(I1)
    + YAML_STATE(I2)
    + YAML_STATE(R);
  return p+s;
}

template<>
void siir_proc_t::update_params (double *p, int n) {
  int m = 0;
  PARAM_SET(Beta1);
  PARAM_SET(Beta2);
  PARAM_SET(gamma);
  PARAM_SET(psi1);
  PARAM_SET(psi2);
  PARAM_SET(sigma12);
  PARAM_SET(sigma21);
  PARAM_SET(omega);
  if (m != n) err("wrong number of parameters!");
}

template<>
void siir_proc_t::update_IVPs (double *p, int n) {
  int m = 0;
  PARAM_SET(pop);
  PARAM_SET(S_0);
  PARAM_SET(I1_0);
  PARAM_SET(I2_0);
  PARAM_SET(R_0);
  if (m != n) err("wrong number of initial-value parameters!");
}

template<>
double siir_proc_t::event_rates (double *rate, int n) const {
  int m = 0;
  double total = 0;
  RATE_CALC(params.Beta1 * state.S * state.I1 / params.pop);
  RATE_CALC(params.Beta2 * state.S * state.I2 / params.pop);
  RATE_CALC(params.gamma * state.I1);
  RATE_CALC(params.gamma * state.I2);
  RATE_CALC(params.psi1 * state.I1);
  RATE_CALC(params.psi2 * state.I2);
  RATE_CALC(params.sigma12 * state.I1);
  RATE_CALC(params.sigma21 * state.I2);
  RATE_CALC(params.omega * state.R);
  if (m != n) err("wrong number of events!");
  return total;
}

template<>
void siir_genealogy_t::rinit (void) {
  double f = params.pop/(params.S_0+params.I1_0+params.I2_0+params.R_0);
state.S = nearbyint(f*params.S_0);
state.I1 = nearbyint(f*params.I1_0);
state.I2 = nearbyint(f*params.I2_0);
state.R = nearbyint(f*params.R_0);
graft(strain1,state.I1);
graft(strain2,state.I2);
}

template<>
void siir_genealogy_t::jump (int event) {
  switch (event) {
  case 0:
      state.S -= 1; state.I1 += 1; birth(strain1,strain1);
      break;
    case 1:
      state.S -= 1; state.I2 += 1; birth(strain2,strain2);
      break;
    case 2:
      state.I1 -= 1; state.R += 1; death(strain1);
      break;
    case 3:
      state.I2 -= 1; state.R += 1; death(strain2);
      break;
    case 4:
      sample(strain1);
      break;
    case 5:
      sample(strain2);
      break;
    case 6:
      state.I1 -= 1; state.I2 += 1; migrate(strain1,strain2);
      break;
    case 7:
      state.I1 += 1; state.I2 -= 1; migrate(strain2,strain1);
      break;
    case 8:
      state.S += 1; state.R -= 1;
      break;
  default:                      // #nocov
    assert(0);                  // #nocov
    break;                      // #nocov
  }
}

GENERICS(SIIR,siir_genealogy_t)
