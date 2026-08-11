// MTBD2: General 2-type linear multitype birth-death-sampling process (BDMM, Stadler psi/r parameterization) (C++)
#include "master.h"
#include "popul_proc.h"
#include "generics.h"
#include "internal.h"

static const int type1 = 1;
static const int type2 = 2;

//! MTBD2 process state.
typedef struct {
  int I1;
  int I2;
  int n12;
  int n21;
} mtbd2_state_t;

//! MTBD2 process parameters.
typedef struct {
  double lambda_1_1;
  double lambda_1_2;
  double lambda_2_1;
  double lambda_2_2;
  double m_1_2;
  double m_2_1;
  double mu_1;
  double mu_2;
  double psi_1;
  double psi_2;
  double r_1;
  double r_2;
  int I1_0;
  int I2_0;
} mtbd2_parameters_t;

using mtbd2_proc_t = popul_proc_t<mtbd2_state_t,mtbd2_parameters_t,12>;
using mtbd2_genealogy_t = master_t<mtbd2_proc_t,2>;

template<>
std::string mtbd2_proc_t::yaml (std::string tab) const {
  std::string t = tab + "  ";
  std::string p = tab + "parameter:\n"
    + YAML_PARAM(lambda_1_1)
    + YAML_PARAM(lambda_1_2)
    + YAML_PARAM(lambda_2_1)
    + YAML_PARAM(lambda_2_2)
    + YAML_PARAM(m_1_2)
    + YAML_PARAM(m_2_1)
    + YAML_PARAM(mu_1)
    + YAML_PARAM(mu_2)
    + YAML_PARAM(psi_1)
    + YAML_PARAM(psi_2)
    + YAML_PARAM(r_1)
    + YAML_PARAM(r_2)
    + YAML_PARAM(I1_0)
    + YAML_PARAM(I2_0);
  std::string s = tab + "state:\n"
    + YAML_STATE(I1)
    + YAML_STATE(I2)
    + YAML_STATE(n12)
    + YAML_STATE(n21);
  return p+s;
}

template<>
void mtbd2_proc_t::update_params (double *p, int n) {
  int m = 0;
  PARAM_SET(lambda_1_1);
  PARAM_SET(lambda_1_2);
  PARAM_SET(lambda_2_1);
  PARAM_SET(lambda_2_2);
  PARAM_SET(m_1_2);
  PARAM_SET(m_2_1);
  PARAM_SET(mu_1);
  PARAM_SET(mu_2);
  PARAM_SET(psi_1);
  PARAM_SET(psi_2);
  PARAM_SET(r_1);
  PARAM_SET(r_2);
  if (m != n) err("wrong number of parameters!");
}

template<>
void mtbd2_proc_t::update_IVPs (double *p, int n) {
  int m = 0;
  PARAM_SET(I1_0);
  PARAM_SET(I2_0);
  if (m != n) err("wrong number of initial-value parameters!");
}

template<>
double mtbd2_proc_t::event_rates (double *rate, int n) const {
  int m = 0;
  double total = 0;
  RATE_CALC(params.lambda_1_1 * state.I1);
  RATE_CALC(params.lambda_1_2 * state.I1);
  RATE_CALC(params.lambda_2_1 * state.I2);
  RATE_CALC(params.lambda_2_2 * state.I2);
  RATE_CALC(params.m_1_2 * state.I1);
  RATE_CALC(params.m_2_1 * state.I2);
  RATE_CALC(params.mu_1 * state.I1);
  RATE_CALC(params.mu_2 * state.I2);
  RATE_CALC(params.r_1 * params.psi_1 * state.I1);
  RATE_CALC(params.r_2 * params.psi_2 * state.I2);
  RATE_CALC((1 - params.r_1) * params.psi_1 * state.I1);
  RATE_CALC((1 - params.r_2) * params.psi_2 * state.I2);
  if (m != n) err("wrong number of events!");
  return total;
}

template<>
void mtbd2_genealogy_t::rinit (void) {
  state.I1 = params.I1_0;
state.I2 = params.I2_0;
state.n12 = 0;
state.n21 = 0;
graft(type1,params.I1_0);
graft(type2,params.I2_0);
}

template<>
void mtbd2_genealogy_t::jump (int event) {
  switch (event) {
  case 0:
      state.I1 += 1; birth(type1,type1);
      break;
    case 1:
      state.I2 += 1; state.n12 += 1; birth(type1,type2);
      break;
    case 2:
      state.I1 += 1; state.n21 += 1; birth(type2,type1);
      break;
    case 3:
      state.I2 += 1; birth(type2,type2);
      break;
    case 4:
      state.I1 -= 1; state.I2 += 1; migrate(type1,type2);
      break;
    case 5:
      state.I2 -= 1; state.I1 += 1; migrate(type2,type1);
      break;
    case 6:
      state.I1 -= 1; death(type1);
      break;
    case 7:
      state.I2 -= 1; death(type2);
      break;
    case 8:
      state.I1 -= 1; sample_death(type1);
      break;
    case 9:
      state.I2 -= 1; sample_death(type2);
      break;
    case 10:
      sample(type1);
      break;
    case 11:
      sample(type2);
      break;
  default:                      // #nocov
    assert(0);                  // #nocov
    break;                      // #nocov
  }
}

GENERICS(MTBD2,mtbd2_genealogy_t)
