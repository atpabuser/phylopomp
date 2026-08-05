// MTBD3: General 3-type linear multitype birth-death-sampling process (BDMM, Stadler psi/r parameterization) (C++)
#include "master.h"
#include "popul_proc.h"
#include "generics.h"
#include "internal.h"

static const int type1 = 1;
static const int type2 = 2;
static const int type3 = 3;

//! MTBD3 process state.
typedef struct {
  int I1;
  int I2;
  int I3;
} mtbd3_state_t;

//! MTBD3 process parameters.
typedef struct {
  double lambda_1_1;
  double lambda_1_2;
  double lambda_1_3;
  double lambda_2_1;
  double lambda_2_2;
  double lambda_2_3;
  double lambda_3_1;
  double lambda_3_2;
  double lambda_3_3;
  double m_1_2;
  double m_1_3;
  double m_2_1;
  double m_2_3;
  double m_3_1;
  double m_3_2;
  double mu_1;
  double mu_2;
  double mu_3;
  double psi_1;
  double psi_2;
  double psi_3;
  double r_1;
  double r_2;
  double r_3;
  int I1_0;
  int I2_0;
  int I3_0;
} mtbd3_parameters_t;

using mtbd3_proc_t = popul_proc_t<mtbd3_state_t,mtbd3_parameters_t,24>;
using mtbd3_genealogy_t = master_t<mtbd3_proc_t,3>;

template<>
std::string mtbd3_proc_t::yaml (std::string tab) const {
  std::string t = tab + "  ";
  std::string p = tab + "parameter:\n"
    + YAML_PARAM(lambda_1_1)
    + YAML_PARAM(lambda_1_2)
    + YAML_PARAM(lambda_1_3)
    + YAML_PARAM(lambda_2_1)
    + YAML_PARAM(lambda_2_2)
    + YAML_PARAM(lambda_2_3)
    + YAML_PARAM(lambda_3_1)
    + YAML_PARAM(lambda_3_2)
    + YAML_PARAM(lambda_3_3)
    + YAML_PARAM(m_1_2)
    + YAML_PARAM(m_1_3)
    + YAML_PARAM(m_2_1)
    + YAML_PARAM(m_2_3)
    + YAML_PARAM(m_3_1)
    + YAML_PARAM(m_3_2)
    + YAML_PARAM(mu_1)
    + YAML_PARAM(mu_2)
    + YAML_PARAM(mu_3)
    + YAML_PARAM(psi_1)
    + YAML_PARAM(psi_2)
    + YAML_PARAM(psi_3)
    + YAML_PARAM(r_1)
    + YAML_PARAM(r_2)
    + YAML_PARAM(r_3)
    + YAML_PARAM(I1_0)
    + YAML_PARAM(I2_0)
    + YAML_PARAM(I3_0);
  std::string s = tab + "state:\n"
    + YAML_STATE(I1)
    + YAML_STATE(I2)
    + YAML_STATE(I3);
  return p+s;
}

template<>
void mtbd3_proc_t::update_params (double *p, int n) {
  int m = 0;
  PARAM_SET(lambda_1_1);
  PARAM_SET(lambda_1_2);
  PARAM_SET(lambda_1_3);
  PARAM_SET(lambda_2_1);
  PARAM_SET(lambda_2_2);
  PARAM_SET(lambda_2_3);
  PARAM_SET(lambda_3_1);
  PARAM_SET(lambda_3_2);
  PARAM_SET(lambda_3_3);
  PARAM_SET(m_1_2);
  PARAM_SET(m_1_3);
  PARAM_SET(m_2_1);
  PARAM_SET(m_2_3);
  PARAM_SET(m_3_1);
  PARAM_SET(m_3_2);
  PARAM_SET(mu_1);
  PARAM_SET(mu_2);
  PARAM_SET(mu_3);
  PARAM_SET(psi_1);
  PARAM_SET(psi_2);
  PARAM_SET(psi_3);
  PARAM_SET(r_1);
  PARAM_SET(r_2);
  PARAM_SET(r_3);
  if (m != n) err("wrong number of parameters!");
}

template<>
void mtbd3_proc_t::update_IVPs (double *p, int n) {
  int m = 0;
  PARAM_SET(I1_0);
  PARAM_SET(I2_0);
  PARAM_SET(I3_0);
  if (m != n) err("wrong number of initial-value parameters!");
}

template<>
double mtbd3_proc_t::event_rates (double *rate, int n) const {
  int m = 0;
  double total = 0;
  RATE_CALC(params.lambda_1_1 * state.I1);
  RATE_CALC(params.lambda_1_2 * state.I1);
  RATE_CALC(params.lambda_1_3 * state.I1);
  RATE_CALC(params.lambda_2_1 * state.I2);
  RATE_CALC(params.lambda_2_2 * state.I2);
  RATE_CALC(params.lambda_2_3 * state.I2);
  RATE_CALC(params.lambda_3_1 * state.I3);
  RATE_CALC(params.lambda_3_2 * state.I3);
  RATE_CALC(params.lambda_3_3 * state.I3);
  RATE_CALC(params.m_1_2 * state.I1);
  RATE_CALC(params.m_1_3 * state.I1);
  RATE_CALC(params.m_2_1 * state.I2);
  RATE_CALC(params.m_2_3 * state.I2);
  RATE_CALC(params.m_3_1 * state.I3);
  RATE_CALC(params.m_3_2 * state.I3);
  RATE_CALC(params.mu_1 * state.I1);
  RATE_CALC(params.mu_2 * state.I2);
  RATE_CALC(params.mu_3 * state.I3);
  RATE_CALC(params.r_1 * params.psi_1 * state.I1);
  RATE_CALC(params.r_2 * params.psi_2 * state.I2);
  RATE_CALC(params.r_3 * params.psi_3 * state.I3);
  RATE_CALC((1 - params.r_1) * params.psi_1 * state.I1);
  RATE_CALC((1 - params.r_2) * params.psi_2 * state.I2);
  RATE_CALC((1 - params.r_3) * params.psi_3 * state.I3);
  if (m != n) err("wrong number of events!");
  return total;
}

template<>
void mtbd3_genealogy_t::rinit (void) {
  state.I1 = params.I1_0;
state.I2 = params.I2_0;
state.I3 = params.I3_0;
graft(type1,params.I1_0);
graft(type2,params.I2_0);
graft(type3,params.I3_0);
}

template<>
void mtbd3_genealogy_t::jump (int event) {
  switch (event) {
  case 0:
      state.I1 += 1; birth(type1,type1);
      break;
    case 1:
      state.I2 += 1; birth(type1,type2);
      break;
    case 2:
      state.I3 += 1; birth(type1,type3);
      break;
    case 3:
      state.I1 += 1; birth(type2,type1);
      break;
    case 4:
      state.I2 += 1; birth(type2,type2);
      break;
    case 5:
      state.I3 += 1; birth(type2,type3);
      break;
    case 6:
      state.I1 += 1; birth(type3,type1);
      break;
    case 7:
      state.I2 += 1; birth(type3,type2);
      break;
    case 8:
      state.I3 += 1; birth(type3,type3);
      break;
    case 9:
      state.I1 -= 1; state.I2 += 1; migrate(type1,type2);
      break;
    case 10:
      state.I1 -= 1; state.I3 += 1; migrate(type1,type3);
      break;
    case 11:
      state.I2 -= 1; state.I1 += 1; migrate(type2,type1);
      break;
    case 12:
      state.I2 -= 1; state.I3 += 1; migrate(type2,type3);
      break;
    case 13:
      state.I3 -= 1; state.I1 += 1; migrate(type3,type1);
      break;
    case 14:
      state.I3 -= 1; state.I2 += 1; migrate(type3,type2);
      break;
    case 15:
      state.I1 -= 1; death(type1);
      break;
    case 16:
      state.I2 -= 1; death(type2);
      break;
    case 17:
      state.I3 -= 1; death(type3);
      break;
    case 18:
      state.I1 -= 1; sample_death(type1);
      break;
    case 19:
      state.I2 -= 1; sample_death(type2);
      break;
    case 20:
      state.I3 -= 1; sample_death(type3);
      break;
    case 21:
      sample(type1);
      break;
    case 22:
      sample(type2);
      break;
    case 23:
      sample(type3);
      break;
  default:                      // #nocov
    assert(0);                  // #nocov
    break;                      // #nocov
  }
}

GENERICS(MTBD3,mtbd3_genealogy_t)
