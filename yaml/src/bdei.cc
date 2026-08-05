// BDEI: Linear birth-death with exposed and infectious classes (C++)
#include "master.h"
#include "popul_proc.h"
#include "generics.h"
#include "internal.h"

static const int exposed = 1;
static const int infectious = 2;

//! BDEI process state.
typedef struct {
  int E;
  int I;
} bdei_state_t;

//! BDEI process parameters.
typedef struct {
  double sigma;
  double lambda;
  double mu;
  double chi;
  double pop;
  double E0;
  double I0;
} bdei_parameters_t;

using bdei_proc_t = popul_proc_t<bdei_state_t,bdei_parameters_t,4>;
using bdei_genealogy_t = master_t<bdei_proc_t,2>;

template<>
std::string bdei_proc_t::yaml (std::string tab) const {
  std::string t = tab + "  ";
  std::string p = tab + "parameter:\n"
    + YAML_PARAM(sigma)
    + YAML_PARAM(lambda)
    + YAML_PARAM(mu)
    + YAML_PARAM(chi)
    + YAML_PARAM(pop)
    + YAML_PARAM(E0)
    + YAML_PARAM(I0);
  std::string s = tab + "state:\n"
    + YAML_STATE(E)
    + YAML_STATE(I);
  return p+s;
}

template<>
void bdei_proc_t::update_params (double *p, int n) {
  int m = 0;
  PARAM_SET(sigma);
  PARAM_SET(lambda);
  PARAM_SET(mu);
  PARAM_SET(chi);
  if (m != n) err("wrong number of parameters!");
}

template<>
void bdei_proc_t::update_IVPs (double *p, int n) {
  int m = 0;
  PARAM_SET(pop);
  PARAM_SET(E0);
  PARAM_SET(I0);
  if (m != n) err("wrong number of initial-value parameters!");
}

template<>
double bdei_proc_t::event_rates (double *rate, int n) const {
  int m = 0;
  double total = 0;
  RATE_CALC(params.sigma * state.E);
  RATE_CALC(params.lambda * state.I);
  RATE_CALC(params.mu * state.I);
  RATE_CALC(params.chi * state.I);
  if (m != n) err("wrong number of events!");
  return total;
}

template<>
void bdei_genealogy_t::rinit (void) {
  double m = params.pop/(params.E0 + params.I0);
state.E = nearbyint(m*params.E0);
state.I = nearbyint(m*params.I0);
graft(exposed, state.E);
graft(infectious, state.I);
}

template<>
void bdei_genealogy_t::jump (int event) {
  switch (event) {
  case 0:
      state.E -= 1; state.I += 1; migrate(exposed, infectious);
      break;
    case 1:
      state.E += 1; birth(infectious, exposed);
      break;
    case 2:
      state.I -= 1; death(infectious);
      break;
    case 3:
      state.I -= 1; sample_death(infectious);
      break;
  default:                      // #nocov
    assert(0);                  // #nocov
    break;                      // #nocov
  }
}

GENERICS(BDEI,bdei_genealogy_t)
