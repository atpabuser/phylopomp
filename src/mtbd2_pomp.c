// mtbd2_pomp.c --- exact-likelihood pomp bridge for MTBD2
//
// Two-type general linear multitype birth-death-sampling process
// (Kuehnert et al. 2016; Vaughan & Stadler 2025), following the King,
// Lin & Ionides (2024) structured Markov genealogy process framework.
//
// Generalizes three already-validated templates in this package:
//   - same-type and cross-type (mixed-fork) birth: src/bdss_pomp.c
//   - migration, always regular (never a singular/observed node): src/bdei_pomp.c
//   - dual-channel sampling (destructive + non-destructive),
//     marginalized at ambiguous terminal leaves: src/seirs_pomp.c
//
// See the companion document mtbd2_filter.tex/.pdf for the derivation
// and an audit of what was checked numerically and against what.
//
// Tip-type-observed mode (parameter 'obstype', 0/1): by default the
// coloring is fully latent and marginalized, including at sample tips
// (Vaughan & Stadler's edge-typed-tree distinction is collapsed into
// the filter). With obstype=1, the true deme at each sample (as known
// from the data, e.g. camel/human host species) is read from userdata
// and used to constrain rather than marginalize the coloring there --
// analogous to KLI's Remark 3 / the "tip type observed" likelihood
// described in bdss_filter_v2.pdf's "Sample metadata" box. The true
// per-sample deme survives gendat()'s default obscure=TRUE (see
// R/mtbd2_pomp.R for why); no special gendat() call is needed here.

#include "pomplink.h"
#include "internal.h"

#define Type1 1
#define Type2 2

static const int mtbd2_nrate = 12;

static inline int mtbd2_random_choice (double n) {
  return floor(R_unif_index(n));
}

static void mtbd2_change_color (double *color, int nsample,
                          int n, int from, int to) {
  int i = -1;
  while (n >= 0 && i < nsample) {
    i++;
    if (!ISNA(color[i]) && nearbyint(color[i]) == from) n--;
  }
  assert(i < nsample);
  assert(n == -1);
  assert(nearbyint(color[i]) == from);
  color[i] = to;
}

#define lambda_11 (__p[__parindex[0]])
#define lambda_12 (__p[__parindex[1]])
#define lambda_21 (__p[__parindex[2]])
#define lambda_22 (__p[__parindex[3]])
#define m12       (__p[__parindex[4]])
#define m21       (__p[__parindex[5]])
#define mu1       (__p[__parindex[6]])
#define mu2       (__p[__parindex[7]])
#define psi1      (__p[__parindex[8]])
#define psi2      (__p[__parindex[9]])
#define r1        (__p[__parindex[10]])
#define r2        (__p[__parindex[11]])
#define POP       (__p[__parindex[12]])
#define I10       (__p[__parindex[13]])
#define I20       (__p[__parindex[14]])
#define OBSTYPE   (__p[__parindex[15]])
#define I1        (__x[__stateindex[0]])
#define I2        (__x[__stateindex[1]])
#define ll        (__x[__stateindex[2]])
#define node      (__x[__stateindex[3]])
#define ell1      (__x[__stateindex[4]])
#define ell2      (__x[__stateindex[5]])
// N12/N21: per-particle running counts of realized 1->2 / 2->1 cross-type
// birth events (spillovers), both regular and singular (branch-point).
// Auxiliary bookkeeping only -- they do not feed back into any rate, pi, or
// ll computation, so they cannot affect the likelihood; they exist purely
// so saved_states()/mtbd2_trajectory_posterior() can report a spillover-
// count posterior (Vaughan & Stadler 2025's headline result) as a byproduct
// of the same pfilter() call, rather than requiring a separate post-hoc
// event-counting pass.
#define N12       (__x[__stateindex[6]])
#define N21       (__x[__stateindex[7]])
#define COLOR     (__x[__stateindex[8]])

#define MTBD2_EVENT_RATES                               \
  mtbd2_event_rates(__x,__p,t,                           \
              __stateindex,__parindex,__covindex,        \
              __covars,rate,logpi,&penalty)

// Event order (12 events):
//  0: birth 1->1 (same-type, within-deme fork mass excluded via 'disc')
//  1: birth 2->2 (same-type)
//  2: birth 1->2, source deme-1 individual untracked  (no recolor)
//  3: birth 1->2, source deme-1 individual tracked     (recolor 1->2)
//  4: birth 2->1, source deme-2 individual untracked  (no recolor)
//  5: birth 2->1, source deme-2 individual tracked     (recolor 2->1)
//  6: migration 1->2, source untracked (no recolor)
//  7: migration 1->2, source tracked    (recolor 1->2)
//  8: migration 2->1, source untracked (no recolor)
//  9: migration 2->1, source tracked    (recolor 2->1)
// 10: death in deme 1
// 11: death in deme 2
// Sampling (both channels) is always singular -- see the 'case 1' block
// in mtbd2_gill -- and contributes only a continuous decay term here.
static double mtbd2_event_rates
(
 double *__x,
 const double *__p,
 double t,
 const int *__stateindex,
 const int *__parindex,
 const int *__covindex,
 const double *__covars,
 double *rate,
 double *logpi,
 double *penalty
 ) {
  double event_rate = 0;
  double alpha, pi, disc;
  *penalty = 0;
  assert(I1 >= ell1 && ell1 >= 0);
  assert(I2 >= ell2 && ell2 >= 0);
  // 0: 1->1 birth (within-type, disc = missing branch-point mass)
  alpha = lambda_11*I1;
  disc = (I1 > 0) ? ell1*(ell1-1)/I1/(I1+1) : 1;
  *penalty += alpha*disc;
  event_rate += (*rate = alpha*(1-disc)); rate++;
  *logpi = 0; logpi++;
  // 1: 2->2 birth (within-type)
  alpha = lambda_22*I2;
  disc = (I2 > 0) ? ell2*(ell2-1)/I2/(I2+1) : 1;
  *penalty += alpha*disc;
  event_rate += (*rate = alpha*(1-disc)); rate++;
  *logpi = 0; logpi++;
  // 2,3: 1->2 birth (cross-type; source deme is 1)
  // Support-safe (uniform) proposal. mtbd2_filter.tex's "Formal
  // derivation-matching audit" confirms the case-2/case-3 terms in
  // mtbd2_gill (log(1-ell2/I2) and log(1-ell1/I1)-log(I2)) already equal
  // log(Phi) exactly (the -logpi[event] subtracted before the switch
  // supplies the -log(pi) half of the Phi/pi ratio) -- that identity holds
  // for *any* valid positive, ell1-summing-to-1 proposal, so pi need not
  // equal Phi itself and can be changed independently of those terms. The
  // audited "naive kernel" pi = (I1-ell1)/I1 (a specific individual drawn
  // uniformly from the SOURCE deme) could hit exactly zero whenever deme 1
  // was fully tracked (I1==ell1), even though Phi for this branch is a
  // deme-2-side quantity (per Corollary "Derivation of the migration
  // weight", reused here for the untracked-parent mixed-birth branch) that
  // stays positive there -- the audit checked algebraic self-consistency,
  // not this boundary/support case (same failure mode identified and fixed
  // in bdei_pomp.c/bdss_pomp.c per bdei_filter_v2.pdf/bdss_filter_v2.pdf
  // Sec.9 "Support-safe choices", both dated after this audit).
  alpha = lambda_12*I1;
  pi = 1/(ell1+1);
  event_rate += (*rate = alpha*pi); rate++;
  *logpi = log(pi); logpi++;
  pi = ell1/(ell1+1);
  event_rate += (*rate = alpha*pi); rate++;
  *logpi = log(pi)-log(ell1); logpi++;
  // 4,5: 2->1 birth (cross-type; source deme is 2) -- mirror of 2,3.
  alpha = lambda_21*I2;
  pi = 1/(ell2+1);
  event_rate += (*rate = alpha*pi); rate++;
  *logpi = log(pi); logpi++;
  pi = ell2/(ell2+1);
  event_rate += (*rate = alpha*pi); rate++;
  *logpi = log(pi)-log(ell2); logpi++;
  // 6,7: migration 1->2 (source deme is 1)
  // Same support-safety fix. mtbd2_filter.tex Sec."Migration: Always
  // Regular", Corollary "Derivation of the migration weight" (and,
  // independently, mtbd_companion.pdf Sec.4.4/Definition 4.2 for the
  // generic model) confirm the true Phi for migrate12 is a deme-2
  // (destination, post-event) quantity -- (n2-ell2)/n2 regular, 1/n2
  // singular -- not a function of deme 1 at all; the case-6/case-7 terms in
  // mtbd2_gill already encode it, verified by the same formal audit. Only
  // the proposal (which *did* use deme 1's occupancy, per the audited
  // naive kernel) needed the support-safety fix, same as events 2-5 above.
  alpha = m12*I1;
  pi = 1/(ell1+1);
  event_rate += (*rate = alpha*pi); rate++;
  *logpi = log(pi); logpi++;
  pi = ell1/(ell1+1);
  event_rate += (*rate = alpha*pi); rate++;
  *logpi = log(pi)-log(ell1); logpi++;
  // 8,9: migration 2->1 (source deme is 2) -- mirror of 6,7.
  alpha = m21*I2;
  pi = 1/(ell2+1);
  event_rate += (*rate = alpha*pi); rate++;
  *logpi = log(pi); logpi++;
  pi = ell2/(ell2+1);
  event_rate += (*rate = alpha*pi); rate++;
  *logpi = log(pi)-log(ell2); logpi++;
  // 10: death in deme 1
  alpha = mu1*I1;
  if (I1 > ell1) {
    event_rate += (*rate = alpha); rate++;
    *logpi = 0; logpi++;
  } else {
    *rate = 0; rate++;
    *logpi = 0; logpi++;
    *penalty += alpha;
  }
  // 11: death in deme 2
  alpha = mu2*I2;
  if (I2 > ell2) {
    event_rate += (*rate = alpha); rate++;
    *logpi = 0; logpi++;
  } else {
    *rate = 0; rate++;
    *logpi = 0; logpi++;
    *penalty += alpha;
  }
  // sampling (both channels; always singular -- decay only)
  *penalty += psi1*I1 + psi2*I2;
  assert(R_FINITE(event_rate));
  return event_rate;
}

void mtbd2_rinit
(
 double *__x,
 const double *__p,
 double t0,
 const int *__stateindex,
 const int *__parindex,
 const int *__covindex,
 const double *__covars
 ){
  double adj = POP/(I10+I20);
  I1 = nearbyint(I10*adj);
  I2 = nearbyint(I20*adj);
  ell1 = 0;
  ell2 = 0;
  N12 = 0;
  N21 = 0;
  ll = 0;
  node = 0;
}

void mtbd2_gill
(
 double *__x,
 const double *__p,
 const int *__stateindex,
 const int *__parindex,
 const int *__covindex,
 const double *__covars,
 double t,
 double dt
 ){
  double tstep = 0.0, tmax = t + dt;
  double *color = &COLOR;
  const int nsample = *get_userdata_int("nsample");
  const int *nodetype = get_userdata_int("nodetype");
  const int *lineage = get_userdata_int("lineage");
  const int *sat = get_userdata_int("saturation");
  const int *index = get_userdata_int("index");
  const int *child = get_userdata_int("child");
  // 'deme' holds the TRUE type at every node (see gendat(obscure=FALSE)),
  // meaningful at sample nodes only; used only when OBSTYPE != 0.
  const int *truedeme = get_userdata_int("deme");
  int obstype = (int) nearbyint(OBSTYPE);

  int parent = (int) nearbyint(node);

#ifndef NDEBUG
  int nnode = *get_userdata_int("nnode");
  assert(parent>=0);
  assert(parent<=nnode);
#endif

  int parlin = lineage[parent];
  int parcol = color[parlin];
  assert(parlin >= 0 && parlin < nsample);

  ll = 0;

  switch (nodetype[parent]) {
  default:
    break;
  case 0:                       // root
    assert(sat[parent]==1);
    {
      int c = child[index[parent]];
      assert(lineage[parent]==lineage[c]);
      if (I1-ell1 + I2-ell2 > 0) {
        double x = (I1-ell1)/(I1-ell1 + I2-ell2);
        if (unif_rand() < x) {
          color[lineage[c]] = Type1;
          ell1 += 1;
          ll -= log(x);
        } else {
          color[lineage[c]] = Type2;
          ell2 += 1;
          ll -= log(1-x);
        }
      } else {
        ll += R_NegInf;
        if (unif_rand() < 0.5) {
          color[lineage[c]] = Type1;
          ell1 += 1; I1 += 1;
        } else {
          color[lineage[c]] = Type2;
          ell2 += 1; I2 += 1;
        }
      }
    }
    break;
  case 1:                       // sample (destructive and/or non-destructive)
    if (obstype) {
      // Tip-type-observed mode (KLI Remark 3; bdss_filter_v2.pdf Sec.
      // "Sample metadata"): the true deme at this sample is known data,
      // not latent. If the particle's own simulated coloring disagrees,
      // this whole coloring history is incompatible with the data: kill
      // the particle (ll = R_NegInf) and force-correct the tracked state
      // so later steps stay internally consistent (mirrors bdei_pomp.c's
      // case 1 self-heal).
      int truetype = truedeme[parent];
      if (parcol != truetype) {
        ll += R_NegInf;
        if (parcol == Type1) {
          ell1 -= 1; ell2 += 1; I1 -= 1; I2 += 1;
        } else if (parcol == Type2) {
          ell2 -= 1; ell1 += 1; I2 -= 1; I1 += 1;
        }
        color[parlin] = truetype;
        parcol = truetype;
      }
    }
    if (sat[parent] == 1) {     // unambiguous non-destructive inline node
      int c = child[index[parent]];
      if (parcol == Type1) {
        color[lineage[c]] = Type1;
        ll += (I1 > 0) ? log((1-r1)*psi1) : R_NegInf;
      } else if (parcol == Type2) {
        color[lineage[c]] = Type2;
        ll += (I2 > 0) ? log((1-r2)*psi2) : R_NegInf;
      } else {
        ll += R_NegInf;
      }
    } else if (sat[parent] == 0) { // ambiguous leaf: marginalize destructive/non-destructive
      if (parcol == Type1) {
        ell1 -= 1;
        ll += (I1 > 0) ? log(psi1) : R_NegInf;
        if (unif_rand() < r1) {        // destructive, proposal prob r1
          ll += log(I1);
          I1 -= 1;
        } else {                        // non-destructive, proposal prob 1-r1
          ll += log(I1-ell1);
        }
      } else if (parcol == Type2) {
        ell2 -= 1;
        ll += (I2 > 0) ? log(psi2) : R_NegInf;
        if (unif_rand() < r2) {
          ll += log(I2);
          I2 -= 1;
        } else {
          ll += log(I2-ell2);
        }
      } else {
        ll += R_NegInf;
      }
    } else {
      assert(0);                // #nocov
      ll += R_NegInf;           // #nocov
    }
    color[parlin] = R_NaReal;
    break;
  case 2:                       // branch point
    assert(sat[parent]==2);
    {
      int c1 = child[index[parent]];
      int c2 = child[index[parent]+1];
      assert(c1 != c2);
      assert(lineage[c1] != lineage[c2]);
      assert(lineage[c1] != parlin || lineage[c2] != parlin);
      assert(lineage[c1] == parlin || lineage[c2] == parlin);

      if (parcol == Type1) {
        double lam = lambda_11 + lambda_12;
        if (lam <= 0) { ll += R_NegInf; break; }
        double p_11 = lambda_11/lam;
        if (unif_rand() < p_11) {
          // 1->1 birth: both children type1
          I1 += 1; ell1 += 1;
          ll += log(2*lam/I1);
          if (I1*(I1-1) <= 0) ll += R_NegInf;
          color[lineage[c1]] = Type1;
          color[lineage[c2]] = Type1;
        } else {
          // 1->2 birth: one type1, one type2
          I2 += 1; ell2 += 1; N12 += 1;
          ll += log(lam/I2);
          if (I1*I2 <= 0) ll += R_NegInf;
          if (unif_rand() < 0.5) {
            color[lineage[c1]] = Type1;
            color[lineage[c2]] = Type2;
          } else {
            color[lineage[c1]] = Type2;
            color[lineage[c2]] = Type1;
          }
          ll -= log(0.5);
        }
      } else if (parcol == Type2) {
        double lam = lambda_21 + lambda_22;
        if (lam <= 0) { ll += R_NegInf; break; }
        double p_22 = lambda_22/lam;
        if (unif_rand() < p_22) {
          // 2->2 birth: both children type2
          I2 += 1; ell2 += 1;
          ll += log(2*lam/I2);
          if (I2*(I2-1) <= 0) ll += R_NegInf;
          color[lineage[c1]] = Type2;
          color[lineage[c2]] = Type2;
        } else {
          // 2->1 birth: one type2, one type1
          I1 += 1; ell1 += 1; N21 += 1;
          ll += log(lam/I1);
          if (I1*I2 <= 0) ll += R_NegInf;
          if (unif_rand() < 0.5) {
            color[lineage[c1]] = Type1;
            color[lineage[c2]] = Type2;
          } else {
            color[lineage[c1]] = Type2;
            color[lineage[c2]] = Type1;
          }
          ll -= log(0.5);
        }
      } else {
        ll += R_NegInf;
      }
    }
    break;
  }

  // continuous portion of filter equation
  if (tmax > t && R_FINITE(ll)) {

    double rate[mtbd2_nrate], logpi[mtbd2_nrate];
    int event;
    double event_rate = 0;
    double penalty = 0;

    event_rate = MTBD2_EVENT_RATES;
    tstep = exp_rand()/event_rate;

    while (t + tstep < tmax) {
      event = rcateg(event_rate,rate,mtbd2_nrate);
      assert(event>=0 && event<mtbd2_nrate);
      ll -= penalty*tstep + logpi[event];
      switch (event) {
      case 0:                   // 1->1 birth
        I1 += 1;
        assert(!ISNAN(ll));
        break;
      case 1:                   // 2->2 birth
        I2 += 1;
        assert(!ISNAN(ll));
        break;
      case 2:                   // 1->2 birth, source untracked
        I2 += 1; N12 += 1;
        ll += log(1-ell2/I2);
        assert(!ISNAN(ll));
        break;
      case 3:                   // 1->2 birth, source tracked -> recolor
        mtbd2_change_color(color,nsample,mtbd2_random_choice(ell1),Type1,Type2);
        ell1 -= 1; ell2 += 1;
        I2 += 1; N12 += 1;
        ll += log(1-ell1/I1)-log(I2);
        assert(!ISNAN(ll));
        break;
      case 4:                   // 2->1 birth, source untracked
        I1 += 1; N21 += 1;
        ll += log(1-ell1/I1);
        assert(!ISNAN(ll));
        break;
      case 5:                   // 2->1 birth, source tracked -> recolor
        mtbd2_change_color(color,nsample,mtbd2_random_choice(ell2),Type2,Type1);
        ell2 -= 1; ell1 += 1;
        I1 += 1; N21 += 1;
        ll += log(1-ell2/I2)-log(I1);
        assert(!ISNAN(ll));
        break;
      case 6:                   // migration 1->2, source untracked
        I1 -= 1; I2 += 1;
        ll += log(1-ell2/I2);
        assert(!ISNAN(ll));
        break;
      case 7:                   // migration 1->2, source tracked -> recolor
        mtbd2_change_color(color,nsample,mtbd2_random_choice(ell1),Type1,Type2);
        ell1 -= 1; ell2 += 1;
        I1 -= 1; I2 += 1;
        ll -= log(I2);
        assert(!ISNAN(ll));
        break;
      case 8:                   // migration 2->1, source untracked
        I2 -= 1; I1 += 1;
        ll += log(1-ell1/I1);
        assert(!ISNAN(ll));
        break;
      case 9:                   // migration 2->1, source tracked -> recolor
        mtbd2_change_color(color,nsample,mtbd2_random_choice(ell2),Type2,Type1);
        ell2 -= 1; ell1 += 1;
        I2 -= 1; I1 += 1;
        ll -= log(I1);
        assert(!ISNAN(ll));
        break;
      case 10:                  // death in deme 1
        assert(I1>=1);
        I1 -= 1;
        break;
      case 11:                  // death in deme 2
        assert(I2>=1);
        I2 -= 1;
        break;
      default:                  // #nocov
        assert(0);              // #nocov
        break;                  // #nocov
      }

      ell1 = nearbyint(ell1);
      ell2 = nearbyint(ell2);

      t += tstep;
      event_rate = MTBD2_EVENT_RATES;
      tstep = exp_rand()/event_rate;

    }
    tstep = tmax - t;
    ll -= penalty*tstep;
  }
  node += 1;
}

# define lik  (__lik[0])

void mtbd2_dmeas
(
 double *__lik,
 const double *__y,
 const double *__x,
 const double *__p,
 int give_log,
 const int *__obsindex,
 const int *__stateindex,
 const int *__parindex,
 const int *__covindex,
 const double *__covars,
 double t
 ) {
  assert(!ISNAN(ll));
  lik = (give_log) ? ll : exp(ll);
}
