##' @name mtbd2_pomp
##' @rdname mtbd2
##' @include mtbd2.R
##' @param x genealogy in \pkg{phylopomp} format.
##' @param pop initial population size (rescales \code{I1_0}, \code{I2_0} proportionally, as in \code{\link{runMTBD2}}'s \code{rinit}).
##' @param obs_tiptype logical; if \code{TRUE}, the true deme (type) at each sample tip is treated as known data and used to constrain the coloring there, rather than marginalized over as latent. This is the "tip type observed" likelihood of KLI Remark 3 (as opposed to the default "tip type unobserved" likelihood, \code{FALSE}), appropriate when the sampled host type (e.g. camel/human) is genuinely part of the data. See Details.
##' @return
##' \code{mtbd2_pomp} returns a \sQuote{pomp} object.
##' @details
##' \code{mtbd2_pomp} constructs a \sQuote{pomp} object containing a given set of data and an MTBD2 model, for exact-likelihood evaluation via \code{\link[pomp]{pfilter}}.
##'
##' \code{obs_tiptype=TRUE} and \code{obs_tiptype=FALSE} are likelihoods of \emph{different} data (the tip-typed vs. the type-unobscured tree) and are not interchangeable; choose according to whether the tip type was actually recorded in the underlying dataset.
##'
##' \code{obs_tiptype=TRUE} constrains the naive per-branch proposal kernel used throughout \pkg{phylopomp}'s coloring simulation (uniform choice among untracked individuals in the source deme) at every one of possibly many sample tips and branch points. On genealogies with many samples this can degenerate at low particle counts \code{Np} (the naive kernel rarely proposes exactly the sequence of colorings compatible with all observed tip types at once, so most particles are discarded and the effective sample size collapses, occasionally yielding an unusably low or \code{-Inf} log-likelihood estimate from a single \code{pfilter} call). Increasing \code{Np} resolves this: a run with tens of samples that fails at \code{Np=4000} typically stabilizes by \code{Np=20000} or more. Always check \code{\link[pomp]{eff_sample_size}} before trusting a single filter run in this mode.
##'
##' Each particle also carries a running count of realized cross-type birth (\eqn{1\to2} and \eqn{2\to1}) events, \code{n12} and \code{n21}, alongside \code{I1}/\code{I2} in every saved state (see \code{\link[pomp]{saved_states}}). These are the natural generalization of Vaughan & Stadler (2025)'s host-to-host spillover counts (e.g. camel-to-human transmissions, whether or not ancestral to the sampled tree) to an arbitrary two-type MTBD application; they cost nothing (they feed back into no rate, proposal, or likelihood computation) and let a single \code{\link[pomp]{pfilter}} call yield a spillover-count posterior directly, in place of the separate stochastic-mapping-plus-trajectory-filter pipeline their approach requires.
##'
##' \code{n12a}/\code{n21a} are the subset of \code{n12}/\code{n21} that involve a currently \emph{tracked} lineage (i.e. that affect the sampled genealogy's own inferred coloring), matching Vaughan & Stadler's distinction between spillover events \emph{ancestral to the sampled dataset} (\code{n12a}/\code{n21a}) and events \emph{in the broader population} (\code{n12}/\code{n21}, which includes \code{n12a}/\code{n21a} as a subset) -- their MERS-CoV application reports both quantities separately.
##'
##' \code{sample_window_1}/\code{sample_window_2}, if given, restrict deme-\eqn{i} sampling to \code{c(first,last)} (in the same time units as \code{x}'s branch lengths): \code{psi_i} applies only within that window, and \code{mu_i} absorbs the rest of \eqn{\delta_i=\mu_i+\psi_i} outside it, matching Vaughan & Stadler (2025)'s time-varying host-specific sampling proportion (their \eqn{s_i(t)}, forced to zero outside the interval between each host type's first and last sample). \eqn{\delta_i} itself does not vary -- only its split between \code{mu_i} and \code{psi_i} does. The default \code{NULL} (sampling active at every \eqn{t}) exactly reproduces the model's previous, window-free behavior; see \code{\link{mtbd2_vs2025_params}} for the corresponding \eqn{(R,\delta,s)}-scale parameterization. A sample observed outside its deme's declared window is a zero-likelihood event (\eqn{\psi_i(t)=0} there), not silently ignored.
##' @param sample_window_1,sample_window_2 \code{NULL}, or a length-2 numeric \code{c(first,last)} giving deme \eqn{i}'s sampling-active window. See Details.
##' @details
##' \code{season_amp_1}, \code{season_period_1}, \code{season_phase_1} apply
##' seasonal (cosine) forcing to \code{lambda_1_1} only (camel-to-camel
##' transmission, e.g. tied to camel calving season), leaving
##' \code{lambda_1_2}/\code{lambda_2_1}/\code{lambda_2_2} constant:
##' \deqn{\lambda_{11}(t) = \lambda_{1,1}\bigl(1+\code{season\_amp\_1}\cdot
##' \cos(2\pi(t-\code{season\_phase\_1})/\code{season\_period\_1})\bigr).}
##' The default \code{season_amp_1=0} makes this identically \code{lambda_1_1}
##' regardless of the period/phase values, exactly reproducing the
##' non-seasonal model. \code{season_amp_1} must lie in \eqn{[0,1)} so
##' \eqn{\lambda_{11}(t)} stays nonnegative. \code{season_period_1} should
##' ordinarily be \emph{fixed} (e.g. at \code{1} if \code{x}'s branch lengths
##' are in years), not included in \code{\link{mtbd2_mif2}}'s or
##' \code{\link{mtbd2_pmcmc}}'s \code{estimate}: unlike amplitude and phase,
##' period is typically poorly identified from a single phylodynamic tree,
##' and letting it float invites a multimodal likelihood surface with peaks
##' at harmonics of the true period. \code{season_phase_1} is only
##' interpretable in calendar terms (e.g. "transmission peaks in December")
##' if \code{t0} is itself anchored to a real calendar time rather than an
##' arbitrary origin -- nothing here enforces that; it is the caller's
##' responsibility when constructing \code{x}.
##' @param season_amp_1 numeric in \eqn{[0,1)}; seasonal amplitude for \code{lambda_1_1}. Default \code{0} (no seasonality). See Details.
##' @param season_period_1 numeric > 0; seasonal period for \code{lambda_1_1}, in the same time units as \code{x}'s branch lengths. See Details.
##' @param season_phase_1 numeric; seasonal phase (time of peak \code{lambda_1_1}) for \code{lambda_1_1}. See Details.
##' @importFrom pomp pomp onestep
##' @export
mtbd2_pomp <- function (
  x,
  lambda_1_1, lambda_1_2, lambda_2_1, lambda_2_2,
  m_1_2, m_2_1,
  mu_1, mu_2,
  psi_1, psi_2,
  r_1, r_2,
  I1_0, I2_0, pop,
  obs_tiptype = FALSE,
  sample_window_1 = NULL, sample_window_2 = NULL,
  season_amp_1 = 0, season_period_1 = 1, season_phase_1 = 0
)
{
  chk_window <- function (w, nm) {
    if (is.null(w)) return(c(-Inf,Inf))
    if (!is.numeric(w) || length(w) != 2L || w[1L] > w[2L])
      pStop(sQuote(nm)," must be NULL or a length-2 numeric c(first,last) ",
        "with first <= last.")
    w
  }
  win1 <- chk_window(sample_window_1,"sample_window_1")
  win2 <- chk_window(sample_window_2,"sample_window_2")
  if (season_amp_1 < 0 || season_amp_1 >= 1)
    pStop(sQuote("season_amp_1")," must lie in [0,1).")
  if (season_period_1 <= 0)
    pStop(sQuote("season_period_1")," must be positive.")
  ## NB: gendat()'s default obscure=TRUE is required here (obscure=FALSE
  ## retains migration nodes as real tree nodes, which mtbd2_gill's
  ## nodetype switch has no case for). The true per-sample deme survives
  ## obscuring regardless: genealogy_t::obscure() zeroes deme only on
  ## black balls and node-level deme, but gendat.cc reads a sample row's
  ## deme from its blue ball, which obscure() never touches. So gi$deme
  ## already carries the true tip type at every nodetype==1 row even
  ## with the default obscure=TRUE.
  x |> gendat() -> gi
  ivps <- structure(c(I1_0,I2_0),names=c("I1_0","I2_0"))
  if (any(ivps < 0))
    pStop(paste(sQuote(names(ivps)),collapse=","),
      " must be nonnegative.")
  if (any(c(r_1,r_2) < 0 | c(r_1,r_2) > 1))
    pStop(sQuote("r_1")," and ",sQuote("r_2")," must lie in [0,1].")
  pomp(
    data=NULL,
    t0=gi$nodetime[1L],
    times=gi$nodetime[-1L],
    params=c(
      lambda_11=lambda_1_1,lambda_12=lambda_1_2,
      lambda_21=lambda_2_1,lambda_22=lambda_2_2,
      m12=m_1_2,m21=m_2_1,
      mu1=mu_1,mu2=mu_2,
      psi1=psi_1,psi2=psi_2,
      r1=r_1,r2=r_2,
      pop=pop,ivps,
      obstype=as.numeric(obs_tiptype),
      first1=win1[1L],last1=win1[2L],
      first2=win2[1L],last2=win2[2L],
      l11amp=season_amp_1,l11period=season_period_1,l11phase=season_phase_1
    ),
    userdata=gi,
    nstatevars=10L + gi$nsample,
    rinit="mtbd2_rinit",
    rprocess=onestep("mtbd2_gill"),
    dmeasure="mtbd2_dmeas",
    statenames=c(
      "I1","I2","ll","node","ell1","ell2","n12","n21","n12a","n21a","color"
    ),
    paramnames=c(
      "lambda_11","lambda_12","lambda_21","lambda_22",
      "m12","m21","mu1","mu2","psi1","psi2","r1","r2",
      "pop","I1_0","I2_0","obstype",
      "first1","last1","first2","last2",
      "l11amp","l11period","l11phase"
    ),
    PACKAGE="phylopomp"
  )
}
