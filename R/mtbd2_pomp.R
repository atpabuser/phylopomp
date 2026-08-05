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
  obs_tiptype = FALSE
)
{
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
      obstype=as.numeric(obs_tiptype)
    ),
    userdata=gi,
    nstatevars=6L + gi$nsample,
    rinit="mtbd2_rinit",
    rprocess=onestep("mtbd2_gill"),
    dmeasure="mtbd2_dmeas",
    statenames=c(
      "I1","I2","ll","node","ell1","ell2","color"
    ),
    paramnames=c(
      "lambda_11","lambda_12","lambda_21","lambda_22",
      "m12","m21","mu1","mu2","psi1","psi2","r1","r2",
      "pop","I1_0","I2_0","obstype"
    ),
    PACKAGE="phylopomp"
  )
}
