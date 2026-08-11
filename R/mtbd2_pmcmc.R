##' Particle-marginal Metropolis-Hastings for MTBD2
##'
##' A minimal particle-marginal Metropolis-Hastings (PMMH) sampler for
##' \code{\link{mtbd2_pomp}}, giving a Bayesian posterior over rate
##' parameters (and, via \code{\link{mtbd2_trajectory_posterior}} or
##' \code{\link[pomp]{saved_states}} on the retained chain, population
##' trajectories and cross-type-birth event counts) comparable in kind to
##' Vaughan & Stadler (2025)'s BEAST2/BDMM-Prime posterior -- \pkg{phylopomp}
##' otherwise ships only \code{\link{mtbd2_mif2}} (a point-estimate/MLE
##' method), not a Bayesian sampler.
##'
##' @name mtbd2_pmcmc
##' @rdname mtbd2_pmcmc
##' @include mtbd2_infer.R
##' @param x genealogy in \pkg{phylopomp} format.
##' @param start named list of starting parameter values, in the same form
##'   as \code{\link{mtbd2_loglik}}'s \code{params}.
##' @param estimate character vector of parameter names (user-facing names,
##'   as in \code{start}) to update by random-walk Metropolis. All must be
##'   strictly positive in \code{start}: the random walk is on the log
##'   scale, exactly as in \code{\link{mtbd2_mif2}}. Parameters not named
##'   here are held fixed at their \code{start} value throughout.
##' @param log_prior function taking a named numeric vector (the current
##'   values of \code{estimate}, in \pkg{phylopomp}'s native rate
##'   parameterization -- i.e. \emph{not} Vaughan & Stadler's
##'   \eqn{(R,\delta,s)} scale; if sampling in that scale, include the
##'   Jacobian of \code{\link{mtbd2_vs2025_params}} in \code{log_prior})
##'   and returning a single log-density value (\code{-Inf} for
##'   zero-density regions, e.g. a hard prior bound). Default: flat
##'   (improper) prior, \code{function(p) 0}.
##' @param proposal_sd numeric; random-walk standard deviation on the log
##'   scale. Either a single value applied to all of \code{estimate}, or a
##'   named vector.
##' @param Np integer; particles per \code{\link[pomp]{pfilter}} call.
##' @param Niter integer; number of Metropolis-Hastings iterations.
##' @param nreps integer; independent \code{pfilter} replicates averaged
##'   (via \code{\link[pomp]{logmeanexp}}) into a single likelihood estimate
##'   per iteration, both at the current and proposed parameter values. This
##'   reduces the variance of the pseudo-marginal likelihood estimator (at
##'   \code{nreps}-fold cost) but does not change its unbiasedness for
##'   \eqn{L} (not \eqn{\log L}), which is all pseudo-marginal MCMC requires;
##'   \code{nreps=1} is the standard PMMH setting.
##' @param verbose logical; report progress and running acceptance rate.
##' @details
##' This is deliberately minimal: a single random-walk block-update (not
##' adaptive, not delayed-acceptance, no parallel tempering), because
##' \pkg{phylopomp}'s SMC likelihood is itself already the expensive part of
##' any such scheme, and Vaughan \& Stadler's own comparison target
##' (BDMM-Prime via BEAST2) is a mature, heavily-optimized sampler that this
##' is not attempting to match in wall-clock efficiency -- only to provide
##' \emph{a} valid Bayesian sampler so a posterior comparison is possible at
##' all. Tune \code{proposal_sd} and \code{Np} by hand from the trace and
##' acceptance rate (aim for roughly 20-40\% acceptance; the \pkg{coda}
##' package's \code{effectiveSize}, if installed, is a convenient external
##' diagnostic on \code{chain}); a poor \code{Np} choice inflates the variance of
##' the likelihood ratio and can either stall the chain (too small) or waste
##' computation (unnecessarily large). As with \code{\link{mtbd2_mif2}},
##' run multiple independent chains from dispersed starting points before
##' trusting convergence.
##' @return
##' \code{mtbd2_pmcmc} returns a list with elements \code{chain} (a
##' \code{\link[tibble]{tibble}}, one row per iteration, with columns for
##' each parameter in \code{estimate} plus \code{loglik}, \code{log_prior},
##' and \code{accept} (logical, whether that iteration's proposal was
##' accepted -- the row's parameter values are the post-acceptance-decision
##' state either way)), \code{acceptance_rate} (the overall fraction
##' accepted), and \code{pomp} (the final-iteration \sQuote{pomp} object,
##' for follow-up \code{\link[pomp]{pfilter}}/\code{\link[pomp]{saved_states}}
##' calls at the chain's last parameter value).
##' @importFrom tibble tibble
##' @importFrom pomp pfilter logLik logmeanexp
##' @importFrom stats rnorm runif
##' @export
mtbd2_pmcmc <- function (
  x, start, estimate,
  log_prior = function (p) 0,
  proposal_sd = 0.02,
  Np = 2000, Niter = 1000, nreps = 1,
  verbose = FALSE
) {
  bad <- setdiff(estimate,names(mtbd2_parmap))
  if (length(bad) > 0L)
    pStop("mtbd2_pmcmc",": not estimable parameter(s): ",
      paste(sQuote(bad),collapse=","),".")
  missing <- setdiff(estimate,names(start))
  if (length(missing) > 0L)
    pStop("mtbd2_pmcmc",": ",paste(sQuote(missing),collapse=",")," absent from ",
      sQuote("start"),".")
  vals <- unlist(start[estimate])
  if (any(vals <= 0))
    pStop("mtbd2_pmcmc",": the random walk is on the log scale, so all of ",
      sQuote("estimate")," must be strictly positive in ",sQuote("start"),
      "; offending: ",paste(sQuote(estimate[vals<=0]),collapse=","),".")
  if (length(proposal_sd) == 1L && is.null(names(proposal_sd))) {
    sds <- stats::setNames(rep(proposal_sd,length(estimate)),estimate)
  } else {
    miss <- setdiff(estimate,names(proposal_sd))
    if (length(miss) > 0L)
      pStop("mtbd2_pmcmc",": ",sQuote("proposal_sd")," lacks entries for ",
        paste(sQuote(miss),collapse=","),".")
    sds <- proposal_sd[estimate]
  }

  ll_estimate <- function (theta) {
    p <- start
    p[estimate] <- as.list(theta)
    po <- do.call(mtbd2_pomp,c(list(x=x),p))
    lls <- replicate(nreps,pomp::logLik(pomp::pfilter(po,Np=Np)))
    ll <- if (all(!is.finite(lls))) -Inf else
      unname(pomp::logmeanexp(lls[is.finite(lls)]))
    list(loglik=ll,pomp=po)
  }

  theta <- unlist(start[estimate])
  cur <- ll_estimate(theta)
  cur_lp <- log_prior(theta)
  if (!is.finite(cur$loglik + cur_lp))
    pStop("mtbd2_pmcmc",": ",sQuote("start")," has zero (or unevaluable) ",
      "posterior density; choose a different starting point.")

  chain <- matrix(NA_real_,nrow=Niter,ncol=length(estimate)+3L,
    dimnames=list(NULL,c(estimate,"loglik","log_prior","accept")))
  naccept <- 0L
  po_last <- cur$pomp

  for (it in seq_len(Niter)) {
    prop <- theta*exp(stats::rnorm(length(theta),0,sds[names(theta)]))
    prop_lp <- log_prior(prop)
    if (is.finite(prop_lp)) {
      prop_fit <- ll_estimate(prop)
      logR <- (prop_fit$loglik+prop_lp) - (cur$loglik+cur_lp)
      accept <- is.finite(logR) && (log(stats::runif(1L)) < logR)
    } else {
      accept <- FALSE
      prop_fit <- list(loglik=-Inf,pomp=NULL)
    }
    if (accept) {
      theta <- prop; cur <- prop_fit; cur_lp <- prop_lp
      po_last <- cur$pomp
      naccept <- naccept+1L
    }
    chain[it,] <- c(theta,cur$loglik,cur_lp,as.numeric(accept))
    if (verbose && it %% max(1L,Niter %/% 20L) == 0L)
      message("mtbd2_pmcmc: iter ",it,"/",Niter,
        ", loglik=",signif(cur$loglik,5),
        ", accept rate=",signif(naccept/it,3))
  }

  list(
    chain=tibble::as_tibble(as.data.frame(chain)),
    acceptance_rate=naccept/Niter,
    pomp=po_last
  )
}
