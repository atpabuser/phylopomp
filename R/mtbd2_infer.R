##' Likelihood evaluation and exploration for the MTBD2 model
##'
##' Replicated log-likelihood evaluation, and log-likelihood curves in one
##' parameter, for a fixed MTBD2 genealogy.
##'
##' @name mtbd2_infer
##' @rdname mtbd2_infer
##' @include mtbd2_pomp.R
##' @param x genealogy in \pkg{phylopomp} format.
##' @param params named list of the arguments of \code{\link{mtbd2_pomp}} other than \code{x} (i.e. \code{lambda_1_1}, \code{lambda_1_2}, \code{lambda_2_1}, \code{lambda_2_2}, \code{m_1_2}, \code{m_2_1}, \code{mu_1}, \code{mu_2}, \code{psi_1}, \code{psi_2}, \code{r_1}, \code{r_2}, \code{I1_0}, \code{I2_0}, \code{pop}, and optionally \code{obs_tiptype}).
##' @param Np integer; number of particles per \code{\link[pomp]{pfilter}} call.
##' @param nreps integer; number of independent replicate filter runs to combine.
##' @details
##' The sequential Monte Carlo filter underlying \code{\link{mtbd2_pomp}} returns an \emph{unbiased} estimate of the likelihood \eqn{L}, not of \eqn{\log L}; by Jensen's inequality the log of a single run is therefore biased \emph{downward}. Replicate runs must accordingly be combined on the likelihood scale, which \code{mtbd2_loglik} does via \code{\link[pomp]{logmeanexp}} rather than by averaging log-likelihoods. The reported standard error is \code{logmeanexp}'s jackknife estimate, and \code{ess} is the effective sample size of the replicate ensemble (\emph{not} the particle-filter ESS within any one run -- see \code{\link[pomp]{eff_sample_size}} for that).
##'
##' A small \code{ess} relative to \code{nreps}, or a standard error large compared to differences between parameter values of interest, indicates that \code{Np} is too small for the tree and parameters at hand: the estimate is then dominated by a few high-weight replicates and should not be trusted. This is especially likely with \code{obs_tiptype=TRUE} on trees with many samples (see \code{\link{mtbd2_pomp}}).
##' @return
##' \code{mtbd2_loglik} returns a named numeric vector with elements \code{loglik}, \code{se}, and \code{ess}.
##' @importFrom pomp pfilter logmeanexp
##' @export
mtbd2_loglik <- function (x, params, Np = 2000, nreps = 5) {
  po <- do.call(mtbd2_pomp,c(list(x=x),params))
  lls <- replicate(nreps,pomp::logLik(pomp::pfilter(po,Np=Np)))
  if (all(!is.finite(lls)))
    return(c(loglik=-Inf,se=NA_real_,ess=0))
  est <- pomp::logmeanexp(lls[is.finite(lls)],se=TRUE,ess=TRUE)
  c(loglik=unname(est[1L]),se=unname(est[2L]),ess=unname(est[3L]))
}

##' @rdname mtbd2_infer
##' @param par character; the name of the single parameter to vary. Must be one of the names of \code{params}.
##' @param values numeric vector of values of \code{par} at which to evaluate the log-likelihood.
##' @param verbose logical; report progress.
##' @details
##' \code{mtbd2_slice} varies \code{par} across \code{values} while holding \emph{every other} parameter fixed at its value in \code{params}. This is a log-likelihood \emph{slice}, not a profile likelihood: a profile would re-maximize over all the other parameters at each grid point, which this does not do. A slice is much cheaper and is the right diagnostic for checking that the likelihood surface peaks near a known truth in a simulation study; it will, however, generally give narrower apparent intervals than a true profile when parameters are correlated, and should not be read as a confidence interval.
##' @return
##' \code{mtbd2_slice} returns a \code{\link[tibble]{tibble}} with columns \code{value}, \code{loglik}, \code{se}, and \code{ess}, one row per element of \code{values}.
##' @importFrom tibble tibble
##' @export
mtbd2_slice <- function (
  x, params, par, values, Np = 2000, nreps = 5, verbose = FALSE
) {
  if (!is.character(par) || length(par) != 1L)
    pStop("mtbd2_slice",": ",sQuote("par")," must be a single parameter name.")
  if (!(par %in% names(params)))
    pStop("mtbd2_slice",": ",sQuote(par)," is not among the names of ",
      sQuote("params"),".")
  res <- matrix(NA_real_,nrow=length(values),ncol=3L,
    dimnames=list(NULL,c("loglik","se","ess")))
  for (i in seq_along(values)) {
    if (verbose) message("mtbd2_slice: ",par,"=",values[i],
      " (",i,"/",length(values),")")
    p <- params
    p[[par]] <- values[i]
    res[i,] <- mtbd2_loglik(x,p,Np=Np,nreps=nreps)
  }
  tibble::tibble(
    value=values,
    loglik=res[,"loglik"],
    se=res[,"se"],
    ess=res[,"ess"]
  )
}

## Map the user-facing argument names of mtbd2_pomp() onto the internal
## parameter names carried by the 'pomp' object it builds.
mtbd2_parmap <- c(
  lambda_1_1="lambda_11", lambda_1_2="lambda_12",
  lambda_2_1="lambda_21", lambda_2_2="lambda_22",
  m_1_2="m12", m_2_1="m21",
  mu_1="mu1", mu_2="mu2",
  psi_1="psi1", psi_2="psi2",
  r_1="r1", r_2="r2",
  I1_0="I1_0", I2_0="I2_0", pop="pop"
)

##' @rdname mtbd2_infer
##' @param start named list of starting parameter values, in the same form as \code{params}.
##' @param estimate character vector of parameter names (user-facing names, as in \code{start}) to be estimated. All must be strictly positive in \code{start}: estimation proceeds on the log scale, so a parameter fixed at zero cannot be estimated (and a parameter that should be allowed to reach zero is better handled by comparing models with and without it).
##' @param rw_sd numeric; the random-walk standard deviation, on the \emph{log} scale, used to perturb the estimated parameters. Either a single value applied to all of \code{estimate}, or a named vector.
##' @param Nmif integer; number of iterated-filtering iterations.
##' @param cooling_fraction_50 numeric; the fraction by which the random-walk standard deviation is reduced after 50 iterations (see \code{\link[pomp]{mif2}}).
##' @details
##' \code{mtbd2_mif2} performs maximum-likelihood estimation by iterated filtering (\code{\link[pomp]{mif2}}), the standard approach for POMP models whose likelihood is available only via a stochastic (particle-filter) estimate, as is the case here. Parameters are estimated on the log scale, so they remain positive; parameters not named in \code{estimate} are held fixed at their values in \code{start}.
##'
##' Iterated filtering returns a \emph{point estimate}, not a posterior. The log-likelihood reported by \code{\link[pomp]{logLik}} on the returned object is the filtering estimate from the final iteration and is both noisy and biased; it should not be used for inference directly. Re-evaluate the likelihood at the returned point estimate with \code{mtbd2_loglik} (which replicates and combines correctly) before reporting it, and check convergence with \code{\link[pomp]{traces}} across several independent starts before trusting the estimate.
##' @return
##' \code{mtbd2_mif2} returns the \sQuote{mif2d_pomp} object produced by \code{\link[pomp]{mif2}}. Use \code{\link[pomp]{coef}} to extract the point estimate (in internal parameter names; see \code{\link{mtbd2_pomp}}), \code{\link[pomp]{traces}} to inspect convergence.
##' @importFrom pomp mif2 pomp parameter_trans rw_sd coef
##' @export
mtbd2_mif2 <- function (
  x, start, estimate, rw_sd = 0.02,
  Nmif = 50, Np = 2000, cooling_fraction_50 = 0.5,
  verbose = FALSE
) {
  bad <- setdiff(estimate,names(mtbd2_parmap))
  if (length(bad) > 0L)
    pStop("mtbd2_mif2",": not estimable parameter(s): ",
      paste(sQuote(bad),collapse=","),".")
  missing <- setdiff(estimate,names(start))
  if (length(missing) > 0L)
    pStop("mtbd2_mif2",": ",paste(sQuote(missing),collapse=",")," absent from ",
      sQuote("start"),".")
  vals <- unlist(start[estimate])
  if (any(vals <= 0))
    pStop("mtbd2_mif2",": estimation is on the log scale, so all of ",
      sQuote("estimate")," must be strictly positive in ",sQuote("start"),
      "; offending: ",
      paste(sQuote(estimate[vals<=0]),collapse=","),".")

  internal <- unname(mtbd2_parmap[estimate])

  po <- do.call(mtbd2_pomp,c(list(x=x),start))
  po <- pomp::pomp(po,
    partrans=pomp::parameter_trans(log=internal),
    paramnames=internal
  )

  if (length(rw_sd) == 1L && is.null(names(rw_sd))) {
    sds <- stats::setNames(rep(rw_sd,length(internal)),internal)
  } else {
    miss <- setdiff(estimate,names(rw_sd))
    if (length(miss) > 0L)
      pStop("mtbd2_mif2",": ",sQuote("rw_sd")," lacks entries for ",
        paste(sQuote(miss),collapse=","),".")
    sds <- stats::setNames(unname(rw_sd[estimate]),internal)
  }
  rwargs <- lapply(sds,identity)

  pomp::mif2(
    po,
    Nmif=Nmif,
    Np=Np,
    cooling.fraction.50=cooling_fraction_50,
    rw.sd=do.call(pomp::rw_sd,rwargs),
    verbose=verbose
  )
}
