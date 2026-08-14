##' Vaughan & Stadler (2025) MTBD2 reparameterization
##'
##' Converts between the \eqn{(R_{ij},\delta_i,s_i,r_i)} parameterization used
##' by Vaughan & Stadler (2025, \emph{MBE} 42(6), \doi{10.1093/molbev/msaf130})
##' for their BDMM-Prime analysis of the Dudas et al. (2018) MERS-CoV
##' camel/human dataset, and \pkg{phylopomp}'s native
##' \eqn{(\lambda_{ij},\mu_i,\psi_i,r_i)} rate parameters used by
##' \code{\link{mtbd2_pomp}}. Also provides their published priors (Table 2 of
##' that paper) as a sampler, for prior-predictive checks or as jumping-off
##' points for \code{\link{mtbd2_mif2}} or a particle-MCMC run.
##'
##' @name mtbd2_vs2025
##' @rdname mtbd2_vs2025
##' @include mtbd2_pomp.R
##' @details
##' Vaughan & Stadler parameterize a two-type (host-type) linear MTBD model by:
##' \itemize{
##'   \item \eqn{R_{ij}}: for \eqn{i=j}, the within-type effective reproductive
##'     number; for \eqn{i\neq j}, the cross-type (here, zoonotic)
##'     transmission analogue.
##'   \item \eqn{\delta_i = \mu_i+\psi_i}: the "become-uninfectious" rate for
##'     type \eqn{i}.
##'   \item \eqn{s_i = \psi_i/\delta_i}: the sampling proportion for type
##'     \eqn{i}.
##'   \item \eqn{r_i}: probability of removal on sampling (they fix
##'     \eqn{r_i=1} throughout; \code{\link{mtbd2_pomp}} allows \eqn{r_i<1}
##'     for sampled-ancestor models, but that regime is outside their paper's
##'     application).
##' }
##' The relationships \eqn{R_{ij}=\lambda_{ij}/\delta_i},
##' \eqn{s_i=\psi_i/\delta_i} invert to
##' \deqn{\lambda_{ij} = R_{ij}\,\delta_i, \qquad \psi_i = \delta_i\,s_i,
##' \qquad \mu_i = \delta_i\,(1-s_i).}
##' Migration (\code{m_1_2}, \code{m_2_1} in \code{\link{mtbd2_pomp}}) plays no
##' role in their parameterization or their MERS application (spillover is
##' carried entirely by the off-diagonal \eqn{R_{ij}}, i.e. by
##' \code{lambda_1_2}/\code{lambda_2_1}); \code{mtbd2_vs2025_params} therefore
##' always returns \code{m_1_2 = m_2_1 = 0}. Their model \emph{also} assumes
##' \eqn{s_i} is a two-piece step function in time, zero outside the interval
##' between the first and last sample of each type -- this reparameterization
##' does not itself implement that (see \code{\link{mtbd2_pomp}}'s
##' documentation and the package's known-limitations notes for the current
##' lack of covariate/skyline support); the \code{s_i} passed in or drawn here
##' is the \emph{within-window} value.
##' @param R_1_1,R_1_2,R_2_1,R_2_2 numeric; reproductive numbers (diagonal:
##'   within-type; off-diagonal: cross-type/zoonotic), i.e. \eqn{R_{ij}}
##'   entries of Table 2.
##' @param delta_1,delta_2 numeric; become-uninfectious rates \eqn{\delta_i}
##'   (same units as the tree's edge lengths, e.g. \eqn{\mathrm{yr}^{-1}} for
##'   a tree measured in years).
##' @param s_1,s_2 numeric in \eqn{[0,1)}; sampling proportions \eqn{s_i}
##'   (within the sampling window).
##' @param r_1,r_2 numeric in \eqn{[0,1]}; removal-on-sampling probabilities.
##'   Default \code{1} matches Vaughan & Stadler's assumption.
##' @return
##' \code{mtbd2_vs2025_params} returns a named list with elements
##' \code{lambda_1_1}, \code{lambda_1_2}, \code{lambda_2_1}, \code{lambda_2_2},
##' \code{m_1_2}, \code{m_2_1}, \code{mu_1}, \code{mu_2}, \code{psi_1},
##' \code{psi_2}, \code{r_1}, \code{r_2} -- i.e., exactly the rate arguments
##' expected by \code{\link{mtbd2_pomp}} (still missing \code{I1_0},
##' \code{I2_0}, \code{pop}, which are not part of their \eqn{(R,\delta,s)}
##' parameterization and must be supplied separately).
##' @export
mtbd2_vs2025_params <- function (
  R_1_1, R_1_2, R_2_1, R_2_2,
  delta_1, delta_2,
  s_1, s_2,
  r_1 = 1, r_2 = 1
) {
  if (any(c(delta_1,delta_2) <= 0))
    pStop("mtbd2_vs2025_params",": ",sQuote("delta_1")," and ",
      sQuote("delta_2")," must be positive.")
  if (any(c(s_1,s_2) < 0 | c(s_1,s_2) >= 1))
    pStop("mtbd2_vs2025_params",": ",sQuote("s_1")," and ",sQuote("s_2"),
      " must lie in [0,1).")
  list(
    lambda_1_1=R_1_1*delta_1, lambda_1_2=R_1_2*delta_1,
    lambda_2_1=R_2_1*delta_2, lambda_2_2=R_2_2*delta_2,
    m_1_2=0, m_2_1=0,
    mu_1=delta_1*(1-s_1), mu_2=delta_2*(1-s_2),
    psi_1=delta_1*s_1, psi_2=delta_2*s_2,
    r_1=r_1, r_2=r_2
  )
}

##' @rdname mtbd2_vs2025
##' @param n integer; number of draws.
##' @param d_1,d_2 numeric; duration in years between the first and last
##'   sample of each type in the dataset at hand (fixed data-derived
##'   constants, not estimated). Used as the lower bound of the bounded
##'   origin-time sampler described below (they take the max over types in
##'   their single-origin application; do the same here by passing
##'   \code{max(d_1,d_2)} if a single scalar is wanted).
##' @details
##' \code{mtbd2_vs2025_priors} draws \code{n} independent parameter sets from
##' the priors of Table 2: \eqn{R_{ii}\sim\mathrm{LogNormal}(0,0.5)},
##' \eqn{R_{ij}\sim\mathrm{Exp}(1)} for \eqn{i\neq j},
##' \eqn{\delta_i\sim\mathrm{LogNormal}(3,0.5)} (per year),
##' \eqn{s_i\sim\mathrm{Unif}(0,0.1)}, \eqn{r_i=1} fixed -- all four confirmed
##' directly against their published MERS-CoV BEAST 2 XML
##' (\code{mers_relaxed_part2.xml} in
##' \url{https://github.com/tgvaughan/MultiTypeTrajectoryAnalyses}), which
##' matches this sampler exactly for these parameters.
##'
##' The origin time \eqn{T} is the one exception: that same XML encodes its
##' prior as BEAST 2's \code{OneOnX}, an \emph{unbounded} Jeffreys-style
##' density \eqn{p(T)\propto 1/T} on \eqn{(0,\infty)} -- not the bounded
##' \eqn{\mathrm{LogUniform}(d_i,20)} this function previously claimed to
##' reproduce. \code{OneOnX} is improper (infinite total mass) and so cannot
##' itself be sampled from; if \code{d_1}/\code{d_2} are supplied, the
##' \code{T} column here instead draws from
##' \eqn{T\sim\mathrm{LogUniform}(\max(d_1,d_2),20)}, which has the same
##' log-uniform \emph{shape} as \code{OneOnX} but is a proper, boundedly
##' truncated stand-in for it -- adequate for prior-predictive checks or as
##' \code{\link{mtbd2_mif2}}/\code{\link{mtbd2_pmcmc}} starting values, but
##' not a literal reproduction of their encoded prior. Each row is also
##' passed through \code{mtbd2_vs2025_params} so the native
##' \eqn{(\lambda,\mu,\psi)} columns are included alongside the
##' \eqn{(R,\delta,s)} columns actually drawn.
##' @return
##' \code{mtbd2_vs2025_priors} returns a \code{\link[tibble]{tibble}} with
##' \code{n} rows, one column per \eqn{(R,\delta,s)} parameter drawn, one
##' column per native \eqn{(\lambda,\mu,\psi,r)} parameter, and (if
##' \code{d_1}/\code{d_2} given) a column \code{T} of origin-time draws.
##' @importFrom tibble tibble
##' @importFrom stats rlnorm rexp runif
##' @export
mtbd2_vs2025_priors <- function (n, d_1 = NULL, d_2 = NULL) {
  R_1_1 <- stats::rlnorm(n,0,0.5)
  R_2_2 <- stats::rlnorm(n,0,0.5)
  R_1_2 <- stats::rexp(n,1)
  R_2_1 <- stats::rexp(n,1)
  delta_1 <- stats::rlnorm(n,3,0.5)
  delta_2 <- stats::rlnorm(n,3,0.5)
  s_1 <- stats::runif(n,0,0.1)
  s_2 <- stats::runif(n,0,0.1)
  native <- mapply(
    mtbd2_vs2025_params,
    R_1_1=R_1_1,R_1_2=R_1_2,R_2_1=R_2_1,R_2_2=R_2_2,
    delta_1=delta_1,delta_2=delta_2,s_1=s_1,s_2=s_2,
    MoreArgs=list(r_1=1,r_2=1),
    SIMPLIFY=FALSE
  )
  out <- tibble::tibble(
    R_1_1=R_1_1,R_1_2=R_1_2,R_2_1=R_2_1,R_2_2=R_2_2,
    delta_1=delta_1,delta_2=delta_2,s_1=s_1,s_2=s_2,
    lambda_1_1=vapply(native,`[[`,numeric(1L),"lambda_1_1"),
    lambda_1_2=vapply(native,`[[`,numeric(1L),"lambda_1_2"),
    lambda_2_1=vapply(native,`[[`,numeric(1L),"lambda_2_1"),
    lambda_2_2=vapply(native,`[[`,numeric(1L),"lambda_2_2"),
    mu_1=vapply(native,`[[`,numeric(1L),"mu_1"),
    mu_2=vapply(native,`[[`,numeric(1L),"mu_2"),
    psi_1=vapply(native,`[[`,numeric(1L),"psi_1"),
    psi_2=vapply(native,`[[`,numeric(1L),"psi_2"),
    r_1=1,r_2=1
  )
  if (!is.null(d_1) || !is.null(d_2)) {
    lo <- max(c(d_1,d_2))
    out$T <- exp(stats::runif(n,log(lo),log(20)))
  }
  out
}

##' @rdname mtbd2_vs2025
##' @param origin numeric vector; candidate process-origin times, each
##'   expressed as a duration \emph{before the most recent sample} in the
##'   genealogy \code{x} (i.e. Vaughan & Stadler's \eqn{T}, in the same time
##'   units as \code{x}'s branch lengths). Each value must exceed the
##'   root-to-most-recent-sample duration of \code{x}, since a single seed
##'   lineage cannot reach the observed root any faster than the tree itself
##'   requires.
##' @details
##' \code{mtbd2_pomp}'s \code{t0} (the process start) is fixed, at
##' construction time, to the root time of \code{x} -- there is no
##' \eqn{T}-like free origin parameter, unlike Vaughan & Stadler's
##' \eqn{T\sim\mathrm{LogUniform}(d_i,20)}. \code{mtbd2_origin_slice}
##' profiles over \code{origin} by rebuilding the \sQuote{pomp} object at
##' each grid point with \code{t0} moved back by that amount (via
##' \code{\link[pomp]{timezero}<-}) and \code{I1_0}/\code{I2_0}
##' lineages grafted there instead of at the root; the intervening,
##' entirely-unobserved period is simulated by the same regular (non-
##' singular) Gillespie loop used everywhere else in the filter, so no
##' C-level change is needed to support this. This gives a log-likelihood
##' \emph{slice} over \eqn{T} (all other parameters held at \code{params}),
##' directly comparable in shape (though not, without many replicate
##' \code{\link[pomp]{mif2}} profiles at each grid point, in scale) to
##' Vaughan & Stadler's marginal posterior for \eqn{T}.
##' @return
##' \code{mtbd2_origin_slice} returns a \code{\link[tibble]{tibble}} with
##' columns \code{origin}, \code{loglik}, \code{se}, and \code{ess}, one row
##' per element of \code{origin}.
##' @importFrom pomp timezero<- timezero
##' @export
mtbd2_origin_slice <- function (x, params, origin, Np = 2000, nreps = 5, verbose = FALSE) {
  gi <- gendat(x)
  root_time <- gi$nodetime[1L]
  last_time <- gi$nodetime[length(gi$nodetime)]
  span <- last_time - root_time
  bad <- origin[origin <= span]
  if (length(bad) > 0L)
    pStop("mtbd2_origin_slice",": every element of ",sQuote("origin"),
      " must exceed the root-to-most-recent-sample duration (",
      signif(span,4),"); offending: ",paste(signif(bad,4),collapse=","),".")
  po0 <- do.call(mtbd2_pomp,c(list(x=x),params))
  res <- matrix(NA_real_,nrow=length(origin),ncol=3L,
    dimnames=list(NULL,c("loglik","se","ess")))
  for (i in seq_along(origin)) {
    if (verbose) message("mtbd2_origin_slice: origin=",origin[i],
      " (",i,"/",length(origin),")")
    po <- po0
    pomp::timezero(po) <- last_time - origin[i]
    lls <- replicate(nreps,pomp::logLik(pomp::pfilter(po,Np=Np)))
    if (all(!is.finite(lls))) {
      res[i,] <- c(-Inf,NA_real_,0)
    } else {
      est <- pomp::logmeanexp(lls[is.finite(lls)],se=TRUE,ess=TRUE)
      res[i,] <- c(unname(est[1L]),unname(est[2L]),unname(est[3L]))
    }
  }
  tibble::tibble(
    origin=origin,
    loglik=res[,"loglik"],
    se=res[,"se"],
    ess=res[,"ess"]
  )
}
