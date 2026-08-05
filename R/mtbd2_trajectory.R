##' MTBD2 population-trajectory posterior, pooled over tree/parameter uncertainty
##'
##' Constructs a posterior over the MTBD2 population trajectory
##' \eqn{(I_1(t), I_2(t))}, pooling particle-filter samples over a supplied
##' collection of (tree, parameter) pairs.
##'
##' @name mtbd2_trajectory
##' @rdname mtbd2_trajectory
##' @include mtbd2_pomp.R
##' @param draws a list of draws, each itself a list with elements \code{x} (a genealogy in \pkg{phylopomp} format) and \code{params} (a named list of the arguments of \code{\link{mtbd2_pomp}}, i.e. \code{lambda_1_1}, \code{lambda_1_2}, \code{lambda_2_1}, \code{lambda_2_2}, \code{m_1_2}, \code{m_2_1}, \code{mu_1}, \code{mu_2}, \code{psi_1}, \code{psi_2}, \code{r_1}, \code{r_2}, \code{I1_0}, \code{I2_0}, \code{pop}, and optionally \code{obs_tiptype}).
##' @param Np integer; number of particles per \code{\link[pomp]{pfilter}} call.
##' @param times numeric vector of times at which to summarize the pooled trajectory. If \code{NULL} (the default), each draw's own observation times are used and results are pooled onto their union without interpolation (see Details).
##' @param probs numeric vector of quantile levels for the summary band (default the median and a 95\% interval).
##' @param verbose logical; report progress per draw.
##' @details
##' This constructs the \pkg{phylopomp} analogue of the trajectory posterior described by Vaughan & Stadler (2025, \emph{Mol.\ Biol.\ Evol.}\ 42, msaf130), \emph{without} performing their full three-step pipeline. In their approach, \strong{Step 1} (MCMC via BDMM-Prime, sampling the tree topology, MTBD parameters, and nucleotide-model parameters jointly from a genomic alignment) is a wholly external inference not implemented in this package -- it requires BEAST2/BDMM-Prime or an equivalent tool applied to sequence data. \strong{Steps 2--3} (stochastic mapping of ancestral types onto tree edges, followed by a bespoke trajectory particle filter) collapse in \pkg{phylopomp} into a single \code{\link[pomp]{pfilter}} call per tree: because the exact King--Lin--Ionides filter jointly simulates the coloring \emph{and} the full population state at every step, \code{\link[pomp]{saved_states}} already returns per-particle \eqn{(I_1(t), I_2(t))} trajectory samples, conditional on that one tree and one parameter vector, with no separate stochastic-mapping stage required.
##'
##' \code{mtbd2_trajectory_posterior} supplies only the \emph{pooling} that Steps 2--3 would perform \emph{across} an externally supplied ensemble of tree/parameter draws, standing in for the ~400 thinned MCMC iterates Vaughan & Stadler use. \strong{The caller is responsible for ensuring \code{draws} is a genuine (approximate) posterior sample} -- e.g., from a BDMM-Prime MCMC run, thinned appropriately -- if the result is to be interpreted as a Bayesian trajectory posterior. Passing repeated draws of a single fixed tree and parameter vector (e.g., for a sensitivity check) is also supported and simply reproduces the single-tree case with within-tree Monte Carlo variability only.
##'
##' Pooling with \code{times = NULL} takes the union of each draw's own observation times and, for each output time, uses \emph{only} the draws whose own time grid includes it exactly (typically only the terminal time is shared across draws with distinct trees). For a pooled summary interpolated onto a common time grid instead, supply \code{times} explicitly; each particle's trajectory is interpolated onto that grid independently via right-continuous (last-observation-carried-forward) step interpolation, matching the piecewise-constant nature of the population trajectory between genealogy events, before pooling across particles and draws.
##' @return A \code{\link[tibble]{tibble}} with columns \code{time}, \code{deme} (\code{1} or \code{2}), \code{ndraws} (the number of draws contributing to that time/deme cell), and one column per requested quantile, named \code{q}\emph{p} for quantile level \emph{p} (e.g. \code{q2.5}, \code{q50}, \code{q97.5}).
##' @importFrom pomp pfilter saved_states
##' @importFrom stats quantile approx
##' @importFrom tibble tibble as_tibble
##' @export
mtbd2_trajectory_posterior <- function (
  draws, Np = 2000, times = NULL,
  probs = c(0.025, 0.5, 0.975),
  verbose = FALSE
) {
  if (!is.list(draws) || length(draws) < 1L)
    pStop("mtbd2_trajectory_posterior",": ",sQuote("draws")," must be a nonempty list.")

  per_draw <- vector(mode="list",length=length(draws))

  for (k in seq_along(draws)) {
    d <- draws[[k]]
    if (is.null(d$x) || is.null(d$params))
      pStop("mtbd2_trajectory_posterior",": draw ",k," must have elements ",
        sQuote("x")," and ",sQuote("params"),".")
    if (verbose) message("mtbd2_trajectory_posterior: draw ",k,"/",length(draws))

    po <- do.call(mtbd2_pomp,c(list(x=d$x),d$params))
    pf <- pomp::pfilter(po,Np=Np,save.states="filter")
    ss <- pomp::saved_states(pf,format="data.frame")
    ss <- ss[ss$name %in% c("I1","I2"),,drop=FALSE]
    ss$deme <- ifelse(ss$name=="I1",1L,2L)
    ss$draw <- k
    ## '.id' distinguishes particles within this draw's pfilter run; make
    ## it globally unique across draws so later grouping is unambiguous.
    ss$particle <- paste(k,ss$.id,sep="_")
    per_draw[[k]] <- ss[,c("time","deme","value","draw","particle")]
  }

  all <- do.call(rbind,per_draw)

  if (!is.null(times)) {
    ## Interpolate each individual PARTICLE's trajectory onto the common
    ## 'times' grid (not the pooled-across-particles subset -- a
    ## multi-root tree yields several rows at time 0 per particle for
    ## distinct roots, which would otherwise make x non-unique for
    ## approx()). Right-continuous step interpolation (last observed
    ## value carried forward) matches the piecewise-constant population
    ## process.
    particles <- unique(all[,c("draw","deme","particle")])
    interp_list <- vector(mode="list",length=nrow(particles))
    for (i in seq_len(nrow(particles))) {
      pk <- particles$particle[i]; dm <- particles$deme[i]; dr <- particles$draw[i]
      sub <- all[all$particle==pk & all$deme==dm,,drop=FALSE]
      sub <- sub[order(sub$time),]
      ## collapse any exact-duplicate times (multiple roots at t=0, or
      ## simultaneous events) by taking the last (latest-processed) value.
      sub <- sub[!duplicated(sub$time,fromLast=TRUE),]
      f1 <- stats::approx(sub$time,sub$value,xout=times,method="constant",
        rule=2,f=0)
      interp_list[[i]] <- data.frame(
        time=times,deme=dm,value=f1$y,draw=dr,particle=pk
      )
    }
    all <- do.call(rbind,interp_list)
  }

  ## quantile summary per (time,deme) cell, pooling across all particles
  ## from all contributing draws.
  key <- interaction(all$time,all$deme,drop=TRUE)
  groups <- split(seq_len(nrow(all)),key)
  qnames <- paste0("q",probs*100)
  rows <- vector(mode="list",length=length(groups))
  for (i in seq_along(groups)) {
    idx <- groups[[i]]
    sub <- all[idx,,drop=FALSE]
    qs <- stats::quantile(sub$value,probs=probs,names=FALSE)
    row <- as.list(qs)
    names(row) <- qnames
    row$time <- sub$time[1L]
    row$deme <- sub$deme[1L]
    row$ndraws <- length(unique(sub$draw))
    rows[[i]] <- as.data.frame(row,stringsAsFactors=FALSE)
  }
  out <- do.call(rbind,rows)
  out <- out[order(out$deme,out$time),c("time","deme","ndraws",qnames)]
  tibble::as_tibble(out)
}
