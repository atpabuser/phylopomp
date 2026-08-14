##' MTBD2 spillover-count posterior, pooled over tree/parameter uncertainty
##'
##' Constructs a posterior over the total number of realized cross-type
##' (\eqn{1\to2} and \eqn{2\to1}) birth events -- e.g. camel-to-human and
##' human-to-camel transmissions -- pooling particle-filter samples over a
##' supplied collection of (tree, parameter) pairs, both in the wider
##' (mostly unsampled) population and restricted to events ancestral to the
##' sampled genealogy.
##'
##' @name mtbd2_spillover
##' @rdname mtbd2_spillover
##' @include mtbd2_pomp.R
##' @param draws a list of draws, each a list with elements \code{x} (a genealogy in \pkg{phylopomp} format) and \code{params} (a named list of the arguments of \code{\link{mtbd2_pomp}}). See \code{\link{mtbd2_trajectory_posterior}}'s \code{draws} for the same convention; the two functions are designed to be run on the same \code{draws} list.
##' @param Np integer; number of particles per \code{\link[pomp]{pfilter}} call.
##' @param probs numeric vector of quantile levels (default the median and a 95\% interval, matching Vaughan & Stadler's reported HPDs).
##' @param verbose logical; report progress per draw.
##' @details
##' Each particle's \code{n12}/\code{n21} (total realized cross-type births,
##' anywhere in the population) and \code{n12a}/\code{n21a} (the subset
##' involving a currently tracked lineage, i.e. ancestral to the sampled
##' genealogy) accumulate over the whole simulated history -- see
##' \code{\link{mtbd2_pomp}}'s Details. Unlike
##' \code{\link{mtbd2_trajectory_posterior}}'s population-size trajectory,
##' a spillover count is a single terminal number per particle, not a time
##' series: this function takes each particle's \emph{final} saved value
##' (the running count as of the last observation) and pools those across
##' all particles of all supplied draws, exactly standing in for the ~400
##' thinned MCMC iterates Vaughan & Stadler pool over in their own
##' analysis. As there, \strong{the caller is responsible for ensuring
##' \code{draws} is a genuine (approximate) posterior sample} if the result
##' is to be interpreted as a Bayesian count posterior; a single repeated
##' (tree, parameter) draw is also supported and simply reports the
##' within-tree Monte Carlo count distribution at that one parameter value
##' (a cheap plug-in approximation that ignores parameter uncertainty --
##' see \code{\link{mtbd2_pmcmc}} for the properly Bayesian alternative).
##' @return
##' A \code{\link[tibble]{tibble}} with one row per count
##' (\code{n12}, \code{n21}, \code{n12a}, \code{n21a}), columns
##' \code{count} (the name), \code{ndraws}, \code{nparticles}, and one
##' column per requested quantile, named \code{q}\emph{p} as in
##' \code{\link{mtbd2_trajectory_posterior}}.
##' @importFrom pomp pfilter saved_states
##' @importFrom stats quantile
##' @importFrom tibble tibble
##' @export
mtbd2_spillover_posterior <- function (
  draws, Np = 2000,
  probs = c(0.025, 0.5, 0.975),
  verbose = FALSE
) {
  if (!is.list(draws) || length(draws) < 1L)
    pStop("mtbd2_spillover_posterior",": ",sQuote("draws")," must be a nonempty list.")

  countnames <- c("n12","n21","n12a","n21a")
  per_draw <- vector(mode="list",length=length(draws))

  for (k in seq_along(draws)) {
    d <- draws[[k]]
    if (is.null(d$x) || is.null(d$params))
      pStop("mtbd2_spillover_posterior",": draw ",k," must have elements ",
        sQuote("x")," and ",sQuote("params"),".")
    if (verbose) message("mtbd2_spillover_posterior: draw ",k,"/",length(draws))

    po <- do.call(mtbd2_pomp,c(list(x=d$x),d$params))
    pf <- pomp::pfilter(po,Np=Np,save.states="filter")
    ss <- pomp::saved_states(pf,format="data.frame")
    ss <- ss[ss$name %in% countnames,,drop=FALSE]
    ## each particle's FINAL saved value (the running count as of the last
    ## observation) -- ss$time is nondecreasing per particle by construction,
    ## so the max-time row per (.id,name) is the terminal count.
    last_time <- max(ss$time)
    ss <- ss[ss$time==last_time,,drop=FALSE]
    ss$draw <- k
    per_draw[[k]] <- ss[,c("name","value","draw",".id")]
  }

  all <- do.call(rbind,per_draw)

  qnames <- paste0("q",probs*100)
  rows <- vector(mode="list",length=length(countnames))
  for (i in seq_along(countnames)) {
    nm <- countnames[i]
    sub <- all[all$name==nm,,drop=FALSE]
    qs <- stats::quantile(sub$value,probs=probs,names=FALSE)
    row <- as.list(qs)
    names(row) <- qnames
    row$count <- nm
    row$ndraws <- length(unique(sub$draw))
    row$nparticles <- nrow(sub)
    rows[[i]] <- as.data.frame(row,stringsAsFactors=FALSE)
  }
  out <- do.call(rbind,rows)
  tibble::tibble(out[,c("count","ndraws","nparticles",qnames)])
}
