##' @name bdss_pomp
##' @rdname bdss
##' @include bdss.R
##' @param x genealogy in \pkg{phylopomp} format.
##' @return
##' \code{bdss_pomp} returns a \sQuote{pomp} object.
##' @details
##' \code{bdss_pomp} constructs a \sQuote{pomp} object containing a given set of data and a BDSS model.
##' @importFrom pomp pomp onestep
##' @export
bdss_pomp <- function (
  x,
  lambda_nn, lambda_ns, lambda_sn, lambda_ss,
  mu, chi,
  N0, S0, pop
)
{
  x |> gendat() -> gi
  ivps <- structure(c(N0,S0),names=c("N0","S0"))
  if (any(ivps < 0))
    pStop(paste(sQuote(names(ivps)),collapse=","),
      " must be nonnegative.")
  pomp(
    data=NULL,
    t0=gi$nodetime[1L],
    times=gi$nodetime[-1L],
    params=c(
      lambda_nn=lambda_nn,lambda_ns=lambda_ns,
      lambda_sn=lambda_sn,lambda_ss=lambda_ss,
      mu=mu,chi=chi,
      pop=pop,ivps
    ),
    userdata=gi,
    nstatevars=6L + gi$nsample,
    rinit="bdss_rinit",
    rprocess=onestep("bdss_gill"),
    dmeasure="bdss_dmeas",
    statenames=c(
      "N","S","ll","node","ell1","ell2","color"
    ),
    paramnames=c(
      "lambda_nn","lambda_ns","lambda_sn","lambda_ss",
      "mu","chi","pop","N0","S0"
    ),
    PACKAGE="phylopomp"
  )
}
