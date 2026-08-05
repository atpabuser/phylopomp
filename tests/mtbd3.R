png(filename="mtbd3-%02d.png",res=100)

options(digits=3)
suppressPackageStartupMessages({
  library(ggplot2)
  library(phylopomp)
})
set.seed(20260724)

## Basic run/continue, with both sampling channels active:
## psi_i = overall per-capita sampling rate for type i,
## r_i   = probability a sampled type-i individual is removed
##         (r_i=1 -> sample_death only; r_i<1 -> some sample() events too).
x <- runMTBD3(
  time=5,
  I1_0=3,
  psi_1=0.2, psi_2=0.2, psi_3=0.2,
  r_1=0.5, r_2=0.5, r_3=0.5
)
stopifnot(
  attr(x, "model") == "MTBD3",
  inherits(x, "gpsim")
)
y <- continueMTBD3(x, time=6)
plot_grid(
  plot(x)+expand_limits(x=6),
  plot(y)+geom_vline(xintercept=5),
  y |>
    curtail(time=5) |>
    plot()+expand_limits(x=6),
  ncol=1
)

## Exercise cross-type births, migration, and both sampling channels
## via the by-name simulate() interface.
simulate(
  "MTBD3",
  time=3,
  lambda_1_1 = 1.5,
  lambda_1_2 = 0.4,
  lambda_2_1 = 0.4,
  lambda_2_2 = 1.5,
  m_1_2 = 0.2,
  m_2_1 = 0.2,
  mu_1 = 0.4, mu_2 = 0.4, mu_3 = 0.4,
  psi_1 = 0.2, psi_2 = 0.2, psi_3 = 0.2,
  r_1 = 0.5, r_2 = 0.5, r_3 = 0.5,
  I1_0 = 3, I2_0 = 1
) |> freeze(445178631) |>
  plot(prune=FALSE,obscure=FALSE)

y |> yaml() -> yaml_out
stopifnot(
  inherits(yaml_out, "gpyaml"),
  gregexpr("lambda",yaml_out) |> regmatches(yaml_out,m=_) |> lengths()==9,
  gregexpr("psi",yaml_out) |> regmatches(yaml_out,m=_) |> lengths()==3,
  gregexpr("r_",yaml_out) |> regmatches(yaml_out,m=_) |> lengths()==3
)

## Reduction to LBDP: zero out types 2 and 3 entirely (no births into
## them, no migration to them, no initial lineages), leaving a pure
## linear birth-death-sampling process in type 1. Rates are kept mild
## (near-critical) so the tree stays small. In lineages(obscure=FALSE),
## 'deme' is an integer code (1=type1, 2=type2, 3=type3); confirm types
## 2 and 3 never carry any lineages.
set.seed(4471)
g_mtbd3 <- runMTBD3(
  time=6,
  lambda_1_1=1, lambda_1_2=0, lambda_1_3=0,
  lambda_2_1=0, lambda_2_2=0, lambda_2_3=0,
  lambda_3_1=0, lambda_3_2=0, lambda_3_3=0,
  m_1_2=0, m_1_3=0, m_2_1=0, m_2_3=0, m_3_1=0, m_3_2=0,
  mu_1=0.3, mu_2=0, mu_3=0,
  psi_1=0.2, psi_2=0, psi_3=0,
  r_1=1, r_2=1, r_3=1,
  I1_0=3, I2_0=0, I3_0=0
)
ld_mtbd3 <- lineages(g_mtbd3,obscure=FALSE)
stopifnot(
  all(ld_mtbd3$lineages[ld_mtbd3$deme %in% c(2,3)]==0)
)

## With the same seed, MTBD3 restricted to type 1 (r_1=1, so all
## sampling is destructive) should draw an RNG stream identical to
## runLBDP with matched (lambda,mu,chi,n0): same nsample() exactly.
set.seed(2718)
a <- runMTBD3(
  time=6, I1_0=3, I2_0=0, I3_0=0,
  lambda_1_1=1, lambda_1_2=0, lambda_1_3=0,
  lambda_2_1=0, lambda_2_2=0, lambda_2_3=0,
  lambda_3_1=0, lambda_3_2=0, lambda_3_3=0,
  m_1_2=0, m_1_3=0, m_2_1=0, m_2_3=0, m_3_1=0, m_3_2=0,
  mu_1=0.3, mu_2=0, mu_3=0,
  psi_1=0.2, psi_2=0, psi_3=0,
  r_1=1, r_2=1, r_3=1
)
set.seed(2718)
b <- runLBDP(time=6, lambda=1, mu=0.3, chi=0.2, psi=0, n0=3)
stopifnot(
  identical(
    getInfo(a,nsample=TRUE)$nsample,
    getInfo(b,nsample=TRUE)$nsample
  )
)

dev.off()
