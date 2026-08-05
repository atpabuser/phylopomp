library(phylopomp)

## Simulate a handful of "posterior draws": each with its own tree AND
## its own parameter jitter, standing in for thinned MCMC iterates from
## an external tree/parameter inference (e.g. BDMM-Prime).
set.seed(101)
theta0 <- list(
  lambda_1_1=1.0, lambda_1_2=0.3, lambda_2_1=0.2, lambda_2_2=0.9,
  m_1_2=0, m_2_1=0, mu_1=0.3, mu_2=0.3, psi_1=0.3, psi_2=0.3,
  r_1=1, r_2=1, I1_0=3, I2_0=2, pop=5
)
draws <- lapply(1:5, function (k) {
  th <- theta0
  th$lambda_1_2 <- th$lambda_1_2 * exp(rnorm(1,0,0.2))
  g <- do.call(
    runMTBD2,
    c(list(time=3), th[c("lambda_1_1","lambda_1_2","lambda_2_1","lambda_2_2",
      "m_1_2","m_2_1","mu_1","mu_2","psi_1","psi_2","r_1","r_2","I1_0","I2_0")])
  )
  list(x=g, params=th)
})

## Pooled population-trajectory posterior across all draws.
post <- mtbd2_trajectory_posterior(draws, Np=1000, times=seq(0,3,by=0.5))
print(post)
