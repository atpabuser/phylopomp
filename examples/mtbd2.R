library(phylopomp)
set.seed(101)

## psi_i is the overall per-capita sampling rate for type i;
## r_i is the probability that a sampled type-i individual is removed
## (r_i=1: destructive/sample_death only; r_i<1: some samples are
## non-destructive/sample, leaving the lineage active).
runMTBD2(
  time=3,
  I1_0=3,
  psi_1=0.2, psi_2=0.2,
  r_1=0.5, r_2=0.5
) |>
  plot()
