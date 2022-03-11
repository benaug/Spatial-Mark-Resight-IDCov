# Spatial-Mark-Resight-IDCov
Spatial Mark Resight samplers in nimble that allows for individual ID covariates (SMR + categorical SCR)

This code contains contributions from Glenn Stauffer (though any problems should be blamed on me).

Disclaimer: Code has been tested thoroughly, but I restructured it to better match the other latent ID SCR samplers on my github. I don't think I screwed anything up in the process, but will retest soon.

More testscripts to come:
1. Show how to model lam0, sigma, or thinning rates (theta.marked, theta.unmarked) as function of ID covs (G.true).
2. Allow marked and unmarked inds to have different population cateogory level frequencies (gamma)
3. Negative Binomial observation model

Disclaimer 2: This model and data structure is a bit tedious.