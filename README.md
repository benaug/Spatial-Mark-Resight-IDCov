# Spatial-Mark-Resight-IDCov
Spatial Mark Resight samplers in nimble that allows for individual ID covariates (SMR + categorical SCR)

This code contains contributions from Glenn Stauffer (though any problems should be blamed on me).

There are 4 testscripts:
1) IDcov Poisson: Poisson observation model, ID covs only provide ID exclusions for unidentified samples.
2) IDcov DF Poisson: ID covs also used as detection function covariates
3) IDcov 2gamma Poisson: separate population proportions for ID covariates between marked and unmarked individuals (e.g. sex-bias in marking process)
4) IDcov DF 2gamma Poisson: 2 + 3

May add negative binomial versions later. Just requires swapping out the observation model, but needs to be done in custom functions as well. You can do generalized SMR by simply appending an explicit marking process.

Disclaimer: Code has been tested thoroughly, but I restructured it to better match the other latent ID SCR samplers on my github. I don't think I screwed anything up in the process, but will retest soon.

Disclaimer 2: This model and data structure is a bit tedious.

Disclaimer 3: This model assumes ID covs are missing at random. I think there is some robustness to violating this assumption, and greater robustness as the number/proportion of marked individuals increases and when fewer ID covs are missing. It also assumes each ID cov is independent of the others, which may also not be true. The ID covs can be structured in a dependent manner with more model modification (e.g., differing sex ratio between adults and subadults).

Disclaimer 4: The base model assumes the marked and unmarked populations have the same population frequencies for the ID covs (gamma). This may not be true, e.g., if you disproportionately mark males or adults. The "2gamma" testscript shows how to relax this assumption. Though you now must estimate the unmarked individual gammas with 0 identified individuals. This can obviously be done, at least with "good" data, because it is what is done in categorical SCR.

Disclaimer 5: The samplers in this repo let you do quite a few things with ID covariates, but what you will actually be able to do depends on your data. I recommend doing a quick simulation evaluation to see if what you want to do is possible with the data you have or expect to collect.