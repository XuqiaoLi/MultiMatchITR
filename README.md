# MultiMatch ITR
The R codes for our paper ***Multiple Treatments Matched Learning for Individualized Treatment Rules***. In our paper, a method based on matching is proposed to
estimate the optimal ITRs with continuous and survival data in multi-arm observational study. 

## Simulation studies
* `simulation_continuous.R` :  scenarios with continuous data for the proposed and compared methods.
* `simulation_survival.R` :  scenarios with survival data for the proposed and compared methods.

## Real data application
* `hcc-simplified.csv` : synthetic hepatocellular carcinoma data with noise.
* `ITR-covgpsmatch-funcs.R` : the functions for implementation, including matching with covariates and generalized propensity scores.
* `hcc-covgpsmatch.R` : the main function.
