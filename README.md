# An extension of the Unified Skew-Normal family of distributions and application to Bayesian binary regression

This repository is linked to the article [Onorati and Liseo (2023). *An extension of the Unified Skew-Normal family of distributions and application to Bayesian binary regression*](https://arxiv.org/abs/2209.03474) and provides the code for the R package `pSUN`.

There are also the codes implementing the empirical examples of the article including the execution \*.Rout files. \*.RData were excluded for memory space limit.

Notice that the package `pSUN` requires the following dependencies:

* GIGrvg
* TruncatedNormal
* matrixStats
* nleqslv
* statmod
* truncnorm

Furthermore, the following packages are required to obtain the final results:

* sns
* arm
* nlme
* boot
* xtable
* matrixStats


### Simulation Study: Coverage Analysis

Go to `SimStudyIS`. Then there are 6 subfolders named `SimStudyIS_n\*p\*`; each of them has an R script named `SimStudyIS_n\*p\*.R`. Once all R scripts are executed, go back to `SimStudyIS` and run `SimStudyIS_getResults.R` in order to obtain the results. Do not touch the file `SimStudyIS_ExeScript.R`.


### Simulation Study: Polya-Gamma vs pSUN with Small Sample Size

Go to `SimStudyPGvs_pSUN`. Then there are 3 subfolder, namely `SimStudyPGvs_pSUN_50n500pGauss`, `SimStudyPGvs_pSUN_50n500pLaplacit`, and `SimStudyPGvs_pSUN_50n500pDirLapl` for the Gaussian prior, independent Laplace (Laplacit) prior, and the Dirichlet-Laplace prior respectively. Inside of all subfolders there is an R script named `SimStudyPGvs_pSUN_50n500p\*.R`; once all R scripts are executed, go back to `SimStudyPGvs_pSUN` and run `SimStudyPGvs_pSUN_50n500p_getResults.R` in order to obtain the results.


### Simulation Study: a Logit Model with Unbalanced Data

Go to `SimStudyLogitUnbalanced`. Then there are 3 subfolder, namely `SimStudyLogitUnbalanced_50n1pGauss`, `SimStudyLogitUnbalanced_200n1pGauss`, and `SimStudyLogitUnbalanced_1000n1pGauss` for the different sample sizes. Inside of all subfolders there is an R script named `SimStudyLogitUnbalanced_\*n1pGauss.R`; once all R scripts are executed, go back to `SimStudyLogitUnbalanced` and run `SimStudyLogitUnbalanced_getResults.R` in order to obtain the results.


### Simulation Study: a Probit Model with Sparse Parameters

Go to `SimStudyProbitSparse_50n500pGauss`. Then there are 3 subfolder, namely `SimStudyProbitSparse_50n500pGauss`, `SimStudyProbitSparse_50n500pLaplacit`, and `SimStudyProbitSparse_50n500pDirLapl` for the different priors. Inside of all subfolders there is an R script named `SimStudyProbitSparse_50n500p\*.R`; once all R scripts are executed, go back to `SimStudyProbitSparse_50n500p` and run `SimStudyProbitSparse_getResults.R` in order to obtain the results.


### A Real Dataset: Parameters Estimation for Cancer SAGE

Go to `CancerSAGE_ParEst` and execute the R script named `CancerSAGE_ParEst.R`. Run the R script `CancerSAGE_getResults.R` to obtain the results.


### A Real Dataset: Model and Variable Selection for Lung Cancer

Go to `LungCancer_ModSel` and execute the R script named `LungCancer_ModSel.R`. Run the R script `LungCancer_getResults.R` to obtain the results.
