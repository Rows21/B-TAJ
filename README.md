# A Bayesian Threshold-Aligned Joint Disease Progression Model (B-TAJ DPM) for Alzheimer’s Disease

<p align="center">
  <a href='https://arxiv.org/abs/2606.18139'>
        <img src='https://img.shields.io/badge/Paper-PDF-green?style=flat&logo=arXiv&logoColor=green' alt='arXiv PDF'>
      </a>
  <a href='https://creativecommons.org/licenses/by-nc/4.0/'>
      <img src='https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg'>
    </a>
  </p>

Official R implementation of B-TAJ DPM, from the following paper:

A Bayesian Threshold-Aligned Joint Disease Progression Model (B-TAJ DPM) for Alzheimer’s Disease.  \
Rong Wu, Duygu Tosun, Isabella Hausle, Margo Heston, and Aaron Wolfe Scheffler, \
University of California San Francisco \
for the Alzheimer’s Disease Neuroimaging Initiative

## Requirements

- R 4.6.1 (R 4.5.3 or later may work, but the current workflow was checked
  with R 4.6.1)
- Rtools for Windows, including `make` and `g++`, because NIMBLE compiles C++
  code
- Pandoc 2.8 or later to render an R Markdown file to HTML or PDF

Install the main R packages with:

```r
install.packages(c(
  "nimble", "matrixStats", "ggplot2", "splines2", "tidyr",
  "dplyr", "mgcv", "MASS", "knitr", "rmarkdown"
))
```

Both commands should return an executable path.

## Repository layout

- [`model/`](./model/) contains data generation, model code, custom samplers,
  and numerical helpers.
- [`simulation/`](./simulation/) contains simulation evaluation scripts.
- [`posterior code/`](./posterior%20code/) contains posterior extraction and
  posterior-predictive utilities, including
  [`posterior predictive distribution.R`](./posterior%20code/posterior%20predictive%20distribution.R).
- [`B-TAJ DPM Tutorial.Rmd`](./B-TAJ%20DPM%20Tutorial.Rmd) is the original
  simulation tutorial.


## Tutorial simulation

The tutorial loads the simulation components from [`model/`](./model/) and
uses `simu()` to:

1. generate longitudinal amyloid and tau observations;
2. fit the simulated data with the B-TAJ model;
3. retain posterior draws for trajectories, covariates, variance parameters,
   correlation, and subject-level random effects.

Open [`B-TAJ DPM Tutorial.Rmd`](./B-TAJ%20DPM%20Tutorial.Rmd) in RStudio or
run:

```r
setwd("E:/Github/B-TAJ")
knitr::knit("B-TAJ DPM Tutorial.Rmd")
```

### Training posterior

The Rmd first runs:

```r
res <- simu(
  varying_coef_code_sim,
  varying_coef_code,
  n_subjects = 1000,
  niter = 2000,
  nburnin = 1000
)
```

This is the training fit. Its posterior samples are kept fixed and are used
as the fitted population-level model.

### New data prediction

The Rmd then creates a new simulation with:

```r
data_1 <- generate_d1(n_subjects = 2, ...)
simudata <- generate_vc(varying_coef_code_sim, data_1, ...)
```

This N=2 dataset is not refit as a new population model. Instead, its
observed amyloid and tau data are passed to
[`posterior_predictions_code`](./posterior%20code/posterior%20predictive%20distribution.R).
The fixed training posterior supplies the population trajectory and variance
parameters, while the new data update the subject-specific `alpha` random
effects.

The Rmd then reports posterior threshold-age predictions. Amyloid predictions
are produced for both new subjects; tau predictions are produced only for
subjects included in `tau_amy_link` and therefore having tau observations.

The Rmd may take several minutes because it compiles NIMBLE models and runs
MCMC twice: once for the N=1000 training simulation and once for the N=2
posterior prediction.


## Citation

If you use this implementation, cite:


## License

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png

© 2025 Rong Wu.
