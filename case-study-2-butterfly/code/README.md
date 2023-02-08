# code

This directory contains all code for raw data analysis, fitting models, and producing the results shown in the manuscript for the butterfly cases tudy. The scripts directly in the `code` directory are data prep scripts used to prepare the butterfly data for fitting the model in both `spAbundance` and `Stan`. The scripts in `spAbundance-code/` and `stan-code`/ contain the specific scripts for fitting the model in `spAbundance` and `Stan`, respectively.   

1. `butterfly-data-prep.R`: script to prepare the raw data from the five butterfly data sets into the format necessary for fitting the integrated community model.

## spAbundance-code

`spAbundance` is an R package that is currently only available via download from GitHub. We fit all models using `v0.1.0` of `spAbundance`. The current in development version of this package can be installed via `R` by using `devtools::install_github("doserjef/spAbundance")`.  See [the package GitHub page](https://github.com/doserjef/spAbundance) for more details. This package will be released on CRAN by the end of May 2023.

1. `main-spAbundance.R`: code to run the integrated community model using the [`spAbundance` R package](https://github.com/doserjef/spAbundance). 
2. `summary-spAbundance.R`: code to summarize results from the `spAbundance` model fit and generate Figure 3 in the manuscript. 


## stan-code

1. `main-stan.R`: code to call `Stan` code using the `rstan` package and fit the integrated community model.
2. `icm.stan`: `Stan` file containing the custom `Stan` code for the integrated community model to assess trends in the butterfly community. 
3. `summary-stan.R`: code to summarize results from the `stan` model fit and generate the corresponding Figure 3 in the manuscript, but using the `stan` model.
