# code

This directory contains all code for raw data analysis, fitting models, and producing the results shown in the manuscript for the butterfly cases tudy. The scripts directly in the `code` directory are data prep scripts used to prepare the butterfly data for fitting the model in both `spAbundance` and `Stan`. The scripts in `spAbundance-code/` and `stan-code`/ contain the specific scripts for fitting the model in `spAbundance` and `Stan`, respectively. Note that the `stan-code` directory is currently empty, but scripts for fitting this model in `Stan` will be added to this directory shortly.

1. `butterfly-data-prep.R`: script to prepare the raw data from the five butterfly data sets into the format necessary for fitting the integrated community model.


## spAbundance-code

`spAbundance` is an R package that is currently only available via download from GitHub. We fit all models using `v0.1.0` of `spAbundance`. The current in development version of this package can be installed via `R` by using `devtools::install_github("doserjef/spAbundance")`.  See [the package GitHub page](https://github.com/doserjef/spAbundance) for more details. This package will be released on CRAN by the end of May 2023.

1. `main-spAbundance.R`: code to run the integrated community model using the [`spAbundance` R package](https://github.com/doserjef/spAbundance). 
2. `summary-spAbundance.R`: code to summarize results from the `spAbundance` model fit and generate Figure 3 in the manuscript. 


## stan-code

This directory is currently empty, but will eventually contain scripts for fitting the integrated community model using custom Stan code.
