# [Integrated community models: A framework combining multi-species data sources to estimate the status, trends, and dynamics of biodiversity](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.14012)

### Elise F. Zipkin, Jeffrey W. Doser, Courtney L. Davis, Wendy Leuenberger, Samuel Ayebare, Kayla L. Davis

### Journal of Animal Ecology

### Code/Data DOI: [![DOI](https://zenodo.org/badge/586894609.svg)](https://zenodo.org/badge/latestdoi/586894609)

### Please contact the first author for questions about the code or data used in the manuscript: Elise F. Zipkin 

---------------------------------

## Abstract

1. Data deficiencies among rare or cryptic species preclude assessment of community-level processes using many existing approaches, limiting our understanding of the trends and stressors for large numbers of species, often within taxonomically related groups. Yet evaluating the dynamics of whole communities, not just common or charismatic species, is critical to understanding and predicting the responses of biodiversity to ongoing environmental pressures.   
2. A recent surge in both citizen science and government-funded data collection efforts has led to a wealth of biodiversity data. However, these data collection programs use a wide range of sampling protocols (from unstructured, opportunistic observations of wildlife to well-structured, design-based programs) and record information at a variety of spatiotemporal scales. As a result, available biodiversity data vary substantially in quantity and quality, which must be carefully reconciled for meaningful ecological analysis. 
3. Hierarchical modeling, including single-species integrated models and hierarchical community models, has improved our ability to assess and predict biodiversity trends and processes. Here, we highlight an ‘integrated community modeling’ framework that combines both data integration and community modeling to simultaneously improve inferences on species- and community-level dynamics. 
4. We illustrate the framework with a series of worked examples. Our three case studies demonstrate the utility of integrated community models in extending the geographic scope to evaluate species distributions and community-level richness patterns, discerning population and community trends over time, and estimating demographic rates and population growth for communities of sympatric species. We implemented these worked examples using multiple software methods through the R platform via packages with a formula-based interface and through development of custom code in JAGS, NIMBLE, and Stan. 
5. Integrated community models provide an exciting approach to model biological and observational processes for multiple species using multiple data types and sources simultaneously, thus accounting for uncertainty and sampling error within a unified framework. By leveraging the combined benefits of both data integration and community modeling, integrated community models can produce valuable information about both common and rare species as well as community-level dynamics, allowing us to evaluate holistically the effects of global change on biodiversity. 

## Repository Directory

This repository is split into three sub-directories, with each sub-directory corersponding to one of the case studies in the manuscript. Each sub-directory contains the necessary code and data to implement the integrated community model. See the `README` files in each sub-directory for further information.  

### [case-study-1-birds](./case-study-1-birds): contains all code, data, and documents related to Case Study 1: Spatial distributions of forest birds across the Northeastern United States.

### [case-study-2-butterfly](./case-study-2-butterfly): contains all code, data, and documents related to Case Study 2: Temporal trends of butterflies in the Midwestern United States.

### [case-study-3-mipm](./case-study-3-mipm): contains all code, data, and documents related to Case Study 3: Estimating species- and community-level demographcic rates and population growth.
