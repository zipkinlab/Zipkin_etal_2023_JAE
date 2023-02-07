# Case Study 2: Temporal trends of butterflies in the Midwestern United States

This repository contains all code related to case study 2 in which we develop an integrated community model to quantify relative abundance trends in ten open-habit-associated butterfly species in the Midwestern U.S. from 2008-2017 using data from five volunteer-based monitoring programs. The directory is split into the following structure. Each sub-directory contains its own README file with details on each of the files in that specific sub-directory. 

+ `code`: contains all scripts for data processing, running the integrated community model in `spAbundance`, and processing of the results. Also contains custom code for running the integrated community model in `Stan`.  
+ `data`: contains data for fitting the model. Note this only contains data from three of the five data sources, as two of the data sources (North American Butterfly Association (NABA) and Ohio Lepidopterists) are propietary. For access to the NABA data, please contact NABA directly (www.naba.org). 
+ `results`: directory containing model results files used to generate Figure 3 in the manuscript.
