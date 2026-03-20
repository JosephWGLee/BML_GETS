# Bayesian Machine Learning approach to modelling gene expression time series 

This repository contains parametric and non-parametric bayesian approach to identifying gene expression profiles which change significantly over time across different experimental groups (e.g., Male, Female, Fasted, Fed) in response to fasting and fed conditions.

# Comparison Overview

Here, we compared and optimised two time-series differential expression tools:

- GPcounts (Non-parametric/Bayesain): Utilises a negative-binomial Gaussian Process framework. Utilizes a negative-binomial Gaussian Process framework. It is implemented in Python using TensorFlow and GPflow.
- maSigPro (Parametric/Frequentist): Applies polynomial regression and stepwise variable selection. It is implemented in R and is commonly used in the wet-lab environment. 

By comparing these tools, this project assesses how each handles the high levels of noise and overdispersion typical of stress-condition RNA-seq datasets and shows differences in the genes identified by each method.

#Contents

 - data: Example toy data (processed and preprocessed)
 - docs: Poster, presentation and full write up
 - notebooks
 - scripts
 - outputs: Example outputs with full data
 - requirements: Requirements for GPcounts python install (Python 3.10, GPflow 2.9, Tensorflow 2.18.0)

#Documentation

Further details regarding the experimental conditions, stress models, and biological interpretation can be found in the project report. Additional information for each package can be found in their respective repositories:

https://github.com/ManchesterBioinference/GPcounts
https://github.com/mjnueda/maSigPro
