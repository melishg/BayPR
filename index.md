---
title: BayPR: Empirical Bayes Inference for Prevalence Ratios
theme: jekyll-theme-cayman
filename: index.md
---

{% for page in site.pages %}
    <a href={{ page.filename }}>{{ page.title }}</a>
{% endfor %}

# BayPR: Empirical Bayes Inference for Prevalence Ratios

This package estimates the prevalence ratio of an exposure in a population of cases and a reference population using samples from each population. This can be applied to assess the disease-causing liability of DNA variants, as described in Collaco et al, 2021.




## Installation & Requirements

All models are fit using the `optimx` function in the [optimx](https://cran.r-project.org/web/packages/optimx/index.html) package. Dependencies can be installed using `install.packages` or using the 'Packages' panel in Rstudio. Any additional dependencies are only used to arrange, plot, and report model output.

Download `source-bayes_prevalence_ratio.r` from Github, and use the `source` function to read its functions into the R workspace.


### Files Included:

  - R Source code:
    - `source-bayes_prevalence_ratio.r`: R code for estimating Empirical Bayes Beta Binomial Models
    - `simulate_beta_binomial_data.r`: R code for simulating data according to a single-component or two-component beta binomial model.
  - Example Datasets:
    - `single_bb_example.csv`: an example dataset generated from a single-component beta binomial model
    - `mixture_bb_example.csv`: an example dataset generated from a two-component beta binomial model
  - R Markdown Report: Note that both .Rmd files are required.
    - `Main-bayes_prevalence_ratio.Rmd`
    - `report_body.Rmd`




## Usage:

`source-bayes_prevalence_ratio.r` contains three main functions:

  - `bb.prevalence.ratio.mm`: Estimate a single component beta-binomial model using the Method of Moments (MM)
  - `bb.prevalence.ratio.mml`: Estimate a single component beta-binomial model using the Maximum Marginal Likelihood (MML).
  - `bb.prevalence.ratio.mml`: Estimate a two-component beta-binomial mixture model using the method of Maximum Marginal Likelihood (MML)
  
These are the functions which fit each model: these only rely on the `optimx` package for model fitting. These return the estimated model parameters as well as posterior summaries for the prevalence ratio for each variant.


It is strongly recommended to fit the single and two-component beta binomial models using maximum marginal likelihood using `bb.prevalence.ratio.mml` and `bb.prevalence.ratio.mml` before attempting to use the R Markdown report. Using `verbose = TRUE` will return the optimization object from `optimx`, allowing for more detailed inspection of results.

After fitting the models using these functions, the results of all three models can be compared by using the `Main-bayes_prevalence_ratio.Rmd`: this report runs the models, saves a spreadsheet of results, and attempts to plot the results. Plots may not be well behaved in the case of numerical overflow (values of `Inf` when the result is larger than what can be represented).

Two simulated datasets are provided in Github: users can create additional simulated datasets. Users can also assess the effect of varying sample size ratios by modification of the simulation code.




## Input file:

The input file should be saved in the directory in the file path specified in line 51 of Main-bayes_prevalence_ratio.Rmd
The input file should be a .csv file with 6 columns: n.1, y.1, n.0, y.0, event.label, flag

  - n.1 (integer): the total number of exposed among cases
    - e.g. Alleles genotyped among cases 
  - y.1 (integer): the total number of events among cases
    - e.g. Allele count among cases
  - n.0 (integer): the total number of exposed among reference sample
    - e.g. Alleles genotyped among reference sample 
  - y.0 (integer): the total number of events among reference sample
    - e.g. Allele count among reference sample    
  - event.label (character string): a uniquely identified variant ID, containing only letters, numbers, periods, and underscores. Used for plotting and arranging results. Long values may distort labeling of plots. 
  - flag ("exclude" or "include"): A variable indicating whether to include or exclude a variant while running the script.

Output file:
- The output file will be saved to the directory in the file path specified in line 54 of Main-bayes_prevalence_ratio.Rmd
- The file name will be as indicated in line 60 of Main-bayes_prevalence_ratio.Rmd




## Output

The report will create a .csv file for results. The name of this file will be based on the name of the `dataset.name` parameter.

  - id: (from input file) a uniquely identified variant ID
  - y.0: (from input file) the allele count among controls 
  - n.0: (from input file) the total number of alleles genotyped among controls 
  - y.1: (from input file) the allele count among cases
  - n.1: (from input file) the total number of alleles genotyped among cases 
  - t.k: The total number of observed variants 
  - r.k: The relative sample size in cases relative to controls 
  - theta.k.hat: Transformed prevalence ratio (see documentation)
  - pi.0: The sample prevalence in Reference Sample
  - pi.1: The sample prevalence in the Case sample
  - gamma.hat: Direct Estimate (ratio of sample prevalence values)
  - ll.1: Log Likelihood - Component 1
  - ll.2: Log Likelihood - Component 2
  - log.mlr.2.1: Log Marginal Likelihood Ratio of component 2 vs. 1
  - mlr.2.1: Marginal likelihood ratio of component 2 vs. 1
  - posterior.odds.2.1: Posterior odds of component 2 vs. 1
  - posterior.probability.2: Posterior probability of component 2
  - mean_of_posterior_medians: internal - for plotting results
  - plot.order: internal - for plotting results
  - plot.group: internal - for plotting results
  - parameter: Model parameter being estimated
  - model: Which model results are reported
    - BB-MM 1 Component, Method of Moments Estimate
    - BB-MML: 1 Component, Maximum Marginal Likelihood Estimate
    - BB Mixture: 2 Component Mixture, Maximum Marginal Likelihood Estimate
  - posterior_mean: Posterior mean, if it exists (see documentation)
  - posterior_variance: Posterior variance, if it exists (see documentation)
  - pq.0.01: posterior quantile 0.01
  - pq.0.025: posterior quantile 0.025
  - pq.0.05: posterior quantile 0.05
  - pq.0.1: posterior quantile 0.1
  - pq.0.25: posterior quantile 0.25
  - pq.0.5: posterior quantile 0.5
  - pq.0.75: posterior quantile 0.75
  - pq.0.9: posterior quantile 0.9
  - pq.0.95: posterior quantile 0.95
  - pq.0.975: posterior quantile 0.975
  - pq.0.99: posterior quantile 0.99
  - prior_mean: Prior mean for the model
  - prior_variance: Prior variance for the model
