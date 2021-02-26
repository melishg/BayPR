Assessment of Disease Liability of Genetic Variants:
================
Josh Betz
2021-02-25 01:25



  - [Assessing Disease Liability from Population-level
    Data](#assessing-disease-liability-from-population-level-data)
      - [Notation](#notation)
      - [Measuring Associations with Disease
        Status](#measuring-associations-with-disease-status)
      - [Relative Risk, Odds Ratio, and Prevalence
        Ratio](#relative-risk-odds-ratio-and-prevalence-ratio)
  - [Statistical Inference for the Prevalence
    Ratio](#statistical-inference-for-the-prevalence-ratio)
      - [Binomial Likelihood](#binomial-likelihood)
          - [Obtaining the Prevalence
            Ratio](#obtaining-the-prevalence-ratio)
      - [The Beta-Binomial Model](#the-beta-binomial-model)
          - [Bayes and Shrinkage](#bayes-and-shrinkage)
          - [The Beta Prime Distribution](#the-beta-prime-distribution)
          - [The Effect of \(r_{k}\) on the
            Prior](#the-effect-of-r_k-on-the-prior)
      - [Empirical Bayes Estimation](#empirical-bayes-estimation)
      - [Using a Mixture Model](#using-a-mixture-model)
  - [Simulated Example](#simulated-example)
      - [Single Beta-Binomial
        Component](#single-beta-binomial-component)
      - [Two Component Beta-Binomial
        Mixture](#two-component-beta-binomial-mixture)

One way of assessing the disease liability of genetic variants is
assessing the frequency of variants in case and reference populations.
Statistical inference allows us to use the prevalence of variants in
samples from various populations to compare the prevalence of variants
across these populations. The number of possible variants may be very
large, and the frequency of some variants can be very small. These rare
variants pose a particular challenge to traditional statistical
approaches: these methods may perform poorly or fail when few or no
variants are observed in either population. In order to address the
issue of rare variants, a Bayesian statistical approach can be used.

# Assessing Disease Liability from Population-level Data

In this setting, we have a sample of \(n_{1k}\) individuals from one
population (i.e. cases of a disease, such as Cystic Fibrosis) and and
\(n_{0k}\) individuals assessed for the variant have been sampled from
another population (i.e. some reference population, such as GnomAD),
where \(k = 1, 2, \ldots, K\) indexes the \(K\) variants of interest.
The number of individuals assessed for a given variant in each sample.

The sample sizes for each variant \(n_{0k}\) and \(n_{1k}\) are fixed by
the design of the study: the number of variants of type \(k\) observed
in cases (denoted \(Y_{1k}\)) and the number of variants of type \(k\)
observed in a sample from a reference population (denoted \(Y_{0k}\))
are random. We would like to use the frequency of variant \(k\) in each
population (\(Y_{1k}/n_{1k}\) in cases and \(Y_{0k}/n_{0k}\) in the
reference population) to compare the population prevalence of the
variant between cases and the reference population.

## Notation

  - \(k = 1, 2, \ldots, K\): variant being assessed for disease
    liability
  - \(\pi_{0k}\): Prevalence of variant \(k\) in reference population
  - \(Y_{0k}\): Number of individuals with for variant \(k\) in
    reference sample
  - \(n_{0k}\): Number of individuals assessed for variant \(k\) in
    reference sample
  - \(\pi_{1k}\): Prevalence of variant \(k\) in reference population
  - \(Y_{1k}\): Number of individuals with for variant \(k\) variant in
    reference sample
  - \(n_{1k}\): Number of individuals assessed for variant \(k\) in
    reference sample
  - \(r_{k} = n_{1}/n_{0}\): Relative sample size - reference sample
    size relative to case sample size
  - \(T_{k} = Y_{0k} + Y_{1k}\): Total count of variant \(k\) in both
    samples
  - \(\gamma_{k} = \pi_{1k}/\pi_{0k}\): Population prevalence ratio for
    variant \(k\) in cases relative to the reference population

## Measuring Associations with Disease Status

Measures of association between a binary outcome and binary exposure
include the relative risk of disease in those exposed relative to those
unexposed \((RR_{D \vert E})\), and the odds ratio \((OR)\). In this
setting, the outcome is disease prevalence, and the exposure is the
presence of the genetic variant of interest. Often we are most
interested in the relative risk, a ratio of the risk of disease \((D+)\)
among exposed \((E+)\) to those unexposed \((E-)\):

\[ RR_{D \vert E} = \frac{Pr\{\text{Disease} \vert \text{Exposed}\}}{Pr\{\text{Disease} \vert \text{Unexposed}\}} = \frac{Pr\{D+ \vert E+\}}{Pr\{D+ \vert E-\}} = \frac{Pr\{D+, E+\}/Pr\{E+\}}{Pr\{D+, E-\}/Pr\{E-\}} \]

The relative risk conditions on a particular exposure, whereas the
prevalence ratio conditions on the disease state. With genetic variants,
there may be many possible types of exposures, each with different
potential effects. Let’s denote the variant by \(V\) and we are
interested in a particular variant, denoted \(k\):

\[ RR = \frac{Pr\{D+ \vert E+\}}{Pr\{D+ \vert E-\}} = \frac{Pr\{D+ \vert V=k\}}{\sum_{v\neq k}Pr\{D+ \vert V=v\}} \]

Case-control studies are often more feasible than prospective studies,
especially in the setting of rare diseases. In these studies, it is
often easier to get probabilities or odds of exposure given the disease
state rather than on the exposure status. The odds ratio \((OR)\)
measures the relative odds of disease among those exposed relative to
those unexposed.

\[ OR_{D \vert E} = \frac{Pr\{D+ \vert E+\}/Pr\{D- \vert E+\}}{Pr\{D+ \vert E-\}/Pr\{D- \vert E-\}} = \frac{Pr\{D+, E+\}Pr\{D-,E-\}}{Pr\{D+, E-\}Pr\{D-, E+\}} \]

The odds ratio is symmetric: the odds ratio of disease given exposure is
the same as the odds ratio of exposure given disease.

\[ OR_{D \vert E} = \frac{Pr\{D+, E+\}Pr\{D-,E-\}}{Pr\{D+, E-\}Pr\{D-, E+\}} = \frac{Pr\{E+ \vert D+\}/Pr\{E- \vert D+\}}{Pr\{E+ \vert D-\}/Pr\{E- \vert D-\}} = OR_{E \vert D} \]

When the prevalence of the disease is low, the \(OR \approx RR\).

Another measure of association is the prevalence ratio of the exposure
in those with the disease relative to those without
\((PR_{E \vert D})\):

\[PR_{E \vert D} = \frac{Pr\{\text{Exposed} | \text{Disease}\}}{Pr\{\text{Exposed}|\text{No Disease}\}}\]

If the exposure is a genetic variant that’s considered a potential
diagnostic marker or test \(T\), the prevalence ratio of the exposure in
those with the disease to those without is the same as the *positive
diagnostic likelihood ratio*:

\[LR_{+} = \frac{Pr\{T+ \vert D+\}}{Pr\{T+ \vert D-\}} = \frac{\text{sensitivity}}{1\,-\,\text{specificity}}\]

Values of these quantities near 1 indicate negligible association
between an outcome and an exposure, values below 1 indicate a negative
association between the outcome and the exposure (i.e. the outcome is
less common among the exposed), and values above 1 indicate a positive
association between the outcome and the exposure (i.e. the outcome is
more common among the exposed).

## Relative Risk, Odds Ratio, and Prevalence Ratio

Given an appropriate samples from populations of interest, we can
calculate these quantities as follows from a (2x2) table:

| Exposure  |     Disease      |    No Disease    |    Total     |
| :-------: | :--------------: | :--------------: | :----------: |
|  Exposed  | \(Pr\{E+, D+\}\) | \(Pr\{E+, D-\}\) | \(Pr\{E+\}\) |
| Unexposed | \(Pr\{E-, D+\}\) | \(Pr\{E+, D-\}\) | \(Pr\{E-\}\) |
|   Total   |   \(Pr\{D+\}\)   |   \(Pr\{D-\}\)   |              |

These are commonly denoted with letters to simplify calculations:

| Exposure  |  Disease  | No Disease |   Total   |
| :-------: | :-------: | :--------: | :-------: |
|  Exposed  |   \(A\)   |   \(B\)    | \(A + B\) |
| Unexposed |   \(C\)   |   \(D\)    | \(C + D\) |
|   Total   | \(A + C\) | \(B + D\)  |           |

Using the notation above, this becomes:

|    Variant     |          Cases          |        Controls         |               Total               |
| :------------: | :---------------------: | :---------------------: | :-------------------------------: |
| Variant \(k\)  |     \(Y_{0k} = A\)      |     \(Y_{1k} = B\)      |         \(T_{k} = A + B\)         |
| Other Variants | \(n_{0k} - Y_{0k} = C\) | \(n_{1k} - Y_{1k} = D\) | \(n_{0} + n_{1} - T_{k} = C + D\) |
|     Total      |   \(n_{0k} = A + C\)    |   \(n_{1k} = B + D\)    |                                   |

\[ RR_{D \vert E} = \frac{Pr\{D+ \vert E+\}}{Pr\{D+ \vert E-\}} =  \frac{A/(A + B)}{C/(C + D)} = \frac{A(C+D)}{C(A+B)} = \frac{AC + AD}{AC + BC} = \frac{Y_{0k}/T_{k}}{(n_{0k} - Y_{0k})/(n_{0} + n_{1} - T_{k})} \]

\[ OR_{D \vert E} = \frac{A/B}{C/D} = \frac{AD}{BC} = \frac{Y_{0k}/Y_{1k}}{(n_{0k} - Y_{0k})/(n_{1k} - Y_{1k})} \]

\[ PR_{E \vert D} = \frac{Pr\{E+ \vert D+\}}{Pr\{E+ \vert D-\}}  = \frac{A/(A+C)}{B/(B+D)} = \frac{A(B+D)}{B(A+C)} = \frac{Y_{0k}/n_{0k}}{Y_{1k}/n_{1k}} \]

# Statistical Inference for the Prevalence Ratio

## Binomial Likelihood

When we have a fixed sample size of \(n\) independent individuals, each
from a population with prevalence \(\pi\), we can model the number of
prevalent cases in this sample using the binomial likelihood model:

\[Y \sim Bin(n, \pi):\,Pr\{Y = y|n, \pi\} = {n \choose y} \pi^{y} (1 - \pi)^{n-y}\]

This binomial likelihood can be approximated using a Poisson
distribution with rate parameter \(\lambda = n\pi\) if \(n\), the number
of observations, is large and \(\pi_{k}\), the prevalence of the
variant, is small.

\[Y \sim Poisson(\lambda = n\pi):\,Pr\{Y = y|n, \pi\} = \frac{(n\pi)^{y}e^{-n\pi}}{y!}\]

If we observe \(Y_{1k}\) variants of type \(k\) in cases, \(Y_{0k}\)
variants of type \(k\) in controls, let \(T_{k}\) denote their sum. We
can model the proportion of type \(k\) variants that arose from cases
out of the the total number of type \(k\) variants observed, denoted
\(\theta_{k}\), using a binomial distribution:

\[Y_{0k} \sim Poisson(n_{0k}\pi_{0k}),\, Y_{1k} \sim Poisson(n_{1k}\pi_{1k});\quad Y_{0k} \perp Y_{1k}; \quad \left(Y_{1k} \vert Y_{0k} + Y_{1k} = T_{k}\right) \sim Bin(T_{k}, \theta_{k}),\]

where
\(\theta_{k} = (n_{1k}\pi_{1k})/(n_{1k}\pi_{1k} + n_{0k}\pi_{0k})\), the
rate of occurrences in cases divided by the total rate of occurrences
(the sum of the rates in cases and controls).

Note that since \(\theta_{k}\) is a proportion:
\(0 \le \theta_{k} \le1\). We are interested in obtaining the
*prevalence ratio* \(\gamma_{k} = \pi_{1k}/\pi_{0k}\), where
\(0 < \gamma_{k}\). We can infer about this quantity by using a
transformation.

### Obtaining the Prevalence Ratio

We can divide the numerator and denominator by \(\pi_{0k}\) to
parameterize this distribution by the prevalence ratio \(\gamma_{k}\):

\[\theta_{k} = \frac{n_{1k}\pi_{1k}/\pi_{0k}}{n_{1k}\pi_{1k}/\pi_{1k} + n_{0k}\pi_{0k}/\pi_{1k}} = \frac{n_{1k}\gamma_{k}}{n_{1k}\gamma_{k} + n_{0k}}\]

Let \(r_{k} = n_{1k}/n_{0k}\) denote the ratio of the sample size of
cases relative to controls. Dividing the numerator and denominator again
by \(n_{0k}\), the control sample size, gives:

\[\theta_{k} = \frac{n_{1k}\gamma_{k}/n_{0k}}{n_{1k}\gamma_{k}/n_{0k} + n_{0k}/n_{0k}} = \frac{r_{k}\gamma_{k}}{r_{k}\gamma_{k} + 1}\]

Parameterizing the model this way allows the modeling multiple variants,
allowing for variation in sample sizes across different variants, as
long as their ratio of sample sizes for any given variant \(r_{k}\) is
comparable.

Since \(\theta_{k} = r_{k}\gamma_{k}/(r_{k}\gamma_{k} + 1)\), we can
rearrange this to get back to the prevalence ratio:
\(\gamma_{k} = \frac{1}{r_{k}} \frac{\theta_{k}}{(1 - \theta_{k})}\).

So we infer about \(\theta_{k}\), and then use the above transformation
to obtain the prevalence ratio \(\gamma_{k}.\) Note that this
transformation depends on the relative samples sizes
\(r_{k} = n_{1k}/n_{0k}\).

## The Beta-Binomial Model

Bayesian statistics augments the likelihood probability model with a
prior probability model. Our model for the data is a binomial
distribution, which depends on \(n = T_{k}\), the number of variants
observed and \(\pi = \theta_{k}\), the proportion of variants arising
from the population of cases. Instead of viewing each of our \(K\)
variants separately, we could view them as a sample from a population of
genetic variants, and the proportions of variants due to cases in this
population follow a probability distribution called the Beta
distribution.

The shape of the beta distribution is controlled by two parameters,
\(\alpha\) and \(\beta\). Depending on the values these parameters, the
beta distribution can take on many possible shapes, as seen in the
figures below. The mean of the beta distribution is given by
\(E[\theta_{k} \vert \alpha, \beta] = \mu = \alpha/(\alpha + \beta)\),
and its variance is given by
\(var[\theta_{k} \vert \alpha, \beta] = \alpha \beta/\left((\alpha + \beta)^2(\alpha + \beta + 1)\right)\).

![**Figure 1:** Examples of beta distributions obtained by varying the
values of the parameters alpha, indicated by the color, and beta,
indicated by the figure
panel.](Prevalence_Ratio_210225_files/figure-gfm/beta_example-1.png)

**Figure 1:** Examples of beta distributions obtained by varying the
values of the parameters alpha, indicated by the color, and beta,
indicated by the figure panel.

We can also re-parameterize the beta distribution in terms of \(\mu\),
its mean, and \(M = \alpha + \beta\):

\[Pr\{\theta_{k} \vert \mu, M\} = \frac{\Gamma(M)}{\Gamma(\mu M)\Gamma((1-\mu)M)} \theta^{\mu M - 1}(1 - \theta)^{(1-\mu)M - 1}\]

In this parameterization, \(\alpha = \mu M\) and \(\beta = (1-\mu)M\).
We can think of \(M\) as a prior sample size: \(\alpha\) prior variants
among cases and \(\beta\) prior variants among controls. When
parameterized this way, the variance of the beta distribution is given
by \(var[\theta_{k} \vert \mu,\,M] = \mu(1 - \mu)/(M + 1)\).

The combination of a beta distribution for the proportion of variants
arising from cases, and the binomial likelihood for the number of
variants observed in cases, gives a Beta-Binomial Distribution. The
parameters of the beta distribution are called hyperparameters, because
their values determine the distribution of the proportion parameter in
the binomial model.

The probability density function for the prior is:

\[Pr\{\theta_{k} | \alpha,\,\beta\} = \frac{\Gamma(\alpha + \beta)}{\Gamma(\alpha)\Gamma(\beta)} \theta^{\alpha - 1}_{k}(1 - \theta_{k})^{\beta - 1}\]

The likelihood for the observed data is:

\[Pr\{Y_{1k} |T_{1k},\,\theta_{k}\} = {T_{1k} \choose Y_{1k}} \theta_{k}^{Y_{1k}} (1 - \theta_{k})^{Y_{1k} - T_{1k}}\]

Note the similarity between these two functions: both are products of
the terms \(\theta_{k}\) and \((1 - \theta_{k})\). When combined using
Bayes’ Rule, they give a posterior distribution is also a beta
distribution:

\[Pr\{\theta_{k}|Y_{1k},\,T_{k},\,\alpha,\,\beta\} = \frac{\Gamma(\alpha + \beta + T_{k})}{\Gamma(\alpha + Y_{1k})\Gamma(\beta+ Y_{0k})} \theta^{\alpha + Y_{1k} - 1}(1 - \theta)^{\beta + Y_{0k} - 1}\]

This is just a beta distribution with \(\alpha'= \alpha + Y_{1k}\) and
\(\beta'= \beta + Y_{0k}\): this is why we can interpret \(\alpha\) as
the prior number of variants observed in cases, and \(\beta\) as the
prior number of variants observed among controls. The posterior
distribution, and its summaries (such as its mean, variance, and
quantiles), allow us to infer about \(\theta_{k}\), the proportion of
variants due to cases, for each of our \(K\) variants.

### Bayes and Shrinkage

The advantage of the Bayesian approach is that the prior distribution
acts as a stabilizing influence when the number of observed variants is
small, and this influence diminishes as the number of observed variants
becomes larger. The effect of the prior can be more easily understood
when viewed through the parameterization of \(\mu\) and \(M\). The mean
of the posterior distribution is given by
\(E[\theta_{k} \vert \alpha, \beta] = \mu = \alpha/(\alpha + \beta)\).
Re-writing this in terms of \(\mu\) and \(M\) gives:

\[E[\theta_{k} \vert \alpha, \beta] = \frac{\alpha + Y_{1k}}{\alpha + \beta + T_{k}} = \frac{\mu M + Y_{1k}}{M + T_{k}} = \frac{M}{T_{k} + M}\mu + \frac{T_{k}}{T_{k} + M}\left(\frac{Y_{1k}}{T_{k}}\right) = s_{k}\mu + (1 - s_{k})\left(\frac{Y_{1k}}{T_{k}}\right)\]

Since \(M\) acts as a prior sample size, and \(T_{k}\) is the total
number of variants observed, the quantity \(s_{k} = M/(M+T_{K})\)
represents the proportion of the information coming from the prior, and
\((1 - s_{k}) = T_{k}/(M+T_{K})\) represents the proportion of
information coming from the observed data. Here we see the posterior
mean is a combination of the prior mean \(\mu\) and the proportion of
variants among cases to the total number of variants observed
\((Y_{1k}/T_{k})\), each weighted according their sample size
contribution. When \(T_{k}\), the number of observed variants is small,
the posterior mean is ‘pulled’ or ‘shrunken’ towards the prior mean:
this gives more stable estimates. As the number of observed variants
becomes increasingly larger than the prior sample size, the effect of
the prior diminishes.

![**Figure 2:** An example of how the shape of the prior affects the
posterior mean, and how this depends on both the number of observed
variants, the proportion of variants seen in cases, and the prior sample
size. The prior distribution is shown in gray. When the prior sample
size (M) is small, the effect of the prior is smallest, and the
posterior means are very close to the direct estimates. As the prior
sample size increases, direct estimates are ‘pulled’ or ‘shrunken’
towards the prior mean (here, 0.5), with greater shrinkage occurring
when fewer variants are
observed.](Prevalence_Ratio_210225_files/figure-gfm/shrinkage_example-1.png)

**Figure 2:** An example of how the shape of the prior affects the
posterior mean, and how this depends on both the number of observed
variants, the proportion of variants seen in cases, and the prior sample
size. The prior distribution is shown in gray. When the prior sample
size (M) is small, the effect of the prior is smallest, and the
posterior means are very close to the direct estimates. As the prior
sample size increases, direct estimates are ‘pulled’ or ‘shrunken’
towards the prior mean (here, 0.5), with greater shrinkage occurring
when fewer variants are observed.

### The Beta Prime Distribution

The beta distribution describes the distribution of a random variable
over the interval \((0, 1)\), which is convenient for describing a
proportion or probability. However, if we want to infer about the
prevalence ratio \(\gamma_{k}\), we need to transform back to the
prevalence ratio scale:

\[\frac{\theta_{k}}{r_{k}(1-\theta_{k})} = \gamma_{k}\]

Notice that when \(r_{k} = 1\) (sample sizes are identical between cases
and controls), this transformation is from the probability scale to the
odds scale. If we have a beta random variable, and transform it to the
odds scale, this results in a variable with the beta prime distribution.
The shape of this distribution is governed by the same parameters,
\(\alpha\) and \(\beta\): If \(X \sim Beta(\alpha,\,\beta)\) then
\(Y = X/(1-X) \sim Beta\,Prime(\alpha,\,\beta)\).

If \((Y \vert \alpha, \beta) \sim Beta\,Prime(\alpha, \beta)\), its mean
is given by \(E[Y \vert \alpha, \beta] = \alpha/(\beta - 1)\), and its
variance is given by
\(var[Y \vert \alpha, \beta] = \alpha(\alpha + \beta - 1)/\left((\beta - 2)(\beta - 1)^2\right)\).
Note that if \(\beta \le 1\), the mean of this distribution is not
finite, and if \(\beta \le 2\) its variance is not finite.

For a beta-prime posterior distribution of prevalence ratios, we can
plug in our posterior values of \(\alpha'= \alpha + Y_{1k}\) and
\(\beta'= \beta + Y_{0k}\). When \(r_{k}\neq 1\) (the sample size
differs between cases and controls), the posterior becomes a scaled
version of the beta prime distribution. The mean is scaled by
\(1/r_{k}\), and the variance is scaled by \(1/r_{k}^2\):

\[E[\gamma_{k} \vert Y_{1k}, T_{k}, \alpha, \beta] = \frac{1}{r_{k}} \frac{\alpha + Y_{1k}}{(\beta + Y_{0k} - 1)} = \frac{\frac{\alpha}{n_{1k}} + \frac{Y_{1k}}{n_{1k}}}{\frac{\beta - 1}{n_{0k}} + \frac{Y_{0k}}{n_{0k}}}\]

By adding \(\alpha/n_{1k}\) to the numerator, and \((\beta-1)/n_{0k}\)
to the denominator, this stabilizes the prevalence ratio when the number
of observed variants is smaller, and this stabilizing factor decreases
as the sample size becomes larger. This can also be seen when we
paramaterize instead using \(\mu\) and \(M\):

\[E[\gamma_{k} \vert Y_{1k}, T_{k}, \mu, M] = \frac{\mu\frac{M}{n_{1k}} + \frac{Y_{1k}}{n_{1k}}}{(1 - \mu)\frac{M}{n_{0k}} - \frac{1}{n_{0k}} + \frac{Y_{0k}}{n_{0k}}}\]

The numerator is pushed towards the prior mean \(\mu\), and this effect
diminishes as \(n_{1k}\), the sample size among cases, becomes large
compared to the prior sample size \(M\). The denominator is pushed
towards \((1 - \mu)\), and this effect diminishes as \(n_{0k}\), the
sample size among controls, becomes large compared to the prior sample
size \(M\).

![**Figure 3:** Examples of two beta distributions, modeling the
proportion of variants observed in cases, and their corresponding
distributions of prevalence
ratios.](Prevalence_Ratio_210225_files/figure-gfm/Theta_vs_Gamma_Scale-1.png)

**Figure 3:** Examples of two beta distributions, modeling the
proportion of variants observed in cases, and their corresponding
distributions of prevalence ratios.

### The Effect of \(r_{k}\) on the Prior

Note that the transformation to the prevalence ratio scale depends on
\(r_{k}\), the ratio of sample sizes in cases \((n_{1})\) relative to
controls \((n_{0})\). Note how the same beta model could represent
populations of variants with lower, equal, or higher prevalence in cases
relative to controls, depending on the value of \(r_{k}:\)

![**Figure 4:** Example of how one beta distribution, could represent
lower, equal, or higher prevalence between cases and controls, depending
on the relative sample size between cases and
controls.](Prevalence_Ratio_210225_files/figure-gfm/Variation_in_r_k-1.png)

**Figure 4:** Example of how one beta distribution, could represent
lower, equal, or higher prevalence between cases and controls, depending
on the relative sample size between cases and controls.

This also illustrates the importance of the assumption about the
variation \(r_{k}\) across all \(K\) variants: if there is appreciable
variation in this ratio across variants, then variants will essentially
be ‘pulled’ or ‘shrunken’ in different directions.

## Empirical Bayes Estimation

Up until now, we have been treating the parameters of our beta prior
distribution, \(\alpha\) and \(\beta\), as known quantities, and showing
how their values ‘pull’ or ‘shrink’ the direct estimates towards the
prior mean on on the proportion (or \(\theta_{k}\)) scale, and ‘pull’ or
‘shrink’ the numerator and denominator on the ratio scale. But how do we
determine suitable values for these ‘hyperparameters?’

Rather than either supplying exact values for these parameters, or
specifying a prior distribution over these parameters, we can find the
values of these parameters that maximize the *marginal likelihood* of
the observed data, averaging over the parameter \(\theta_{k}\):
\(Pr\{Y_{1k} \vert T_{k}, \, \mu, \,M\}\)

\[Pr\{Y_{1k} \vert T_{k}, \, \mu, \,M\} = f(Y_{1k} \vert T_{k}, \, \mu, \,M) =
  \frac{\Gamma(T_{k} + 1)}
  {\Gamma\left(Y_{1k} + 1 \right) \Gamma\left(T_{k} - Y_{1k} + 1\right)}
  \frac{\Gamma(Y_{1k} + \mu M) \Gamma\left(T_{k} - Y_{1k} + (1 - \mu)M \right)}
  {\Gamma\left(T_{k} + M \right)}
  \frac{\Gamma\left(M \right)}
  {\Gamma(\mu M) \Gamma\left((1-\mu)M\right)}\]

where \(\Gamma\) denotes the gamma function.

## Using a Mixture Model

Instead of assuming that all variants are represented by one beta
distribution that describes the proportion of variants among cases, we
can relax this assumption by assuming that the observed data were
generated from a mixture of different beta distributions.

Instead of our model having only two parameters, \(\mu\) and \(M\), our
mixture model will have 5 parameters: \(\mu_{1}\) and \(M_{1}\), which
control the shape of one beta distribution, \(\mu_{2}\) and \(M_{2}\),
which control the shape of the second beta distribution, and
\(\epsilon\), which is the proportion of variants arising from the
second beta distribution:

\[f(Y_{1k} \vert T_{k}, \, \mu_{1}, \,M_{1}, \, \mu_{2}, \,M_{2}, \, \epsilon) = \epsilon f(Y_{1k} \vert T_{k}, \, \mu_{1}, \,M_{1}) + (1 - \epsilon)f(Y_{1k} \vert T_{k}, \, \mu_{2}, \,M_{2})\]

Note that two different parameter vectors
\((\mu_{1} = a, \,M_{1} = b, \, \mu_{2} = c, \,M_{2} = d, \epsilon = e)\)
and
\((\mu_{1} = c, \,M_{1} = d, \, \mu_{2} = a, \,M_{2} = b, \epsilon = 1-e)\)
result in identical values of the mixture model: for this reason, the
constraint \(\epsilon < 0.5\) is imposed to identify a unique solution.
Estimates of these parameters are obtained by maximizing the marginal
likelihood of the mixture.

Instead of each direct estimate being ‘pulled’ or ‘shrunken’ in the same
direction, each distribution will exert a different ‘pull’ on the data,
with the ‘pull’ being related to the relative compatibility between each
model component and the data.

From this mixture distribution, in addition to obtaining posterior
means, variances, and quantiles, we can additionally obtain:

  - The marginal likelihood ratio:
    \(MLR_{2/1} = Pr\{Y_{1k} \vert T_{k}, \, \mu_{2}, \,M_{2}\}/Pr\{Y_{1k} \vert T_{k}, \, \mu_{1}, \,M_{1}\}\)
  - The posterior odds of belonging to component 2 vs. 1:
    \(PO_{2k}=(\epsilon/(1 - \epsilon))MLR_{2/1}\)
  - The posterior probability of belonging to component 2:
    \(P_{2k} = PO_{2k}/(1 + PO_{2k})\)

![**Figure 5:** An example of a mixture prior, and the effect of each
mixture component on the direct estimate. When the prior sample size (M)
is small, the effect of the prior is smallest, and the posterior means
are very close to the direct estimates. As the prior sample size
increases, direct estimates are ‘pulled’ or ‘shrunken’ towards each of
the priors. The ‘pull’ of each prior depends on the marginal likelihood
ratio: the relative compatibility between the direct estimate and each
prior probability
component.](Prevalence_Ratio_210225_files/figure-gfm/shrinkage_mix_example-1.png)

**Figure 5:** An example of a mixture prior, and the effect of each
mixture component on the direct estimate. When the prior sample size (M)
is small, the effect of the prior is smallest, and the posterior means
are very close to the direct estimates. As the prior sample size
increases, direct estimates are ‘pulled’ or ‘shrunken’ towards each of
the priors. The ‘pull’ of each prior depends on the marginal likelihood
ratio: the relative compatibility between the direct estimate and each
prior probability component.

# Simulated Example

We can construct simulated data to illustrate the use of the beta
binomial model in practice. In this first simulation, we will assume
that variants differ in their prevalence ratios, which are sampled from
a single beta distribution. Here we are assuming that we have \(n_0\) =
2.510^{5} individuals sampled from a reference population, and \(n_1\) =
500 individuals sampled from a population of cases: this gives \(r_{k}\)
= 0.002.

## Single Beta-Binomial Component

Here we are creating a population of variants whose mean prevalence
ratio is 1 and the variance is 0.025: this translates to a beta
distribution with parameters \(\alpha\) = 40.082 and \(\beta\) =
2.004210^{4}. We are observing a sample of 200 variants from this
population. The soft

![](Prevalence_Ratio_210225_files/figure-gfm/Single_BB-1.png)<!-- -->

**Figure 6:** A hypothetical population of variants whose mean
prevalence ratio is 1 and whose variance of prevalence ratios is 0.025:
this is shown on the prevalence ratio scale (right) and on the
transformed scale (left). This corresponds to a beta distribution with
parameters \(\alpha\) = 40.082 and \(\beta\) = 2.004210^{4}.

With this sample of variants, we can fit the empirical Bayes beta
binomial model to this data:

``` r
bb.mml <-
  bb.prevalence.ratio.mml(
    data = sim.data,
    y.0 = "y.0", n.0 = "n.0",
    y.1 = "y.1", n.1 = "n.1",
    id = "event"
  )
```

The software gives fitted values of \(\hat{\alpha}\) = 65.6411536
(actual value 40.082) and \(\hat{\beta}\) = 3.368775910^{4} (actual
value 2.004210^{4}), which corresponds to a prior mean for the
prevalence ratio of 0.974 (actual value 1) and a prior variance for the
prevalence ratio of 0.014 (actual value 0.025). The result also contains
a `data.frame` summarizing the posterior distribution of the prevalence
ratio for each variant. We can use the 95% equal tail credible intervals
to infer about the prevalence ratio for each variant.

| event |    n.0 |   y.0 | n.1 | y.1 | True Prevalence Ratio | Sample Estimate | EB 95% ETCI Lower | EB Posterior Median | EB 95% ETCI Upper |
| ----: | -----: | ----: | --: | --: | --------------------: | --------------: | ----------------: | ------------------: | ----------------: |
|     1 | 250000 | 10473 | 500 |  19 |                  0.83 |            0.91 |              0.76 |                0.95 |              1.17 |
|     2 | 250000 | 13413 | 500 |  31 |                  1.00 |            1.16 |              0.83 |                1.02 |              1.24 |
|     3 | 250000 | 12525 | 500 |  18 |                  0.77 |            0.72 |              0.72 |                0.90 |              1.11 |
|     4 | 250000 | 12193 | 500 |  20 |                  1.05 |            0.82 |              0.75 |                0.93 |              1.14 |
|     5 | 250000 | 13881 | 500 |  17 |                  0.83 |            0.61 |              0.69 |                0.87 |              1.07 |
|     6 | 250000 | 12091 | 500 |  25 |                  1.15 |            1.03 |              0.80 |                0.99 |              1.20 |
|     7 | 250000 | 15848 | 500 |  33 |                  1.14 |            1.04 |              0.81 |                0.99 |              1.20 |
|     8 | 250000 | 13001 | 500 |  25 |                  0.91 |            0.96 |              0.78 |                0.97 |              1.18 |
|     9 | 250000 | 12050 | 500 |  24 |                  1.23 |            1.00 |              0.79 |                0.98 |              1.19 |
|    10 | 250000 | 12599 | 500 |  28 |                  0.95 |            1.11 |              0.82 |                1.01 |              1.23 |

Empirical Bayes estimates using the single component Beta Binomial
model: the Sample Estimate uses only the sample prevalence data for each
variant. The Empirical Bayes model pools information across variants to
obtain estimates that perform well when observed counts are sparse, even
when evaluated by frequentist criteria (interval coverage, mean squared
error).

![](Prevalence_Ratio_210225_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

**Figure 6:** A sample of simulated variants from a single beta
distribution: the true value of the prevalence ratio for each variant is
shown in open circles. The sample estimate using only the data from each
variant individually, is shown in filled purple circles. The Empirical
Bayes estimate (the posterior median) is shown in green filled squares.
The Empirical Bayes 95% equal-tail credible interval is shown using a
line, with color and shading to indicate whether or not the estimate
covered the true prevalence ratio for that variant.

## Two Component Beta-Binomial Mixture

Rather than assuming all variants arise from one population of variants,
whose prevalence ratios are described by a single beta binomial
distribution, a mixture model allows for two populations of variants,
with 25% coming from a population of variants whose average prevalence
ratio is 3.

![](Prevalence_Ratio_210225_files/figure-gfm/Mixture_Model-1.png)<!-- -->

**Figure 7:** A hypothetical population of variants who arise from two
separate subpopulations, each described by a beta distribution. In one
subpopulation, the mean prevalence ratio is 1 and whose variance of
prevalence ratios is 0.025. In the other subpopulation, the mean
prevalence ratio is 3 and whose variance of prevalence ratios is 1.5.
This is shown on the prevalence ratio scale (right) and on the
transformed scale (left).

With our sample of variants, we can fit the empirical Bayes beta
binomial model to this data. This searches over a grid of possible
starting values to try and avoid convergence to local maxima. The grid
of values to search can be customized.

``` r
bb.mml.mix <-
  bb.2c.mixture.grid.search(
    data = sim.data.2,
    y.0 = "y.0", n.0 = "n.0",
    y.1 = "y.1", n.1 = "n.1",
    id = "event"
  )
```

For the first subpopulation, we get \(\hat{\alpha}_{1}\) = 35.264
(actual value 40.082) and \(\hat{\beta}_{1}\) = 1.761255810^{4} (actual
value 2.004210^{4}), which corresponds to a prior mean for the
prevalence ratio of 1.001 (actual value 1) and a prior variance for the
prevalence ratio of 0.028 (actual value 0.025).

For the second subpopulation, we get \(\hat{\alpha}_{2}\) = 15.217
(actual value 6.042) and \(\hat{\beta}_{2}\) = 2699.507 (actual value
1008), which corresponds to a prior mean for the prevalence ratio of
2.819 (actual value 3) and a prior variance for the prevalence ratio of
0.526 (actual value 1.5). The mixture model estimates that 23% of
variants come from this second population (actual value 25%).

The result also contains a `data.frame` summarizing the posterior
distribution of the prevalence ratio for each variant. We can use the
95% equal tail credible intervals to infer about the prevalence ratio
for each variant, as well as the marginal likelihood ratios or class
probabilities.

| event |    n.0 |   y.0 | n.1 | y.1 | True Prevalence Ratio | Sample Estimate | EB 95% ETCI Lower | EB Posterior Median | EB 95% ETCI Upper | Marginal LR - 2:1 | Probability - 2 |
| ----: | -----: | ----: | --: | --: | --------------------: | --------------: | ----------------: | ------------------: | ----------------: | ----------------: | --------------: |
|     1 | 250000 | 13467 | 500 |  20 |                  0.87 |            0.74 |              0.67 |                0.88 |              1.14 |              0.00 |            0.00 |
|     2 | 250000 | 12990 | 500 |  34 |                  1.10 |            1.31 |              0.88 |                1.13 |              1.45 |              0.09 |            0.02 |
|     3 | 250000 | 12274 | 500 |  16 |                  0.80 |            0.65 |              0.64 |                0.85 |              1.11 |              0.00 |            0.00 |
|     4 | 250000 | 12110 | 500 |  24 |                  0.86 |            0.99 |              0.76 |                0.99 |              1.27 |              0.01 |            0.00 |
|     5 | 250000 | 11566 | 500 |  51 |                  2.36 |            2.20 |              1.73 |                2.30 |              2.91 |            402.27 |            0.99 |
|     6 | 250000 | 12393 | 500 |  22 |                  1.03 |            0.89 |              0.72 |                0.95 |              1.22 |              0.00 |            0.00 |
|     7 | 250000 | 13705 | 500 |  37 |                  1.04 |            1.35 |              0.90 |                1.15 |              1.49 |              0.13 |            0.02 |
|     8 | 250000 | 14467 | 500 |  21 |                  0.70 |            0.73 |              0.66 |                0.87 |              1.12 |              0.00 |            0.00 |
|     9 | 250000 |  9581 | 500 |  21 |                  1.02 |            1.10 |              0.78 |                1.03 |              1.33 |              0.02 |            0.00 |
|    10 | 250000 | 11520 | 500 |  19 |                  0.93 |            0.82 |              0.70 |                0.93 |              1.20 |              0.00 |            0.00 |

Empirical Bayes estimates using the two component Beta Binomial mixture
model: the Sample Estimate uses only the sample prevalence data for each
variant. The Empirical Bayes model pools information across variants to
obtain estimates that perform well when observed counts are sparse, even
when evaluated by frequentist criteria (interval coverage, mean squared
error).
