# Joint Modeling for longitudinal and time-to-event data with cuRe fraction (JMcuR)

**JMcuR** is a R-package dealing joint models for longitudinal and time-to-event data when the population of interest includes a cure fraction. 
Estimation step uses JAGS and predictions using Bayesian or Frequentist approach. 
This package is built from the JMbayes package version 0.4-1 implemented by Dimitris Rizopoulos.

More general details about these joint longitudinal-cure models can be found in the following reference:

Barbieri, A., & Legrand, C. (2019). Joint longitudinal and time-to-event cure models for the assessment of being cured. *Statistical Methods in Medical Research*. https://doi.org/10.1177/0962280219853599

To try the current development version from github, use:

```{r} 
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")}
devtools::install_github("AntoineBbi/JMcuR")
 ```
**Warning:** JMcuR package requires JAGS software (http://mcmc-jags.sourceforge.net/). 
