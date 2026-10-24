DISCLAIMER: This is the code that I used to perform analysis relevant for my bachelor thesis. It is not particularly clean and it is not properly packaged, but the conclusions may be of some use for the curious.

# Comparsion of Markov-Switching GARCH (MS-GARCH) and The Stochastic Volatility (SV) for the financial volatility forecasting.
The goal is to compare the Bayesian fit of two discrete volatility based on its out-of-sample performance. The data is daily price of S&P500, with hourly data used to compute realized volatility for some metrics. MS-GARCH and SV, in their simplest forms, are estimated and sampled with R packages:

- A [fork](https://github.com/achervinskii/MSGARCH) of [MSGARCH](https://github.com/keblu/MSGARCH/)
- [stochvol](https://gregorkastner.github.io/stochvol/)

Draws from predictive distribution are evaluated against the realization in (pseudo-) out-of-sample comparison using different versions of of (quasi-) log-likelihood.
## Usage
The MSGARCH package has to be installed from the [fork](https://github.com/achervinskii/MSGARCH) that allows to extract conditional variances draws from the predictive distribution. Other dependencies can be installed from CRAN.

The scripts that were used to create outputs are located in the `scripts` folder. The folder `R` contains the code that was initially intended to be re-used across scripts. More specifically, 
- Out-of-sample forecasting exercise for both models is implemented in [compare_fit_sv_msgarch.R](./script/compare_fit_sv_msgarch.R)
- A primitive data fetching is implemented in [fetching_data.R](./R/fetching_data.R). It relies on Yahoo Finance endpoint, and thus its usage is subject to [Yahoo Developer API Terms of Use](https://legal.yahoo.com/us/en/yahoo/terms/product-atos/apiforydn/index.html)

## Poster
The [Poster](./ThesisPoster.pdf) contains a brief description of methodology, results and conclusions.

## Conclusions
- MS-GARCH struggles to identify the regime switches. Regimes are unstable and not persistent in MCMC draws. Based on this alone, SV seems more appropriate.
- The out-of-sample comparison is not informative, likely because the assumption of conditional normality is not appropriate for the data, in both SV and MS-GARCH. Because of that, the values of metrics used to compare the models are driven by outliers, as the Gaussian likelihood evaluates to an extreme number.

## Limitations and further research
Among others, the largest limitation is using the conditional normal distribution. In fact, the packages support other distributions, and a back-of-the-envelope trial yield more encouraging results for MS-GARCH. Besides, the authors of MS-GARCH package use non-Gaussian distribution when demonstrating the capabilities in jstatsoft article, which potentially points to a systematic problem with using the Gaussian likelihood in MS-GARCH. The latter may be a topic to pursue further.

