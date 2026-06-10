DISCLAIMER: This project is in an early development phase. It uses a (yet unpublished) fork of the MSGARCH package.  

# Comparsion of Markov-Switching GARCH (MS-GARCH) and The Stochastic Volatility (SV) for the financial volatility forecasting

It is a part of my Bachelor thesis project. The goal is to compare the Bayesian fit of two discrete volatility models using a consistent score. Both models are estimated and sampled with R packages:

- [MSGARCH](https://github.com/keblu/MSGARCH/)
- [stochvol](https://gregorkastner.github.io/stochvol/)

## Usage

~~To setup install packages from renv.lock into the local environment's library, use~~

```r
install.packages("renv")
renv::restore()
```

Forecasting routine for both models is implemented in [compare_fit_sv_msgarch.R](./script/compare_fit_sv_msgarch.R)

A primitive data fetching is implemented in [fetching_data.R](./R/fetching_data.R). It relies on Yahoo Finance endpoint, and thus its usage is subject to [Yahoo Developer API Terms of Use](https://legal.yahoo.com/us/en/yahoo/terms/product-atos/apiforydn/index.html)

## At a glance

soon...


