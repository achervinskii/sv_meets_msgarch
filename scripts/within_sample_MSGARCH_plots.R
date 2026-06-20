# The goal is to create plots for the in-sample bayesian fit
# The states are of interest, as they are the basis of the forecasting ability in MSGARCH

library(stochvol)
library(ggplot2)
library(MSGARCH)
library(dplyr)
library(coda)
library(patchwork)
library(tidyr)

set.seed(1)

# data fetching. not fetched if exists per ticker
source("./R/fetching_data.R")
start_date <- "2016-01-01"
end_date <- "2026-06-13"

ticker <- "^GSPC"
df <- fetch_daily_prices(ticker, start_date, end_date)

# contsurct the return series
y <- df[["close"]]
selected <- !is.na(y)
y <- y[selected]
lret <- diff(log(y))
lret_dates <- df$datetime[selected]
lret_dates <- lret_dates[-length(lret_dates)]
ndraws <- 10000

spec <- CreateSpec(
  variance.spec = list(model = c("sGARCH", "sGARCH")),
  distribution.spec = list(distribution = c("norm", "norm"))
)

fit_ml <- FitML(spec = spec, data = lret)
ctr <- list(
  nburn = ndraws,
  nmcmc = ndraws,
  nthin = 1,
  par0 = fit_ml$par
)
fit_mcmc <- FitMCMC(spec = spec, data = lret, ctr = ctr)
mcmc_par <- fit_mcmc$par
avg_par <- colMeans(mcmc_par)
state <- State(fit_mcmc)
# exclude the last smoothed, consider the second regime
statesm <- state$SmoothProb[1:dim(state$SmoothProb)[1] - 1, , 2]
# collapse the simulations by averaging
if (!is.null(dim(statesm))) {
  probs <- rowMeans(statesm)
} else {
  probs <- statesm
}

# combine in df
df_plot <- data.frame(
  r = lret / (12 * sd(lret)) + 0.5,
  t = lret_dates,
  probs = probs
)

cols <- list(
  msgarch = adjustcolor("red"),
  sv = adjustcolor("#007712"),
  hl = adjustcolor("#ff6d00")
)

plot <- ggplot(
  df_plot[, ],
  aes(x = t, y = probs)
) +
  geom_point(aes(x = t, y = r, color = "return"),
    size = 0.6, alpha = 0.7
  ) +
  geom_line(aes(color = "P(s_t = 2)"), alpha = 0.7) +
  scale_color_manual(values = c(
    "return" = "black",
    "P(s_t = 2)" = cols$hl
  )) +
  labs(
    y = "Probability of the second state",
    x = "Date"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
plot

ggsave("output/msgarch_states_1.png", plot, width = 10, height = 5, dpi = 300)


# run 4 chains to obtain diagnostics
n <- 4
chains_ms <- vector("list", length = 4)
for (i in 1:4) {
  fit <- FitMCMC(spec = spec, data = lret, ctr = ctr)
  pars <- fit$par
  chains_ms[[i]] <- pars
  print(i)
}

coda_to_df <- function(chains_list) {
  # convert the coda objects to the dataframe and plot
  mcmclist <- as.mcmc.list(chains_list)
  effectiveSize(mcmclist)
  mcmcmatrix <- as.matrix(mcmclist, iters = TRUE, chains = TRUE)
  mcmcdf <- data.frame(mcmcmatrix) |>
    pivot_longer(cols = (!ITER & !CHAIN), values_to = "y", names_to = "par")
  mcmcdf
}

mcmc_ms_df <- coda_to_df(chains_ms)
mcmc_ms_df$ITER <- mcmc_ms_df$ITER + 10000
mcmc_ms_df <- mcmc_ms_df |>
  mutate(par = case_when(
    par == "alpha0_1" ~ "omega_1",
    par == "alpha0_2" ~ "omega_2",
    par == "alpha1_1" ~ "alpha_1",
    par == "alpha1_2" ~ "alpha_2",
    TRUE ~ par))

pars <- unique(mcmc_ms_df$par)
aligner <- c(rbind(c(1:3), c(4:6)), 7, 8)
pars <- pars[aligner]

plot_param_mcmc <- function(df, par_name) {
  out <- ggplot(
    filter(df, par == par_name),
    aes(x = ITER, y = y, color = as.factor(CHAIN))
  ) +
    # scale_colour_discrete(palette = scales::pal_brewer(palette = "Dark2")) +
    scale_colour_viridis_d()+
    xlab(NULL) +
    ylab(par_name) +
    geom_line(alpha = 0.8, linewidth = 0.4) +
    theme_bw() +
    labs(color = "") +
    theme(legend.position = "none")
}
plots <- lapply(pars, FUN = \(x) plot_param_mcmc(mcmc_ms_df, x))
plots_combined <- wrap_plots(plots, ncol = 2)
ggsave("output/mcmc_msgarch.png", plots_combined, width = 7, height = 7, dpi = 500)

chains_sv <- vector("list", length = 4)
for (i in 1:4) {
  fit <- svsample(
    lret,
    draws = ndraws,
    burnin = ndraws,
    n_chain = 1,
    quiet = TRUE
  )
  pars <- fit$para
  chains_sv[[i]] <- pars[[1]][, 1:3]
  print(i)
}

mcmc_sv_df <- coda_to_df(chains_sv)
pars <- unique(mcmc_sv_df$par)
plots <- lapply(pars, FUN = \(x) plot_param_mcmc(mcmc_sv_df, x))
plots_combined <- wrap_plots(plots, ncol = 2)
ggsave("output/mcmc_sv.png", plots_combined, width = 7, height = 5, dpi = 500)
