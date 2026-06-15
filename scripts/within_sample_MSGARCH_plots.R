# The goal is to create plots for the in-sample bayesian fit
# The states are of interest, as they are the basis of the forecasting ability in MSGARCH

library(ggplot2)
library(MSGARCH)
library(dplyr)
library(patchwork)

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
ndraws <- 5000

spec <- CreateSpec(
  variance.spec = list(model = c("sGARCH", "sGARCH")),
  distribution.spec = list(distribution = c("norm", "norm"))
)

ctr <- list(
  nburn = ndraws,
  nmcmc = ndraws,
  nthin = 1
)
fit_mcmc <- FitMCMC(spec = spec, data = lret, ctr = ctr)
summary(fit_mcmc)
mcmc_par <- fit_mcmc$par
avg_par <- colMeans(mcmc_par)
# state <- State(object = spec, par = avg_par, data = lret)
state <- State(fit_mcmc)
# exclude the last smoothed, consider the second regime
statesm <- state$SmoothProb[1:dim(state$SmoothProb)[1]-1, , 2]
# collapse the simulations by averaging
if (!is.null(dim(statesm))){
  probs <- rowMeans(statesm)
} else {
  probs <- statesm
}

# probs <- state$SmoothProb[seq_along(lret),1,2]
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
  geom_point(aes(x = t, y = r, color = "r_t"),
    size = 0.6, alpha = 0.7
  ) +
  geom_line(aes(color = "P(S_t = 2)"), alpha = 0.7) +
  scale_color_manual(values = c(
    "r_t" = "black", 
    "P(S_t = 2)" = cols$msgarch)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
plot

ggsave("output/msgarch_states_1.png", plot, width = 10, height = 5, dpi = 300)
