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
lret_dates <- lret_dates[-1]
ndraws <- 5000

ctr <- list(
  nburn = ndraws,
  nmcmc = ndraws,
  nthin = 1
)
spec <- CreateSpec(
  variance.spec = list(model = c("sGARCH", "sGARCH")),
  distribution.spec = list(distribution = c("norm", "norm"))
)
fit_mcmc <- FitMCMC(spec = spec, data = lret, ctr = ctr)
state <- State(object = fit_mcmc)
summary(fit_mcmc)
# exclude the last smoothed, consider the second regime
statesm <- state$SmoothProb[seq_along(lret), , 2]
dim(statesm)
# collapse the simulations by averaging
probs <- rowMeans(statesm)
# get the 95 quantiles
quantiles <- apply(statesm, 1, FUN = function(x) quantile(x, c(0.025, 0.975)))
# combine in df
df_plot <- data.frame(
  r = lret / (6 * sd(lret)) + 0.5,
  t = lret_dates,
  probs = probs,
  lower = quantiles[1, ],
  upper = quantiles[2, ]
)

cols <- list(
  msgarch = adjustcolor("red"),
  sv = adjustcolor("#007712"),
  hl = adjustcolor("#ff6d00")
)

plot <- ggplot(
  df_plot[, ],
  aes(x = t, y = probs, ymin = lower, ymax = upper)
) +
  # geom_ribbon(aes(fill = ".95 band")) +
  geom_point(aes(x = t, y = r, color = "Return obs."),
    size = 0.6, alpha = 0.7
  ) +
  geom_line(aes(color = "Smoothed prob."), alpha = 0.7) +
  scale_color_manual(values = c("black", cols$msgarch)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
plot

ggsave("output/msgarch_states_1.png", plot, width = 10, height = 5, dpi = 300)
