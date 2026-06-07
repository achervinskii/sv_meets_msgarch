set.seed(1)

library(ggplot2)
library(stochvol)
library(MSGARCH)
library(future.apply)
library(progressr)
handlers(handler_progress(format = "[:bar] :percent :eta :message"))

# data fetching. not fetched if exists per ticker
source("./R/fetching_data.R")
start_date <- "2016-01-01"
end_date <- "2026-01-01"

ticker <- "AAPL"
filename <- paste0("./data/", ticker, ".csv")

if (!file.exists(filename)) {
  df <- fetch_daily_prices(
    ticker = ticker, start_date = start_date,
    end_date = end_date
  )
  write.csv(df, file = filename)
} else {
  df <- read.csv(filename)
  df$datetime <- as.POSIXct(df$datetime, tz = "UTC")
}

# test run on the whole sample
y <- df[["close"]]
lret <- diff(log(y))
lr <- length(lret)


sv_model <- function(draws = 1000, burnin = 1000, n_chains = 1) {
  list(
    draws = draws,
    burnin = burnin,
    n_chains = n_chains,
    type = "sv"
  )
}

ms_model <- function(draws = 1000, burnin = 1000, nsim = 10000) {
  list(
    ctr = list(
      nburn = burnin,
      nmcmc = draws + burnin,
      nthin = 1
    ),
    nsim = nsim,
    type = "ms"
  )
}

predict_model <- function(model, y, n_ahead) {
  # predicts with one chain always
  if (model$type == "sv") {
    fit <- svsample(
      y,
      draws = model$draws,
      burnin = model$burnin,
      n_chain = model$n_chains,
      quiet = TRUE
    )
    mc_pred <- predict(fit, steps = n_ahead)
    # both are M x h matrices, with row corresponding to the sample
    mc_sigma <- mc_pred$vol[[1]]
    mc_y <- mc_pred$y[[1]]

    result <- list(
      sigma = mc_sigma,
      y = mc_y
    )
  }

  if (model$type == "ms") {
    # create a specification for msgarch
    spec <- CreateSpec(
      variance.spec = list(model = c("sGARCH", "sGARCH")),
      distribution.spec = list(distribution = c("norm", "norm"))
    )
    fit <- FitMCMC(spec = spec, data = y, ctr = model$ctr)
    mc_pred <- predict(
      object = fit,
      nahead = n_ahead,
      do.return.draw = TRUE
    )
    mc_y <- t(mc_pred$draw)
    mc_sigma <- matrix(NA_real_, nrow = length(mc_y), ncol = n_ahead)
    result <- list(
      sigma = mc_sigma,
      y = mc_y
    )
  }

  result
}

loglik <- function(model, test_y, mc_y, mc_sigma) {
  if (model$type == "sv") {
    result <- log(mean(dnorm(test_y, sd = mc_sigma)))
  }
  if (model$type == "ms") {
    bw <- bw.SJ(mc_y) # verify
    result <- log(mean(dnorm((test_y - mc_y) / bw)) / bw)
  }
  result
}


run_backtest_from <- function(model, last_obs_idx, n_ahead, bt_idx = NULL) {
  train_y <- lret[1:last_obs_idx]
  pred <- predict_model(model, train_y, n_ahead)
  pred_y <- pred$y
  pred_sigma <- pred$sigma

  log_scores_i <- numeric(n_ahead)
  var_95_i <- matrix(NA_real_, n_ahead, 2)

  for (h in 1:n_ahead) {
    log_scores_i[h] <- loglik(
      model,
      test_y = lret[last_obs_idx + h],
      mc_y = pred_y[, h],
      mc_sigma = pred_sigma[, h]
    )
    var_95_i[h, ] <- quantile(
      pred_y[, h],
      probs = c(0.025, 0.975)
    )
  }

  list(
    idx = bt_idx,
    log_scores = log_scores_i,
    var_95 = var_95_i
  )
}

# hyperparames for backtesting
n_ahead <- 5
backtest_length <- 500
draws <- 1000
n_chains <- 1
burnin <- 1000
# derived params for backtesting
# number of the data points used for test
n_bt <- backtest_length - n_ahead + 1
last_train_idx_start <- lr - backtest_length - 1
last_train_idx_end <- lr - n_ahead - 1

sv <- sv_model()

plan(multisession, workers = 12)

# preallocation
log_scores_sv <- matrix(
  NA_real_,
  nrow = n_bt,
  ncol = n_ahead
)

var_95_sv <- array(
  NA_real_,
  dim = c(n_bt, n_ahead, 2)
)


xs <- seq_len(n_bt)
with_progress({
  p <- progressor(along = xs)
  results <- future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      run_backtest_from(sv, last_train_idx_start + x - 1, n_ahead, x)
    },
    future.seed = TRUE
  )
})

for (res in results) {
  idx <- res$idx
  log_scores_sv[idx, ] <- res$log_scores
  var_95_sv[idx, , ] <- res$var_95
}

msgarch <- ms_model()

# preallocation
log_scores_ms <- matrix(
  NA_real_,
  nrow = n_bt,
  ncol = n_ahead
)

var_95_ms <- array(
  NA_real_,
  dim = c(n_bt, n_ahead, 2)
)

xs <- seq_len(n_bt)
with_progress({
  p <- progressor(along = xs)
  results <- future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      run_backtest_from(msgarch, last_train_idx_start + x - 1, n_ahead, x)
    },
    future.seed = TRUE
  )
})

for (res in results) {
  idx <- res$idx
  log_scores_ms[idx, ] <- res$log_scores
  var_95_ms[idx, , ] <- res$var_95
}


horizon_to_plot <- 5
test_start_idx <- last_train_idx_start + horizon_to_plot
test_end_idx <- last_train_idx_end + horizon_to_plot

test_data <- lret[test_start_idx:test_end_idx]
df <- data.frame(
  t = seq_along(test_data),
  y = test_data,
  var_lb_sv = var_95_sv[, horizon_to_plot, 1],
  var_ub_sv = var_95_sv[, horizon_to_plot, 2],
  var_lb_ms = var_95_ms[, horizon_to_plot, 1],
  var_ub_ms = var_95_ms[, horizon_to_plot, 2]
)
# plot vars arojnd the return series
plot <- ggplot(df, aes(x = t)) +
  geom_line(aes(y = y), color = "black") +
  geom_ribbon(aes(ymin = var_lb_sv, ymax = var_ub_sv),
    fill = "blue", alpha = 0.1
  ) +
  geom_line(aes(y = var_lb_sv), color = "blue") +
  geom_line(aes(y = var_ub_sv), color = "blue") +
  geom_ribbon(aes(ymin = var_lb_ms, ymax = var_ub_ms),
    fill = "red", alpha = 0.1
  ) +
  geom_line(aes(y = var_lb_ms), color = "red") +
  geom_line(aes(y = var_ub_ms), color = "red") +
  theme_minimal()

plot


# plot the cumulative log score difference
# the n_ahead distance is shown in gradient form black to red
