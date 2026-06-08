set.seed(1)

library(ggplot2)
library(stochvol)
library(MSGARCH)
library(future.apply)
library(progressr)
library(dplyr)
handlers(handler_progress(format = "[:bar] :percent :eta :message"))

# data fetching. not fetched if exists per ticker
source("./R/fetching_data.R")
start_date <- "2016-01-01"
end_date <- "2026-01-01"

ticker <- "^GSPC"

# contsurct the return series
y <- df[["close"]]
lret <- diff(log(y))

# compute realized volatility on hourly data
df_hourly <- fetch_hourly_prices()


sv_model <- function(draws = 10000, burnin = 10000, n_chains = 1) {
  list(
    draws = draws,
    burnin = burnin,
    n_chains = n_chains,
    type = "sv"
  )
}

ms_model <- function(draws = 10000, burnin = 10000) {
  list(
    ctr = list(
      nburn = burnin,
      nmcmc = draws + burnin,
      nthin = 5
    ),
    # if nsim is set then for each posterior draw nsim observations are sampled
    # to me it seems more appripriate to sample one simulation per mcmc path
    # it also possibly causes memory issues, look up the code
    # nsim = nsim,
    type = "msgarch"
  )
}

# TODO: a proper method dispatch is needed
# I have not yet figured out how the dispatch works when parallelized
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

  if (model$type == "msgarch") {
    # create a specification for msgarch
    spec <- CreateSpec(
      variance.spec = list(model = c("sGARCH", "sGARCH")),
      distribution.spec = list(distribution = c("norm", "norm"))
    )
    fit <- FitMCMC(spec = spec, data = y, ctr = model$ctr)
    mc_pred <- predict(
      object = fit,
      nahead = n_ahead,
      do.return.draw = TRUE,
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

loglik_sv <- function(test_y, mc_sigma) {
  log(mean(dnorm(test_y, sd = mc_sigma)))
}

loglik_kde <- function(test_y, mc_y, df, bw) {
  z <- (test_y - mc_y) / bw
  if (df == -1) { # gaussian kde
    result <- log(mean(dnorm(z) / bw))
  } else { # t kde
    result <- log(mean(dt(z, df = df) / bw))
  }
  result
}

run_backtest_from <- function(model, last_obs_idx, n_ahead,
                              bt_idx = NULL, df = NULL) {
  train_y <- lret[1:last_obs_idx]
  pred <- predict_model(model, train_y, n_ahead)
  pred_y <- pred$y
  pred_sigma <- pred$sigma
  bw <- bw.SJ(train_y)

  if (!is.null(df)) {
    # prepare for recieing log score for every df
    log_scores_i <- array(NA_real_, dim = c(n_ahead, length(df)))
  } else {
    log_scores_i <- array(NA_real_, dim = c(n_ahead, 1))
  }

  # VaR alsways has two values and computed from mc y
  var_95_i <- matrix(NA_real_, n_ahead, 2)

  for (h in 1:n_ahead) {
    test_y <- lret[last_obs_idx + h]

    if (!is.null(df)) {
      # do KDEs
      for (df_idx in seq_along(df)) {
        log_scores_i[h, df_idx] <- loglik_kde(
          test_y = test_y,
          mc_y = pred_y[, h],
          df = df[df_idx],
          bw = bw
        )
      }
    } else {
      log_scores_i[h, 1] <- loglik_sv(
        test_y = test_y,
        mc_sigma = pred_sigma[, h]
      )
    }

    # quantile is the same for both models
    var_95_i[h, ] <- quantile(
      pred_y[, h],
      probs = c(0.025, 0.975)
    )
  }

  # the return, assumed to be a worker specific
  list(
    idx = bt_idx,
    log_scores = log_scores_i,
    var_95 = var_95_i
  )
}

# hyperparames for backtesting
n_ahead <- 10
backtest_length <- 500
msgarch_df <- c(-1, 2, 3)
sv_df <- c(-1, 2, 3)
# derived params for backtesting
# number of the data points used for test used for array indexing
n_bt <- backtest_length - n_ahead + 1
# the moving index of the last training observation
n_total <- length(lret)
last_train_idx_start <- n_total - backtest_length - 1
last_train_idx_end <- n_total - n_ahead - 1
# parallelize
plan(multisession, workers = 6)

# instance of sv
sv <- sv_model()

# preallocation
log_scores_sv <- array(
  NA_real_,
  dim = c(n_bt, n_ahead, length(sv_df))
)

var_95_sv <- array(
  NA_real_,
  dim = c(n_bt, n_ahead, 2)
)

# split the training sets between the workers
# the packages don't support online learning anyways
xs <- seq_len(n_bt)
# with a progress bar
with_progress({
  p <- progressor(along = xs)
  results_sv <- future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      run_backtest_from(sv, last_train_idx_start + x - 1,
        n_ahead,
        bt_idx = x,
        df = sv_df
      )
    },
    # TODO: verify this argument. It is pasted from a warning.
    future.seed = TRUE
  )
})

# gather the results from workers
for (res in results_sv) {
  idx <- res$idx
  log_scores_sv[idx, , ] <- res$log_scores
  var_95_sv[idx, , ] <- res$var_95
}

plan(multisession, workers = 6)
# repeat for ms garch
msgarch <- ms_model()

# preallocation
log_scores_ms <- array(
  NA_real_,
  dim = c(n_bt, n_ahead, length(msgarch_df))
)

var_95_ms <- array(
  NA_real_,
  dim = c(n_bt, n_ahead, 2)
)

xs <- seq_len(n_bt)
with_progress({
  p <- progressor(along = xs)
  results_ms <- future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      run_backtest_from(msgarch, last_train_idx_start + x - 1,
        n_ahead,
        bt_idx = x,
        df = msgarch_df
      )
    },
    future.seed = TRUE
  )
})

for (res in results_ms) {
  idx <- res$idx
  log_scores_ms[idx, , ] <- res$log_scores
  var_95_ms[idx, , ] <- res$var_95
}

# gather everything into a dataframe


# plot of the VaR with 5 day horizon
horizon_to_plot <- 5
test_start_idx <- last_train_idx_start + horizon_to_plot
test_end_idx <- last_train_idx_end + horizon_to_plot
test_data <- lret[test_start_idx:test_end_idx]

var_df <- data.frame(
  t = seq_along(test_data),
  y = test_data,
  var_lb_sv = var_95_sv[, horizon_to_plot, 1],
  var_ub_sv = var_95_sv[, horizon_to_plot, 2],
  var_lb_ms = var_95_ms[, horizon_to_plot, 1],
  var_ub_ms = var_95_ms[, horizon_to_plot, 2]
)

# plot vars arojnd the return series
plot <- ggplot(var_df, aes(x = t)) +
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

ggsave("output/value_at_risk_ms_sv.png",
  plot,
  width = 10, height = 5, dpi = 300
)


# plot the cumulative log score difference
# the n_ahead distance is shown in gradient form black to red
m <- log_scores_ms
s <- log_scores_sv

df_m <- as.data.frame(as.table(m))
df_s <- as.data.frame(as.table(s))
df_m$model <- "msgarch"
df_s$model <- "sv"
df_long <- rbind(df_m, df_s)


names(df_long) <- c("bt_idx", "h", "df_kde", "loglik", "model")
# convert the factors to integers
df_long$bt_idx <- as.integer(df_long$bt_idx)
df_long$h <- as.integer(df_long$h)
df_long$df_kde <- as.integer(df_long$df_kde)
# set df_kde to NULL for sv
# remap the numeric factors back to the actual DoFs
msgarch_idx <- df_long$model == "msgarch"
sv_idx <- df_long$model == "sv"

df_long$df_kde <- ifelse(
  df_long$model == "msgarch",
  msgarch_df[df_long$df_kde],
  sv_df[df_long$df_kde]
)

df_long$df_kde <- factor(
  df_long$df_kde,
  levels = sort(unique(df_long$df_kde)),
  labels = paste0("v = ", sort(unique(df_long$df_kde)))
)
head(df_long)

# align to the prediction period
# subtract 1 since t is indexed from 1, add 1 since we use the return series
df_long$relative_time <- df_long$bt_idx - 1 + df_long$h + 1
df_long$absolute_time <- df_long$relative_time + last_train_idx_start
df_long$datetime <- df$datetime[df_long$absolute_time]

# save df_long, to avoid recomputing when debuging
file_df_long <- "./output/df_long.rda"
save(df_long, file = file_df_long)
# load(file_df_long)

# compute the cumulative predictions within each tuple (model, horizon, and DoF)
# cumulation only across the time dimension
df_cumulative <- df_long |>
  group_by(model, h, df_kde) |>
  mutate(cum_loglik = order_by(relative_time, cumsum(loglik))) |>
  ungroup()

# compute difference between different models
df_diff <- df_cumulative |>
  group_by(bt_idx, h, df_kde, datetime) |>
  summarize(
    cum_diff = cum_loglik[model == "msgarch"] -
      cum_loglik[model == "sv"],
    .groups = "drop"
  )

# plot
diff_plot <- ggplot(df_diff,
  aes(x = datetime, y = cum_diff, color = h, group = h)
) +
  geom_line(aes(color = h)) +
  scale_color_gradient(low = "green", high = "red") +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-15, 15)) +
  theme_minimal() +
  facet_grid(rows = vars(df_kde))

ggsave("output/cum_loglik_ms_minus_sv_across_kde.png",
  diff_plot,
  width = 10, height = 12, dpi = 300
)
