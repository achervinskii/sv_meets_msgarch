is_stationary_msgarch <- function(w, a, b, p) {
  n_regimes <- length(w)
  res <- TRUE
  if (!(length(a) == n_regimes && length(b) == n_regimes)) res <- FALSE
  if (!all(dim(p) == c(n_regimes, n_regimes))) res <- FALSE
  if (!all(a + b < 1)) res <- FALSE
  res
}

is_stationary_sv <- function(eta_0, eta_1, sigma_u) {
  res <- TRUE
  if (!all(lengths(list(eta_0, eta_1, sigma_u)) == 1)) res <- FALSE
  if (!(abs(eta_1) < 1)) res <- FALSE
  if (sigma_u <= 0) res <- FALSE
  res
}

get_quantities_msgarch <- function(w, a, b, p) {
  stationary_d_mc <- function(p) {
    # find the stationary probabilities by finding the fixed point
    eig <- eigen(t(p))
    # the closest to one
    idx <- which.min(abs(eig$values - 1))
    eigv <- Re(eig$vectors[, idx])
    eigv / sum(eigv)
  }
  # find stationary distribution
  stationary_d <- stationary_d_mc(p)
  # find the unconditional within-regime variances
  regime_vars <- w / (1 - a - b)
  list(
    stationary_d = stationary_d,
    regime_vars = regime_vars
  )
}

get_quantities_sv <- function(eta_0, eta_1, sigma_u) {
  expected_h <- eta_0 / (1 - eta_1)
  var_h <- sigma_u**2 / (1 - eta_1**2)
  uncond_var <- exp(expected_h)
  list(
    expected_h = expected_h,
    var_h = var_h,
    uncond_var = uncond_var
  )
}

sim_msgarch <- function(w, a, b, p, len = 1000) {
  n_regimes <- length(w)
  regimes <- 1:n_regimes
  # verify that the arguments align
  stopifnot(is_stationary_msgarch(w, a, b, p))
  # get statistics
  quants <- get_quantities_msgarch(w, a, b, p)
  # pre-generated errors
  eps <- rnorm(len)
  # create vectors for the output
  s <- vector(mode = "integer", length = len)
  sigmasq <- vector(mode = "double", length = len)
  y <- vector(mode = "double", length = len)
  s[1] <- sample(regimes, size = 1, prob = quants$stationary_d)
  sigmasq[1] <- quants$regime_vars[s[1]]
  y[1] <- sqrt(sigmasq[1]) * eps[1]
  for (i in 2:len) {
    reg <- sample(regimes, size = 1, prob = p[s[i - 1], ])
    s[i] <- reg
    sigmasq[i] <- w[reg] + a[reg] * sigmasq[i - 1] + b[reg] * y[i - 1]**2
    y[i] <- sqrt(sigmasq[i]) * eps[i]
  }
  list(
    y = y,
    s = s,
    sigmasq = sigmasq,
    stationary_d = quants$stationary_d
  )
}

sim_sv <- function(eta_0, eta_1, sigma_u, len = 1000) {
  stopifnot(is_stationary_sv(eta_0, eta_1, sigma_u))
  # pre-generated errors
  eps <- rnorm(len)
  u <- rnorm(len, 0, sigma_u)
  h <- vector(mode = "double", length = len)
  y <- vector(mode = "double", length = len)
  quants <- get_quantities_sv(eta_0, eta_1, sigma_u)
  # draw the first log-variance from the stationary distribution
  h[1] <- rnorm(1, mean = sqrt(quants$expected_h), sd = sqrt(quants$var_h))
  y[1] <- exp(h[1] / 2) * eps[1]
  for (i in 2:len) {
    h[i] <- eta_0 + eta_1 * h[i - 1] + u[i]
    y[i] <- exp(h[i] / 2) * eps[i]
  }
  list(
    y = y,
    h = h
  )
}

match_sv_to_msgarch_2 <- function(w, a, b, p, nchains = 1000, burnin = 1000) {
  var_t <- vector(mode = "double", length = nchains)
  var_t_plus_1 <- vector(mode = "double", length = nchains)
  for (i in 1:nchains) {
    asim_msgarch <- sim_msgarch(w, a, b, p, burnin + 2)
    var_t[i] <- asim_msgarch$sigmasq[burnin + 1]
    var_t_plus_1[i] <- asim_msgarch$sigmasq[burnin + 2]
  }
  # estimate a linear model
  logvar_t <- log(var_t)
  logvar_t_plus_1 <- log(var_t_plus_1)
  phi <- cov(logvar_t, logvar_t_plus_1) / var(logvar_t_plus_1)
  mu <- mean(logvar_t_plus_1) - phi * mean(logvar_t)
  res <- logvar_t_plus_1 - mu - phi * logvar_t
  sigma_u <- sqrt(sum(res^2) / (nchains - 2))
  list(
    mu = mu,
    phi = phi,
    sigma_u = sigma_u
  )
}
