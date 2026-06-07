set.seed(1)
source("./R/simulation_and_plotting.R")
library(ggplot2)
library(patchwork)

cols <- list(
  msgarch = adjustcolor("red"),
  sv = adjustcolor("#007712"),
  hl = adjustcolor("#ff6d00")
)
alp <- 0.6
lw <- 0.6
ymax <- 15

w <- c(1, 6)
a <- c(0.2, 0.2)
b <- c(0.6, 0.6)
p11 <- 0.95
p22 <- 0.98
p <- matrix(c(
  p11, 1 - p11,
  1 - p22, p22
), nrow = 2)

# number of observation from the tail to select
n <- 300
# length of each of the simulations
len <- 10000

# determine parameters of a matched sv
sv_pars <- match_sv_to_msgarch_2(w, a, b, p, nchains = 10000, burnin = 100)

# simulate SV
asim_sv <- sim_sv(sv_pars$mu, sv_pars$phi, sv_pars$sigma_u, len = len)
# simulate MSGARCH again
asim_msgarch <- sim_msgarch(w, a, b, p, len = len)

y_msgarch <- asim_msgarch$y
s_msgarch <- asim_msgarch$s
y_sv <- asim_sv$y

# compare the last n obs.
y_msgarch_tail <- tail(asim_msgarch$y, n)
s_msgarch_tail <- tail(asim_msgarch$s, n)
y_sv_tail <- tail(asim_sv$y, n)

df <- data.frame(
  t = seq_along(y_msgarch_tail),
  y_msgarch = y_msgarch_tail,
  y_sv = y_sv_tail,
  regime = s_msgarch_tail
)
is_second <- df$regime == 2
r <- rle(is_second)
ends <- cumsum(r$lengths) + 0.5
starts <- ends - r$lengths + 1 - 0.5
idx <- which(r$values)
second_regime_df <- data.frame(
  tmin = starts[idx],
  tmax = ends[idx]
)
plot <- ggplot(df, aes(x = t)) +
  geom_rect(
    data = second_regime_df,
    aes(xmin = tmin, xmax = tmax, ymin = -ymax, ymax = ymax),
    fill = cols$hl,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_line(aes(y = y_msgarch),
    color = cols$msgarch, linewidth = lw, alpha = alp
  ) +
  geom_line(aes(y = y_sv), color = cols$sv, linewidth = lw, alpha = alp) +
  scale_y_continuous(limits = c(-ymax, ymax)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  )

burnin <- 1000
y_dim <- length(y_sv)
stopifnot(y_dim == length(y_msgarch))
# create a plot with densities for the entire sample
df_dens <- data.frame(
  y = c(y_msgarch[burnin:y_dim], y_sv[burnin:y_dim]),
  group = rep(c("msgarch", "sv"), each = (y_dim - burnin + 1))
)

nmu <- mean(df_dens$y)
nsigma <- sd(df_dens$y)
# plot densities
plot_dens <- ggplot(df_dens, aes(x = y, colour = group)) +
  geom_line(stat = "density", alpha = alp, linewidth = lw) +
  geom_function(
    aes(color = "norm"),
    fun = dnorm,
    args = list(mean = nmu, sd = nsigma),
    linewidth = lw,
    linetype = 2,
    alpha = alp,
  ) +
  scale_x_continuous(limits = c(-ymax, ymax)) +
  scale_colour_manual(
    name = NULL,
    values = c(cols$msgarch, cols$sv, "blue"),
    labels = c("MS-GARCH", "SV", "Normal density"),
    breaks = c("msgarch", "sv", "norm")
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = c(0.6, 0.9)
  )

plot_comb <- plot + plot_dens + plot_layout(widths = c(3, 1))

ggsave("output/plot2.png", plot_comb, width = 10, height = 5, dpi = 300)
