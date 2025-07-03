##################################################################
### Generating figures for Bayesian spatio-temporal SAE model  ###
### Author: Elliot S. Shannon                                  ###
### Email: shann125@msu.edu                                    ###
### Affiliation: Michigan State University                     ###
### Date Created: July 02, 2025                                ###
### Date Modified: July 03, 2025                               ###
### R version: 4.3.1 (2023-06-16)                              ###
##################################################################

rm(list = ls())

######################
### Load Libraries ###
######################

library(sf) # sf_1.0-19
library(ggplot2) # ggplot2_3.5.1
library(dplyr) # dplyr_1.1.4
library(stringr) # stringr_1.5.1
library(scico) # scico_1.5.0
library(ggnewscale) # ggnewscale_0.5.0

#################################
### Function to get quantiles ###
#################################

quant <- function(x) {
  c(
    quantile(x, prob = c(0.5, 0.025, 0.975), na.rm = TRUE),
    mean(x, na.rm = TRUE),
    sd(x, na.rm = TRUE)
  )
}

########################
### Load in the data ###
########################

data <- readRDS("data.rds")

########################
### Prepare the data ###
########################

counties <- data$counties
plots <- data$plots
y <- plots$y # response
x <- plots$x # covariates
W <- data$W # adjacency matrix

years <- c(2008:2021)

T <- length(years)
J <- nrow(W)
N <- length(y)

#################################
### Load in posterior samples ###
#################################

samples <- readRDS("samples.rds")
n.samples <- dim(samples$mu.samples)[2]

###############################
### Get posterior quantiles ###
###############################

mu_quants <- apply(samples$mu.samples, 1, FUN = quant) %>%
  round(digits = 3)
beta_quants <- apply(samples$beta.samples, 1, FUN = quant) %>%
  round(digits = 3)
eta_quants <- apply(samples$eta.samples, 1, FUN = quant) %>%
  round(digits = 3)
u_quants <- apply(samples$u.samples, 1, FUN = quant) %>%
  round(digits = 3)
sigma_sq_quants <- apply(samples$sigma_sq.samples, 1, FUN = quant) %>%
  round(digits = 3)
Sigma_beta_quants <- apply(samples$Sigma_beta.samples, c(1, 2), FUN = quant) %>%
  round(digits = 3)
beta_0_quants <- apply(samples$beta_0.samples, 1, FUN = quant) %>%
  round(digits = 3)
tau_sq_eta_quants <- quant(samples$tau_sq_eta.samples) %>%
  round(digits = 3)
rho_eta_quants <- quant(samples$rho_eta.samples) %>%
  round(digits = 3)
tau_sq_w_quants <- apply(samples$tau_sq_w.samples, 1, FUN = quant) %>%
  round(digits = 3)
rho_w_quants <- quant(samples$rho_w.samples) %>%
  round(digits = 5)

#####################################
### Space-varying coefficient map ###
#####################################

svc_fig <- ggplot() +
  geom_sf(
    data = counties,
    lwd = 0.1,
    color = "white",
    aes(fill = eta_quants[4, ])
  ) +
  scale_fill_scico(NULL, palette = "navia", direction = -1) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  ggtitle("Posterior Mean\nSpace-Varying Coefficient") +
  theme_bw() +
  theme(
    text = element_text(size = 12, family = "Serif"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.key.width = unit(1, "cm"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

svc_fig

ggsave(
  plot = svc_fig,
  filename = "./figs/svc.png",
  device = "png",
  width = 4,
  height = 7,
  units = "in"
)

################################
### Posterior mean of mu map ###
################################

mu_posts_df <- cbind(
  counties,
  "year" = rep(years[1], times = J),
  "mu_mean" = mu_quants[4, 1:J],
  "mu_sd" = mu_quants[5, 1:J]
)

for (t in 2:T) {
  mu_posts_df <- rbind(
    mu_posts_df,
    cbind(
      counties,
      "year" = rep(years[t], times = J),
      "mu_mean" = mu_quants[4, ((t - 1) * J + 1):(t * J)],
      "mu_sd" = mu_quants[5, ((t - 1) * J + 1):(t * J)]
    )
  )
}

mu_fig <- ggplot() +
  geom_sf(
    data = mu_posts_df,
    lwd = 0.1,
    color = "white",
    aes(fill = mu_mean)
  ) +
  facet_wrap(~year, nrow = 2) +
  scale_fill_scico("Carbon (Mg/ha)", palette = "bamako", direction = -1) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  theme_bw() +
  ggtitle("Posterior Mean Carbon Density") +
  theme(
    text = element_text(size = 12, family = "Serif"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.key.width = unit(1, "cm"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

mu_fig

ggsave(
  plot = mu_fig,
  filename = "./figs/mu_mean.png",
  device = "png",
  width = 10,
  height = 7,
  units = "in"
)

##############################
### Posterior sd of mu map ###
##############################

sd_fig <- ggplot() +
  geom_sf(
    data = mu_posts_df,
    lwd = 0.1,
    color = "white",
    aes(fill = mu_sd)
  ) +
  facet_wrap(~year, nrow = 2) +
  scale_fill_scico(
    "Carbon (Mg/ha)",
    palette = "lipari",
    direction = -1
  ) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  theme_bw() +
  ggtitle("Posterior SD Carbon Density") +
  theme(
    text = element_text(size = 12, family = "Serif"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.key.width = unit(1, "cm"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

sd_fig

ggsave(
  plot = sd_fig,
  filename = "./figs/mu_sd.png",
  device = "png",
  width = 10,
  height = 7,
  units = "in"
)

###############################
### Posterior mean of u map ###
###############################

u_posts_df <- cbind(
  counties,
  "year" = rep(years[1], times = J),
  "u_mean" = u_quants[4, 1:J]
)

for (t in 2:T) {
  u_posts_df <- rbind(
    u_posts_df,
    cbind(
      counties,
      "year" = rep(years[t], times = J),
      "u_mean" = u_quants[4, ((t - 1) * J + 1):(t * J)]
    )
  )
}

u_fig <- ggplot() +
  geom_sf(
    data = u_posts_df,
    lwd = 0.1,
    color = "white",
    aes(fill = u_mean)
  ) +
  facet_wrap(~year, nrow = 2) +
  scale_fill_scico(
    NULL,
    palette = "cork",
    direction = 1,
    limits = range(-max(abs(u_quants[4, ])), max(abs(u_quants[4, ])))
  ) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  ggtitle("Posterior Mean Spatio-temporal Intercept") +
  theme_bw() +
  theme(
    text = element_text(size = 12, family = "Serif"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.key.width = unit(1, "cm"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

u_fig

ggsave(
  plot = u_fig,
  filename = "./figs/u.png",
  device = "png",
  width = 10,
  height = 7,
  units = "in"
)

#######################################
### Posterior mean of mu trends map ###
#######################################

mu.samples <- array(NA, dim = c(J, T, n.samples))
for (i in 1:n.samples) {
  mu.samples[,, i] <- matrix(samples$mu.samples[, i], nrow = J, ncol = T)
}

trend_coef <- apply(
  mu.samples,
  c(1, 3),
  function(row) return(coefficients(lm(row ~ years))[2])
)

trend_quants <- apply(trend_coef, 1, FUN = quant)

which_sig <- which(trend_quants[2, ] * trend_quants[3, ] > 0)

trend_df <- cbind(
  counties,
  "trend" = trend_quants[4, ],
  "sig_trend" = rep(0, times = J)
)
trend_df[which_sig, "sig_trend"] <- trend_quants[4, which_sig]

trend_fig <- ggplot() +
  geom_sf(
    data = trend_df,
    lwd = 0.1,
    color = "white",
    aes(fill = trend)
  ) +
  scale_fill_scico(
    "Carbon (Mg/ha/year)",
    palette = "vik",
    direction = -1,
    limits = range(-max(abs(trend_df$trend)), max(abs(trend_df$trend)))
  ) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  ggtitle("Posterior Mean\nCarbon Density Trends") +
  theme_bw() +
  theme(
    text = element_text(size = 12, family = "Serif"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.key.width = unit(1, "cm"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

trend_fig

ggsave(
  plot = trend_fig,
  filename = "./figs/trend.png",
  device = "png",
  width = 4,
  height = 7,
  units = "in"
)

###################################################
### Posterior mean of significant mu trends map ###
###################################################

sig_trend_fig <- ggplot() +
  geom_sf(
    data = trend_df,
    lwd = 0.1,
    color = "white",
    aes(fill = sig_trend)
  ) +
  scale_fill_scico(
    "Carbon (Mg/ha/year)",
    palette = "vik",
    direction = -1,
    limits = range(-max(abs(trend_df$trend)), max(abs(trend_df$trend))),
    na.value = NA
  ) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  ggtitle("Significant Posterior Mean\nCarbon Density Trends") +
  theme_bw() +
  theme(
    text = element_text(size = 12, family = "Serif"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.key.width = unit(1, "cm"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

sig_trend_fig

ggsave(
  plot = sig_trend_fig,
  filename = "./figs/sig_trend.png",
  device = "png",
  width = 4,
  height = 7,
  units = "in"
)

#######################################
### County-level comparison figures ###
#######################################

direct <- plots %>%
  group_by(t, j) %>%
  summarize("direct_mean" = mean(y), "direct_se" = sd(y) / sqrt(n()), "n" = n())

t <- qt(1 - 0.05 / 2, direct$n - 1)
direct[which(direct$direct_se == 0), c("direct_se")] <- NA
direct$lower <- direct$direct_mean - t * direct$direct_se
direct$upper <- direct$direct_mean + t * direct$direct_se

direct_mean <- matrix(NA, nrow = J, ncol = T)
direct_mean[cbind(direct$j, direct$t)] <- direct$direct_mean
direct_se <- matrix(NA, nrow = J, ncol = T)
direct_se[cbind(direct$j, direct$t)] <- direct$direct_se
direct_lower <- matrix(NA, nrow = J, ncol = T)
direct_lower[cbind(direct$j, direct$t)] <- direct$lower
direct_upper <- matrix(NA, nrow = J, ncol = T)
direct_upper[cbind(direct$j, direct$t)] <- direct$upper

mu_mean <- matrix(mu_quants[4, ], nrow = J, ncol = T)
mu_lower <- matrix(mu_quants[2, ], nrow = J, ncol = T)
mu_upper <- matrix(mu_quants[3, ], nrow = J, ncol = T)

x_full <- data$x_full
x_mat <- matrix(NA, nrow = J, ncol = T)
x_mat[cbind(x_full$j, x_full$t)] <- x_full$x

n_jt <- plots %>% group_by(t, j) %>% summarize("n" = n())
n_jt_mat <- matrix(0, nrow = J, ncol = T)
n_jt_mat[cbind(n_jt$j, n_jt$t)] <- n_jt$n

for (i in 1:J) {
  county_fig <- ggplot() +
    geom_point(
      aes(x = years - 0.1, y = direct_mean[i, ], color = "Direct"),
      size = 1
    ) +
    geom_errorbar(aes(
      x = years - 0.1,
      ymin = direct_lower[i, ],
      ymax = direct_upper[i, ],
      width = 0.1,
      color = "Direct"
    )) +
    geom_point(
      aes(x = years + 0.1, y = mu_mean[i, ], color = "Model"),
      size = 1
    ) +
    geom_errorbar(aes(
      x = years + 0.1,
      ymin = mu_lower[i, ],
      ymax = mu_upper[i, ],
      width = 0.1,
      color = "Model"
    )) +
    scale_color_manual(
      "Estimator",
      values = c("Direct" = "#E54E21", "Model" = "#0A9F9D")
    ) +
    new_scale_color() +
    geom_line(
      aes(
        x = years,
        y = (x_mat[i, ] / 100) *
          (max(
            mu_upper[i, ],
            direct_upper[i, ],
            mu_mean[i, ],
            direct_mean[i, ],
            na.rm = TRUE
          ) -
            min(
              mu_lower[i, ],
              direct_lower[i, ],
              mu_mean[i, ],
              direct_mean[i, ],
              na.rm = TRUE
            )) +
          min(
            mu_lower[i, ],
            direct_lower[i, ],
            mu_mean[i, ],
            direct_mean[i, ],
            na.rm = TRUE
          ),
        color = ""
      ),
      size = 1,
      alpha = 0.5,
      group = 1
    ) +
    scale_color_manual("TCC", values = c("#117733"), breaks = "") +
    geom_vline(xintercept = c(2007.5, years + 0.5), alpha = 0.1) +
    geom_text(aes(
      x = years,
      y = max(
        mu_upper[i, ],
        direct_upper[i, ],
        mu_mean[i, ],
        direct_mean[i, ],
        na.rm = TRUE
      ) +
        0.05 *
          (range(
            mu_upper[i, ],
            mu_lower[i, ],
            direct_upper[i, ],
            direct_lower[i, ],
            mu_mean[i, ],
            direct_mean[i, ],
            na.rm = TRUE
          )[2] -
            range(
              mu_upper[i, ],
              mu_lower[i, ],
              direct_upper[i, ],
              direct_lower[i, ],
              mu_mean[i, ],
              direct_mean[i, ],
              na.rm = TRUE
            )[1]),
      label = as.character(n_jt_mat[i, ]),
      family = "Serif"
    )) +
    scale_y_continuous(
      name = "Carbon (Mg/ha)",
      sec.axis = sec_axis(
        ~ (. -
          min(
            mu_lower[i, ],
            direct_lower[i, ],
            mu_mean[i, ],
            direct_mean[i, ],
            na.rm = TRUE
          )) /
          (max(
            mu_upper[i, ],
            direct_upper[i, ],
            mu_mean[i, ],
            direct_mean[i, ],
            na.rm = TRUE
          ) -
            min(
              mu_lower[i, ],
              direct_lower[i, ],
              mu_mean[i, ],
              direct_mean[i, ],
              na.rm = TRUE
            )) *
          100,
        name = "TCC (%)",
        breaks = c(0, 25, 50, 75, 100)
      )
    ) +
    scale_x_continuous(breaks = years) +
    theme_bw() +
    labs(
      x = element_blank(),
      title = paste0(counties$COUNTYNM[i], ", ", counties$STATENM[i])
    ) +
    theme(
      text = element_text(size = 14, family = "Serif"),
      panel.grid.major.x = element_blank(),
      legend.title.position = "top",
      legend.title = element_text(hjust = 0.5),
      legend.position = 'bottom',
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  file.lab <- paste0(
    counties$COUNTYNM[i],
    "_",
    counties$STATENM[i]
  )
  file.lab <- str_replace_all(file.lab, " ", "_") # For names with spaces.
  file.lab <- gsub("[[:punct:]]", "_", file.lab)

  ggsave(
    plot = county_fig,
    filename = paste0("./figs/counties/", file.lab, ".png"),
    device = "png",
    width = 8,
    height = 5,
    units = "in"
  )
}
