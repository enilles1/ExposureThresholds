# Estimate per-exposure risk

library(truncnorm)
library(ggplot2)
library(dplyr)
library(minpack.lm)
library(patchwork)
library(bbmle)
library(MASS)
library(Matrix)

# Load test negative data
data_in <- read.csv("TND_data.csv")

# Infection risk function
# Edited to give a smoother function
infection_risk_per_exposure <- function(t) {
  1 / (1 + exp(1.5 * (t - 3)))
}

# Scaled logit risk function from Middleton et al
# Simple function for exposure = 0 or 1
scaled_logit_fn <- function(t, beta0, beta1, lambda) {
  lambda / (1 + exp(beta0 *(t - beta1)))
}

# Model with repeat exposures
# (Not identifiable)
infection_risk_given_titer <- function(t, exposures, beta0, beta1, lambda) {
  per_exposure_risk <- scaled_logit_fn(t, beta0, beta1, lambda)
  1 - (1 - per_exposure_risk)^exposures
}

# Model with repeat exposures and scaling term
infection_risk_scaled <- function(t, exposures, exposure_scaling, beta0, beta1, lambda) {
  per_exposure_risk <- scaled_logit_fn(t, beta0, beta1, lambda)
  1 - (1 - per_exposure_risk)^(exposures*exposure_scaling)
}

logit_risk <- function(t, beta0, beta1, lambda) {
  lambda / (1 + exp(beta0 * (t - beta1)))
}




# Plotting function for titer-based distributions
fit_plot <- function() {
  
  # Infection risk vs titres
  df <- data.frame(
    Titer = data_in$s_titer_geom,
    Infection = data_in$PCR_result,
    Transmission_val = data_in$percent_positive
  ) 

  # Swap out values
  titer_density <- density(df$Titer, from = 0.1, to = 5)
  titer_density_df <- data.frame(
    Titer = titer_density$x,
    Density = titer_density$y / max(titer_density$y)
  )
  
  # Fit scaled logit model across all transmission groups
  group_df <- df
  
  # Starting parameter guesses
  start_vals <- list(exposure_scaling = 10, beta0 = 1.5, beta1 = 3, lambda = 0.5)

  infection_loglik <- function(exposure_scaling, beta0, beta1, lambda) {
    p <- 1 - (1 - logit_risk(group_df$Titer, beta0, beta1, lambda))^(group_df$Transmission_val * exposure_scaling)
    p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
    -sum(dbinom(group_df$Infection, size = 1, prob = p, log = TRUE))
  }
  
  
  fit_mle <- mle2(infection_loglik,
                  start = list(exposure_scaling = 0.1, beta0 = 0.1, beta1 = 5, lambda = 0.5),
                  method = "L-BFGS-B",
                  lower = c(exposure_scaling = 0, beta0 = -Inf, beta1 = -Inf, lambda = 0),
                  upper = c(exposure_scaling = Inf, beta0 = Inf, beta1 = Inf, lambda = 1))
  
  # If fitting fails, return empty data frame
  if (is.null(fit)) return(data.frame())
  
  # Create prediction range for single exposure
  titer_seq <- seq(min(group_df$Titer), max(group_df$Titer), length.out = 200)
  
  # Draw 1000 bootstrap samples from the parameter covariance
  cov_matrix <- vcov(fit_mle)
  cov_pd <- as.matrix(nearPD(cov_matrix)$mat)
  
  param_samples <- mvrnorm(n = 1000, mu = coef(fit_mle), Sigma = cov_pd)
  
  # Safe wrapper for prediction
  safe_predict <- function(pars) {
    tryCatch({
      preds <- infection_risk_scaled(
        t = titer_seq,
        exposures = 1,
        exposure_scaling = pars["exposure_scaling"],
        beta0 = pars["beta0"],
        beta1 = pars["beta1"],
        lambda = pars["lambda"]
      )
      if (any(!is.finite(preds))) return(rep(NA, length(titer_seq)))
      preds
    }, error = function(e) rep(NA, length(titer_seq)))
  }
  
  # Apply safe prediction
  pred_matrix <- apply(param_samples, 1, safe_predict)
  
  # Drop columns with NA
  pred_matrix <- pred_matrix[, colSums(is.na(pred_matrix)) == 0]
  
  # Compute 95% prediction intervals
  pred_mean <- rowMeans(pred_matrix)
  pred_lower <- apply(pred_matrix, 1, quantile, probs = 0.025)
  pred_upper <- apply(pred_matrix, 1, quantile, probs = 0.975)
  
  # Output data with prediction intervals
  df_pred <- data.frame(
    Titer = titer_seq,
    RelativeRisk = pred_mean / pred_mean[1],
    Lower = pred_lower / pred_mean[1],
    Upper = pred_upper / pred_mean[1],
    Transmission = "Per-exposure risk"
  )
  
  # Update ggplot with ribbon
  ggplot(df_pred, aes(x = Titer, y = RelativeRisk, color = Transmission)) +
    geom_area(data = titer_density_df, aes(x = Titer, y = Density), inherit.aes = FALSE,
              fill = "lightgray", alpha = 0.4) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#ff7f0e", alpha = 0.2, color = NA) +
    geom_line(size = 1.2) +
    scale_color_manual(values = c("Per-exposure risk" = "#ff7f0e")) +
    scale_fill_manual(values = c("Per-exposure risk" = "#ff7f0e")) +
    labs(x = "Anti-Spike Antibody Titer", y = "Probability of Infection") +
    coord_cartesian(xlim = c(0, 5), ylim = c(0, 1)) +
    theme_light(base_size = 10) +
    theme(
      axis.text.x = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 11),
      legend.title = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.5),
      panel.background = element_blank(),
      legend.position = c(0.98, 0.5),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = "gray"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.5, "cm")
    ) +
    guides(fill = guide_none(), color = guide_legend(title = NULL))

}

# Generate plots for each distribution
g1 <- fit_plot() + theme(plot.tag.position = "topleft")

total_combined <- g1
total_combined

ggsave("fitted_model_data.png", total_combined, width = 8, height = 6, units = "in", dpi = 300)
