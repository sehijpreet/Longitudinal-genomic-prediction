rm(list=ls())

library(dplyr)

setwd('')

list.files()


# Load libraries
library(ggplot2)

# Define Legendre function
generate_legendre <- function(x, degree) {
  legendre_basis <- data.frame(P0 = rep(1, length(x)))
  if (degree >= 1) legendre_basis$P1 <- x
  if (degree >= 2) legendre_basis$P2 <- 0.5 * (3 * x^2 - 1)
  if (degree >= 3) legendre_basis$P3 <- 0.5 * (5 * x^3 - 3 * x)
  if (degree >= 4) legendre_basis$P4 <- (35 * x^4 - 30 * x^2 + 3) / 8
  return(legendre_basis)
}

# Compute polynomial values
compute_polynomial_values <- function(coefficients, basis) {
  coefficients[1] * basis$P0 +
    coefficients[2] * basis$P1 +
    coefficients[3] * basis$P2 +
    coefficients[4] * basis$P3
}
Y <- read.csv('leg_coef.csv')
head(Y)



Y <- read.csv('leg_coef.csv')
PC11M1 <- read.csv('Predictions_CV11_M1.csv')
PC11M2 <- read.csv('Predictions_CV11_M2.csv')
PC10M1 <- read.csv('Predictions_CV10_M1.csv')
PC10M2 <- read.csv('Predictions_CV10_M2.csv')
PC9M1  <- read.csv('Predictions_CV9_M1.csv')
PC9M2  <- read.csv('Predictions_CV9_M2.csv')

# Prepare lookup for all datasets
pred_list <- list(
  C11_M1 = PC11M1,
  C11_M2 = PC11M2,
  C10_M1 = PC10M1,
  C10_M2 = PC10M2,
  C9_M1  = PC9M1,
  C9_M2  = PC9M2
)


x_values <- seq(-1, 1, length.out = 8)
legendre_basis <- generate_legendre(x_values, degree = 3)

# Loop over each environment-model pair -------------------------------------

all_results <- list()

for (name in names(pred_list)) {
  cat("Processing:", name, "\n")
  
  # Split environment and model info
  env <- strsplit(name, "_")[[1]][1]
  model <- strsplit(name, "_")[[1]][2]
  
  # Filter coefficients
  YN <- as.matrix(Y[Y$Source == env, c("V1", "V2", "V3", "V4")])
  PN <- as.matrix(pred_list[[name]][Y$Source == env, c("Trait_1", "Trait_2", "Trait_3", "Trait_4")])
  
  # Compute polynomial values
  y_obs_list <- lapply(1:nrow(YN), function(i) compute_polynomial_values(as.numeric(YN[i, ]), legendre_basis))
  y_pred_list <- lapply(1:nrow(PN), function(i) compute_polynomial_values(as.numeric(PN[i, ]), legendre_basis))
  
  # Convert to matrices
  O_leg <- do.call(cbind, y_obs_list)
  P_leg <- do.call(cbind, y_pred_list)
  
  # Compute correlations across weeks
  corr_values <- sapply(1:nrow(O_leg), function(i) cor(O_leg[i, ], P_leg[i, ], use = "complete.obs"))
  
  # Store results in a dataframe
  df <- data.frame(
    Environment = env,
    Model = model,
    Week = 1:length(corr_values),
    Correlation = corr_values
  )
  
  all_results[[name]] <- df
}

# Combine and save -----------------------------------------------------------
final_df <- do.call(rbind, all_results)


# Create a combined Environment_Model label
final_df <- final_df %>%
  mutate(Env_Model = paste(Environment, Model, sep = "_"))

# Pivot to wide format: weeks as rows, Env_Model as columns
wide_df <- final_df %>%
  select(Week, Env_Model, Correlation) %>%
  pivot_wider(
    names_from = Env_Model,
    values_from = Correlation
  ) %>%
  arrange(Week)

# View result
print(wide_df)
write.csv(wide_df, "Correlation_Results.csv", row.names = FALSE)

