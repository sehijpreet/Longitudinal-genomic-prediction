rm(list=ls())
library(dplyr)

setwd('')

list.files()

data <- read.csv('blue_fil.csv')
head(data)
table(data$Season)

seasons <- unique(data$Season)

for (i in seq_along(seasons)) {
  assign(paste0("data", i), subset(data, Season == seasons[i]))
}

data <- subset(data, TraitDate <= 12)
data1[1:20,]

genotypes1 <- unique(data1$Genotype)
genotypes2 <- unique(data2$Genotype)
genotypes3 <- unique(data3$Genotype)
genotypes4 <- unique(data4$Genotype)
genotypes5 <- unique(data5$Genotype)
genotypes6 <- unique(data6$Genotype)
genotypes7 <- unique(data7$Genotype)
genotypes8 <- unique(data8$Genotype)
genotypes9 <- unique(data9$Genotype)
genotypes10 <- unique(data10$Genotype)
genotypes11 <- unique(data11$Genotype)


library(ggplot2)
ggplot(subset(data9,  Genotype %in% unique(data9$Genotype)[1:5]), aes(x = TraitDate, y = cum_emmean, group = Genotype, color = Genotype)) + 
  geom_line() + 
  theme_minimal() +
  labs(title = "Cumulative Value Over Time", x = "Trait Date", y = "Cumulative Value")

library(dplyr)
library(tidyr)

data1 <- subset(data1, TraitDate <= 10)
data2 <- subset(data2, TraitDate <= 10)
data3 <- subset(data3, TraitDate <= 10)
data4 <- subset(data4, TraitDate <= 10)
data5 <- subset(data5, TraitDate <= 9)
data6 <- subset(data6, TraitDate <= 9)
data7 <- subset(data7, TraitDate <= 9)
data8 <- subset(data8, TraitDate <= 9)
data9 <- subset(data9, TraitDate <= 9)
data10 <- subset(data10, TraitDate <= 8)
data11 <- subset(data11, TraitDate <= 8)

genotypes1 <- unique(data1$Genotype)
genotypes2 <- unique(data2$Genotype)
genotypes3 <- unique(data3$Genotype)
genotypes4 <- unique(data4$Genotype)
genotypes5 <- unique(data5$Genotype)
genotypes6 <- unique(data6$Genotype)
genotypes7 <- unique(data7$Genotype)
genotypes8 <- unique(data8$Genotype)
genotypes9 <- unique(data9$Genotype)
genotypes10 <- unique(data10$Genotype)
genotypes11 <- unique(data11$Genotype)



# Define the generate_legendre function

generate_legendre <- function(x, degree) {
  legendre_basis <- data.frame(P0 = rep(1, length(x)))
  
  if (degree >= 1) 
  {
    legendre_basis$P1 <- x
  }
  
  if (degree >= 2) 
  {
    legendre_basis$P2 <- 0.5 * (3 * x^2 - 1)
  }
  
  if (degree >= 3)
  {
    legendre_basis$P3 <- 0.5 * (5 * x^3 - 3 * x)
  }
  
  if (degree >= 4) 
  {
    legendre_basis$P4 <- (35 * x^4 - 30 * x^2 + 3) / 8
  }
  if (degree >= 5) 
  {
    legendre_basis$P5 <- (63 * x^5 - 70 * x^3 + 15 * x) / 8
  }
  return(legendre_basis)
}


coefficients_matrix1 <- matrix(NA, nrow = length(genotypes1), ncol = 4)
coefficients_matrix2 <- matrix(NA, nrow = length(genotypes2), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix3 <- matrix(NA, nrow = length(genotypes3), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix4 <- matrix(NA, nrow = length(genotypes4), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix5 <- matrix(NA, nrow = length(genotypes5), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix6 <- matrix(NA, nrow = length(genotypes6), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix7 <- matrix(NA, nrow = length(genotypes7), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix8 <- matrix(NA, nrow = length(genotypes8), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix9 <- matrix(NA, nrow = length(genotypes9), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix10 <- matrix(NA, nrow = length(genotypes10), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3
coefficients_matrix11 <- matrix(NA, nrow = length(genotypes11), ncol = 4)  # Adjusted size to 4 for P0, P1, P2, P3


mse_per_x1 <- vector("list", length(genotypes1))
mse_per_x2 <- vector("list", length(genotypes2))
mse_per_x3 <- vector("list", length(genotypes3))
mse_per_x4 <- vector("list", length(genotypes4))
mse_per_x5 <- vector("list", length(genotypes5))
mse_per_x6 <- vector("list", length(genotypes6))
mse_per_x7 <- vector("list", length(genotypes7))
mse_per_x8 <- vector("list", length(genotypes8))
mse_per_x9 <- vector("list", length(genotypes9))
mse_per_x10 <- vector("list", length(genotypes10))
mse_per_x11 <- vector("list", length(genotypes11))


# ------------------ AIC Storage Lists ------------------
aic_per_x1  <- vector("list", length(genotypes1))
aic_per_x2  <- vector("list", length(genotypes2))
aic_per_x3  <- vector("list", length(genotypes3))
aic_per_x4  <- vector("list", length(genotypes4))
aic_per_x5  <- vector("list", length(genotypes5))
aic_per_x6  <- vector("list", length(genotypes6))
aic_per_x7  <- vector("list", length(genotypes7))
aic_per_x8  <- vector("list", length(genotypes8))
aic_per_x9  <- vector("list", length(genotypes9))
aic_per_x10 <- vector("list", length(genotypes10))
aic_per_x11 <- vector("list", length(genotypes11))

# ------------------ BIC Storage Lists ------------------
bic_per_x1  <- vector("list", length(genotypes1))
bic_per_x2  <- vector("list", length(genotypes2))
bic_per_x3  <- vector("list", length(genotypes3))
bic_per_x4  <- vector("list", length(genotypes4))
bic_per_x5  <- vector("list", length(genotypes5))
bic_per_x6  <- vector("list", length(genotypes6))
bic_per_x7  <- vector("list", length(genotypes7))
bic_per_x8  <- vector("list", length(genotypes8))
bic_per_x9  <- vector("list", length(genotypes9))
bic_per_x10 <- vector("list", length(genotypes10))
bic_per_x11 <- vector("list", length(genotypes11))



# Function to scale to [-1, 1]
scale_to_unit <- function(x) {
  (2 * (x - min(x)) / (max(x) - min(x))) - 1
}

# Function to reverse the scaling (if needed)
rescale_from_unit <- function(x_scaled, original_x) {
  min_x <- min(original_x)
  max_x <- max(original_x)
  ((x_scaled + 1) / 2) * (max_x - min_x) + min_x
}




fit_legendre_model <- function(data, genotypes, degree = 3, color_set,
                               season_label = "", plot_fit = TRUE) {
  
  # Initialize storage
  plot_initialized <- FALSE
  coef_matrix <- matrix(NA, nrow = length(genotypes), ncol = degree + 1)
  mse_list <- vector("list", length(genotypes))
  aic_values <- numeric(length(genotypes))
  bic_values <- numeric(length(genotypes))
  
  # Setup for plotting range
  all_trait_dates <- sort(unique(data$TraitDate))
  y_range <- range(data$cum_emmean, na.rm = TRUE)
  
  # Loop through genotypes
  for (i in seq_along(genotypes)) {
    genotype <- genotypes[i]
    subset_data <- data[data$Genotype == genotype & data$TraitDate %in% all_trait_dates, ]
    
    x <- subset_data$TraitDate
    x_scaled <- scale_to_unit(x)
    y <- subset_data$cum_emmean
    
    # Generate Legendre basis
    legendre_basis <- generate_legendre(x_scaled, degree)
    
    # Fit model (no intercept)
    model <- lm(y ~ . - 1, data = legendre_basis)
    
    # Store coefficients and metrics
    coef_matrix[i, ] <- coef(model)
    y_pred <- predict(model, newdata = legendre_basis)
    
    squared_errors <- (y_pred - y)^2
    
    mse_list[[i]] <- data.frame(
      Genotype = genotype,
      TraitDate = x,
      TD_scaled = x_scaled,
      SquaredError = squared_errors
    )
    
    # Store model AIC and BIC
    aic_values[i] <- AIC(model)
    bic_values[i] <- BIC(model)
    
    # Plot fits (if enabled)
    if (plot_fit) {
      if (!plot_initialized) {
        plot(x_scaled, y_pred, col = color_set[i], xlab = "Scaled TraitDate",
             ylab = "cum_emmean",
             main = paste("Legendre Polynomial Fitting", season_label),
             ylim = y_range)
        plot_initialized <- TRUE
      }
      lines(x_scaled, y_pred, col = color_set[i], lwd = 2)
    }
  }
  
  # Combine squared errors from all genotypes
  mse_df <- do.call(rbind, mse_list)
  
  # Weighted combined MSE across all genotypes
  total_obs <- sum(sapply(genotypes, function(g) nrow(data[data$Genotype == g, ])))
  combined_mse_per_trait_date <- aggregate(
    SquaredError ~ TraitDate,
    data = mse_df,
    FUN = function(x) sum(x) / total_obs
  )
  combined_mse <- sum(mse_df$SquaredError) / total_obs
  
  # Weighted AIC & BIC (by number of observations per genotype)
  obs_per_geno <- sapply(genotypes, function(g) nrow(data[data$Genotype == g, ]))
  combined_aic <- sum(aic_values * obs_per_geno) / sum(obs_per_geno)
  combined_bic <- sum(bic_values * obs_per_geno) / sum(obs_per_geno)
  
  # Combine per-genotype summary
  metrics_df <- data.frame(
    Genotype = genotypes,
    AIC = aic_values,
    BIC = bic_values,
    MSE = sapply(mse_list, function(df) mean(df$SquaredError))
  )
  
  # Return all results
  list(
    coefficients = coef_matrix,
    metrics = metrics_df,
    mse_df = mse_df,
    combined_mse_per_trait_date = combined_mse_per_trait_date,
    combined_mse = combined_mse,
    combined_aic = combined_aic,
    combined_bic = combined_bic
  )
}




colors <- c('green', 'cornflowerblue', 'goldenrod2', 'lightpink3', 'maroon2', 
            'olivedrab4', 'seagreen', 'navyblue', 'mediumpurple3', 'tomato1', 
            'ivory4', 'lemonchiffon4', 'hotpink4', 'royalblue')

results1 <- fit_legendre_model(data1, genotypes1, degree = 3, color_set = colors, season_label = "2013-14")
results2 <- fit_legendre_model(data2, genotypes2, degree = 3, color_set = colors, season_label = "2014-15")
results3 <- fit_legendre_model(data3, genotypes3, degree = 3, color_set = colors, season_label = "2015-16")
results4 <- fit_legendre_model(data4, genotypes4, degree = 3, color_set = colors, season_label = "2016-17")
results5 <- fit_legendre_model(data5, genotypes5, degree = 3, color_set = colors, season_label = "2017-18")
results6 <- fit_legendre_model(data6, genotypes6, degree = 3, color_set = colors, season_label = "2018-19")
results7 <- fit_legendre_model(data7, genotypes7, degree = 3, color_set = colors, season_label = "2019-20")
results8 <- fit_legendre_model(data8, genotypes8, degree = 3, color_set = colors, season_label = "2020-21")
results9 <- fit_legendre_model(data9, genotypes9, degree = 3, color_set = colors, season_label = "2021-22")
results10 <- fit_legendre_model(data10, genotypes10, degree = 3, color_set = colors, season_label = "2022-23")
results11 <- fit_legendre_model(data11, genotypes11, degree = 3, color_set = colors, season_label = "2023-24")





coefficients_matrices <- list(
  results1$coefficients,
  results2$coefficients,
  results3$coefficients,
  results4$coefficients,
  results5$coefficients,
  results6$coefficients,
  results7$coefficients,
  results8$coefficients,
  results9$coefficients,
  results10$coefficients,
  results11$coefficients
)

genotype_lists <- list(
  genotypes1, genotypes2, genotypes3, genotypes4, genotypes5, 
  genotypes6, genotypes7, genotypes8, genotypes9, genotypes10, genotypes11
)

# Add genotype rownames
for (i in seq_along(coefficients_matrices)) {
  rownames(coefficients_matrices[[i]]) <- genotype_lists[[i]]
}


# Add source labels C1 to C11
source_labels <- paste0("C", 1:11)

# Combine with Source column
combined_coefficients <- do.call(rbind, Map(function(mat, label) {
  cbind(mat, Source = label)
}, coefficients_matrices, source_labels))


head(combined_coefficients)

genotype_names <- rownames(combined_coefficients)
combined_coefficients <- as.data.frame(combined_coefficients, stringsAsFactors = FALSE)
combined_coefficients$Genotype <- genotype_names

#write.csv(combined_coefficients, 'leg_coef.csv')
length(unique(combined_coefficients$Genotype))





# Collect all metric data frames from your results
metrics_list <- list(
  results1$metrics, results2$metrics, results3$metrics, results4$metrics,
  results5$metrics, results6$metrics, results7$metrics, results8$metrics,
  results9$metrics, results10$metrics, results11$metrics
)

# Season labels
source_labels <- paste0("C", 1:11)

# Combine all metrics with Source labels
combined_metrics <- do.call(rbind, Map(function(df, label) {
  df$Source <- label
  df
}, metrics_list, source_labels))

# ----------------------------------------------------------
# Compute ONE mean (average) AIC, BIC, MSE per season
# ----------------------------------------------------------
season_summary <- combined_metrics %>%
  group_by(Source) %>%
  summarise(
    Mean_MSE_D3 = mean(MSE, na.rm = TRUE),
    Mean_AIC_D3 = mean(AIC, na.rm = TRUE),
    Mean_BIC_D3 = mean(BIC, na.rm = TRUE)
  )

# ----------------------------------------------------------
# Save and view results
# ----------------------------------------------------------
write.csv(season_summary, "Metrics_D3.csv", row.names = FALSE)
print(season_summary)

