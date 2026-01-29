rm(list=ls())

library(dplyr)

#setwd('')

list.files()

data <- read.csv('leg_coef.csv')

data$CV <- as.numeric(sub("C", "", data$Source))

data$CV_11 <- ifelse(data$CV == 11, NA, data$CV)
data$CV_10 <- ifelse(data$CV == 10, NA, data$CV_11)
data$CV_9 <- ifelse(data$CV == 9, NA, data$CV_10)

head(data)
Y <- data

write.csv(Y, 'Y.csv')
load('G.rda')
load('E.rda')
load('GxE.rda')


###Model 2 


library(BGLR)

# -----------------------------
# Settings
# -----------------------------
nIter <- 5000
burnIn <- 1000

# Response matrix (4 traits)
Y2 <- as.matrix(Y[, 2:5])
rownames(Y2) <- Y$Genotype

# ETA structure
ETA <- list(
  list(K = G, model = 'RKHS'),
  list(K = E, model = 'RKHS'),
  list(K= GxE, model='RKHS')
)

####Function

run_bglr_prediction <- function(Y, Y2, cv_column, ETA, nIter, burnIn) {
  
  # Initialize empty list to store predictions for each trait
  TF <- vector("list", ncol(Y2))
  for (i in seq_along(TF)) {
    TF[[i]] <- rep(NA, nrow(Y))  # initialize with NAs
  }
  
  # Loop through traits
  for (trait in 1:ncol(Y2)) {
    yNA <- Y2[, trait]
    
    # Mask out values where that CV column is NA
    yNA[is.na(Y[[cv_column]])] <- NA
    
    # Run BGLR
    fm <- BGLR(y = yNA, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = TRUE)
    
    # Store predicted values for NA positions only
    TF[[trait]][is.na(Y[[cv_column]])] <- fm$yHat[is.na(Y[[cv_column]])]
  }
  
  # Combine predictions
  TF_mat <- do.call(cbind, TF)
  colnames(TF_mat) <- paste0("Trait_", 1:ncol(Y2))
  rownames(TF_mat) <- Y$X
  
  return(TF_mat)
}




# Predict for environment CV = 11
pred_CV11 <- run_bglr_prediction(Y, Y2, "CV_11", ETA, nIter, burnIn)
write.csv(pred_CV11, "Predictions_CV11_M2.csv")


Y_true_11 <- Y[is.na(Y$CV_11), 2:5]  # observed Y values for environment 11
P_pred_11 <- pred_CV11[is.na(Y$CV_11), ]

cor_results_11 <- diag(cor(Y_true_11, P_pred_11))
names(cor_results_11) <- paste0("Trait_", 1:4)
print(cor_results_11)



# Predict for environment CV = 10
pred_CV10 <- run_bglr_prediction(Y, Y2, "CV_10", ETA, nIter, burnIn)
write.csv(pred_CV10, "Predictions_CV10_M2.csv")

Y_true_10 <- Y[Y$CV == 10, 2:5]
P_pred_10 <- pred_CV10[Y$CV == 10, ]
cor_results_10 <- diag(cor(Y_true_10, P_pred_10))

names(cor_results_10) <- paste0("Trait_", 1:4)
print(cor_results_10)


pred_CV9 <- run_bglr_prediction(Y, Y2, "CV_9", ETA, nIter, burnIn)
write.csv(pred_CV9, "Predictions_CV9_M2.csv")

Y_true_9  <- Y[Y$CV == 9,  2:5]
P_pred_9  <- pred_CV9[Y$CV == 9, ]
cor_results_9  <- diag(cor(Y_true_9, P_pred_9))

names(cor_results_9) <- paste0("Trait_", 1:4)
print(cor_results_9)





###Model 2 


library(BGLR)

# -----------------------------
# Settings
# -----------------------------
nIter <- 5000
burnIn <- 1000

# Response matrix (4 traits)
Y2 <- as.matrix(Y[, 2:5])
rownames(Y2) <- Y$Genotype

# ETA structure
ETA <- list(
  list(K = G, model = 'RKHS'),
  list(K = E, model = 'RKHS'),
  list(K= GxE, model='RKHS')
)


####Function

run_bglr_prediction <- function(Y, Y2, cv_column, ETA, nIter, burnIn) {
  
  # Initialize empty list to store predictions for each trait
  TF <- vector("list", ncol(Y2))
  for (i in seq_along(TF)) {
    TF[[i]] <- rep(NA, nrow(Y))  # initialize with NAs
  }
  
  # Loop through traits
  for (trait in 1:ncol(Y2)) {
    yNA <- Y2[, trait]
    
    # Mask out values where that CV column is NA
    yNA[is.na(Y[[cv_column]])] <- NA
    
    # Run BGLR
    fm <- BGLR(y = yNA, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = TRUE)
    
    # Store predicted values for NA positions only
    TF[[trait]][is.na(Y[[cv_column]])] <- fm$yHat[is.na(Y[[cv_column]])]
  }
  
  # Combine predictions
  TF_mat <- do.call(cbind, TF)
  colnames(TF_mat) <- paste0("Trait_", 1:ncol(Y2))
  rownames(TF_mat) <- Y$X
  
  return(TF_mat)
}




# Predict for environment CV = 11
pred_CV11 <- run_bglr_prediction(Y, Y2, "CV_11", ETA, nIter, burnIn)
write.csv(pred_CV11, "Predictions_CV11_M2.csv")


Y_true_11 <- Y[is.na(Y$CV_11), 2:5]  # observed Y values for environment 11
P_pred_11 <- pred_CV11[is.na(Y$CV_11), ]

cor_results_11 <- diag(cor(Y_true_11, P_pred_11))
names(cor_results_11) <- paste0("Trait_", 1:4)
print(cor_results_11)



# Predict for environment CV = 10
pred_CV10 <- run_bglr_prediction(Y, Y2, "CV_10", ETA, nIter, burnIn)
write.csv(pred_CV10, "Predictions_CV10_M2.csv")

Y_true_10 <- Y[Y$CV == 10, 2:5]
P_pred_10 <- pred_CV10[Y$CV == 10, ]
cor_results_10 <- diag(cor(Y_true_10, P_pred_10))

names(cor_results_10) <- paste0("Trait_", 1:4)
print(cor_results_10)


pred_CV9 <- run_bglr_prediction(Y, Y2, "CV_9", ETA, nIter, burnIn)
write.csv(pred_CV9, "Predictions_CV9_M2.csv")

Y_true_9  <- Y[Y$CV == 9,  2:5]
P_pred_9  <- pred_CV9[Y$CV == 9, ]
cor_results_9  <- diag(cor(Y_true_9, P_pred_9))

names(cor_results_9) <- paste0("Trait_", 1:4)
print(cor_results_9)



# Combine all results into a data frame
cor_summary <- data.frame(
  CV = c(9, 10, 11),
  Trait_1 = c(cor_results_9[1],  cor_results_10[1],  cor_results_11[1]),
  Trait_2 = c(cor_results_9[2],  cor_results_10[2],  cor_results_11[2]),
  Trait_3 = c(cor_results_9[3],  cor_results_10[3],  cor_results_11[3]),
  Trait_4 = c(cor_results_9[4],  cor_results_10[4],  cor_results_11[4])
)

print(cor_summary)

write.csv(cor_summary, "Correlation_coef_M2.csv", row.names = FALSE)


####Model 1


# ETA structure for M1 (no GxE term)
ETA_M1 <- list(
  list(K = G, model = 'RKHS'),
  list(K = E, model = 'RKHS')
)

# M1 predictions for CV = 11
pred_CV11_M1 <- run_bglr_prediction(Y, Y2, "CV_11", ETA_M1, nIter, burnIn)
write.csv(pred_CV11_M1, "Predictions_CV11_M1.csv")

# M1 predictions for CV = 10
pred_CV10_M1 <- run_bglr_prediction(Y, Y2, "CV_10", ETA_M1, nIter, burnIn)
write.csv(pred_CV10_M1, "Predictions_CV10_M1.csv")

# M1 predictions for CV = 9
pred_CV9_M1 <- run_bglr_prediction(Y, Y2, "CV_9", ETA_M1, nIter, burnIn)
write.csv(pred_CV9_M1, "Predictions_CV9_M1.csv")


# --- CV = 11 ---
Y_true_11 <- Y[Y$CV == 11, 2:5]
P_pred_11 <- pred_CV11_M1[Y$CV == 11, ]
cor_results_11 <- diag(cor(Y_true_11, P_pred_11))
names(cor_results_11) <- paste0("Trait_", 1:4)

# --- CV = 10 ---
Y_true_10 <- Y[Y$CV == 10, 2:5]
P_pred_10 <- pred_CV10_M1[Y$CV == 10, ]
cor_results_10 <- diag(cor(Y_true_10, P_pred_10))
names(cor_results_10) <- paste0("Trait_", 1:4)

# --- CV = 9 ---
Y_true_9 <- Y[Y$CV == 9, 2:5]
P_pred_9 <- pred_CV9_M1[Y$CV == 9, ]
cor_results_9 <- diag(cor(Y_true_9, P_pred_9))
names(cor_results_9) <- paste0("Trait_", 1:4)


# Combine all correlation results for M1
cor_summary_M1 <- data.frame(
  CV = c(9, 10, 11),
  Trait_1 = c(cor_results_9[1],  cor_results_10[1],  cor_results_11[1]),
  Trait_2 = c(cor_results_9[2],  cor_results_10[2],  cor_results_11[2]),
  Trait_3 = c(cor_results_9[3],  cor_results_10[3],  cor_results_11[3]),
  Trait_4 = c(cor_results_9[4],  cor_results_10[4],  cor_results_11[4])
)

# Save to CSV
write.csv(cor_summary_M1, "Correlation_coef_M1.csv", row.names = FALSE)
print(cor_summary_M1)

