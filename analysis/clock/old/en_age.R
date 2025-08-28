library(glmnet)
library(caret)
library(sesame)
library(SummarizedExperiment)

sesameDataCache()

# Load data (assuming it's in a data frame called 'data')
se <- readRDS('~/Projects/hippoage/data/ADNI/EWAS/20250224_HV_SE.rds')
betas <- assays(se)[[1]]
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Impute missing CpGs
imputeRowMean <- function(mtx) {
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
  return(mtx)
}

betas <- imputeRowMean(betas)
data <- as.data.frame(t(betas))

data$Sex <- coldata$Sex
# data$Hippocampal_Volume_Adjusted <- coldata$M_Hippo_VOL_ASHST1.adjusted
data$Age <- coldata$Age_MRI  

preprocess_age <- function(age) {
  # Preprocesses the age column using the formula F(age) = (age - 20) / (20 + 1).
  return((age - 20) / 21)
}

inverse_preprocess_age <- function(f_age) {
  # Computes the inverse transformation to recover the original age.
  return((f_age * 21) + 20)
}
  
# Split data into features (X) and target (Y)
# X <- data[, !colnames(data) %in% c("Hippocampal_Volume_Adjusted", "Sex")]
X <- data[, !colnames(data) %in% c("Age", "Sex")]
Y <- preprocess_age(data$Age)
sex <- data$Sex

# Create 80/20 train/test split
set.seed(42)
train_indices <- createDataPartition(Y, p = 0.8, list = FALSE)
X_train <- X[train_indices, ]
Y_train <- Y[train_indices]

X_test <- X[-train_indices, ]
Y_test <- Y[-train_indices]
sex_train <- sex[train_indices]

# Create k stratified folds
k_fold <- 2
folds <- createFolds(sex_train, k = k_fold, list = TRUE, returnTrain = FALSE)

# Function to perform EWAS and select top n CpGs
select_top_cpgs <- function(X, Y, sex, n = 4000) {
  metadf <- as.data.frame(cbind(Y,sex))
  X <- t(X)
  ewas <- DML(X, ~sex+Y, meta=metadf) # rows are probes and columns are samples
  smry <- summaryExtractTest(ewas)
  
  # Extract probe IDs for top n CpGs
  helper <- function(smry, feature="Y", top_n) {
    library(dplyr)
    pval = paste0("Pval_",feature)
    cgs <- smry %>% 
      arrange(!!sym(pval)) %>%
      slice_min(order_by = !!sym(pval), n = top_n) %>%
      select(Probe_ID)
    return(cgs)
  }
  
  # Return top n cpgs
  top_cpgs <- helper(smry,top_n=n)
  return(top_cpgs)
}

registerDoParallel(cores = 6)

# Perform cross-validation
cv_results <- lapply(1:k_fold, function(i) {
  train <- unlist(folds[-i])
  val <- unlist(folds[i])
  
  # Perform EWAS and select top n CpGs
  top_cpgs <- unlist(select_top_cpgs(X_train[train, ], Y_train[train], sex_train[train]))
  
  # Prepare data for glmnet
  x <- as.matrix(scale(X_train[train, top_cpgs]))
  y <- Y_train[train]
  
  # Perform cross-validation
  cv_fit <- cv.glmnet(x, y, alpha = 0.5, nfolds = 5, parallel = TRUE, trace.it = 1)
  
  # Evaluate on validation set
  x_val <- as.matrix(scale(X_train[val, top_cpgs], center = attr(x, "scaled:center"), 
                           scale = attr(x, "scaled:scale")))
  y_val <- Y_train[val]
  pred <- predict(cv_fit, newx = x_val, s = "lambda.min")
  mse <- mean((y_val - pred)^2)
  
  list(top_cpgs = top_cpgs, cv_fit = cv_fit, mse = mse)
})

# Evaluate CV performance
print(paste("Mean MSE across folds:", mean(unlist(lapply(cv_results, function(x) x$mse)))))
print(paste("Target Variance:",var(Y_train)))
      
# Find union CpGs across all folds
overlap_cpgs <- Reduce(intersect, lapply(cv_results, function(x) x$top_cpgs))
# union_cpgs <- Reduce(union, lapply(cv_results, function(x) x$top_cpgs))

# Train final model using union CpGs and best parameters
x_final <- as.matrix(X_train[, overlap_cpgs])
y_final <- Y_train

best_lambda <- mean(sapply(cv_results, function(x) x$cv_fit$lambda.min))
final_model <- glmnet(x_final, y_final, alpha = 0.5, lambda = best_lambda)

print(paste("Number of CpGs in final model:", length(coef(final_model)@i)))

# Evaluate on train set
x_train_final <- as.matrix(X_train[, overlap_cpgs])
y_pred_train <- predict(final_model, newx = x_train_final)
train_mse <- mean((Y_train - y_pred_train)^2)
train_cor <- cor(y_pred_train, Y_train)
print(paste("Train MSE:", train_mse, "Train Pearson R:", train_cor))

# Evaluate on test set
x_test_final <- as.matrix(X_test[, overlap_cpgs])
y_pred <- predict(final_model, newx = x_test_final)
test_mse <- mean((Y_test - y_pred)^2)
test_cor <- cor(y_pred, Y_test)
print(paste("Test MSE:", test_mse, "Test Pearson R:", test_cor))


