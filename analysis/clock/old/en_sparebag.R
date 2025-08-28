library(glmnet)
library(dplyr)
library(caret)
library(sesame)
library(doParallel)
library(SummarizedExperiment)
library(glmnetUtils)
library(glmnet)

# registerDoParallel(cores = 6)
sesameDataCache()

# Load data (assuming it's in a data frame called 'data')
se <- readRDS('~/Projects/hippoage/data/ADNI/Methylation/20250227_SPAREBA_SE.rds')
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
data$Age <- coldata$Age_MRI  
data$Sex <- coldata$Sex
data$Leuk <- coldata$Leuk
data$SPAREBAG <- coldata$SPARE_BA.SVM.RBF. - coldata$Age_MRI

# Split data into features (X) and target (Y)
X <- data[, !colnames(data) %in% c("Age", "Sex", "Leuk", "SPAREBAG")]
Y <- data$SPAREBAG
covars <- data[,c('Sex', 'Leuk')]

# Remove unnecessary items
rm(data)
rm(se)
rm(betas)
gc()

# Create 80/20 train/test split
set.seed(42)
train_indices <- createDataPartition(covars[, 'Sex'], p = 0.8, list = FALSE)
X_train <- X[train_indices,]
Y_train <- Y[train_indices]
covars_train <- covars[train_indices,]

X_test <- X[-train_indices,]
Y_test <- Y[-train_indices]

# Function to perform EWAS and select top n CpGs
select_top_cpgs_var <- function(X, threshold) {
  
  # Compute variance for each CpG site
  cpg_variances <- apply(X, 2, var)
  
  # Filter CpGs with variance above a threshold
  filter <- quantile(cpg_variances, threshold)  # Keep top 10% most variable
  high_var_cpgs <- names(cpg_variances[cpg_variances >= filter])
  
  # Return
  return(high_var_cpgs)
}

# Function to perform EWAS and select top n CpGs
select_top_cpgs_ewas <- function(X, Y, covars_train, n = 20000) {
  metadf <- as.data.frame(cbind(Y,covars_train))
  X <- t(X)
  ewas <- DML(X, ~Age+Sex+Leuk+Y, meta=metadf) # rows are probes and columns are samples
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
  return(top_ewas_cpgs)
}

# Select top CpGs based on correlation
select_top_cpgs_corr <- function(X, Y, corr) {
  cor_df <- cor(X,Y)
  top_corr_cpgs <- unlist(names(cor_df[abs(cor_df) >= corr,]))
  
  print(mean(abs(cor_df)))
  
  return(top_corr_cpgs)
}

# Perform EWAS and select top n CpGs
# top_cpgs <- unlist(select_top_cpgs(X_train, Y_train, covars_train))

# Select CpGs with most variance
top_cpgs_var <- select_top_cpgs_var(X_train,0.9)

# Select CpGs with highest correlation with SPAREBA
top_cpgs_corr <- select_top_cpgs_corr(X_train,Y_train,0.2)

# Prepare data for glmnet
top_cpgs <- intersect(top_cpgs_var,top_cpgs_corr)

x <- as.matrix(X_train[,top_cpgs])
y <- Y_train

# Perform cross-validation
cv_model <- cv.glmnet(x, y, type.measure = "mse", alpha = 0.2, nfolds = 10, parallel = TRUE, trace.it = 1)

# Assess cross validation CV metrics
mean(cv_model$cvm)
var(y)

# Get the best parameters
best_lambda.1se <- cv_model$lambda.1se
best_lambda.min <- cv_model$lambda.min

# Evaluate on train set
y_pred_train <- predict(cv_model, newx = x, s = best_lambda.min)
train_mse <- mean((Y_train - y_pred_train)^2)
train_cor <- cor(y_pred_train, Y_train)
print(paste("Train MSE:", train_mse, "Train Pearson R:", train_cor))

# Evaluate on test set
x_test <- as.matrix(X_test[,top_cpgs])
y_pred_test <- predict(cv_model, newx = x_test, s = best_lambda.min)
test_mse <- mean((Y_test - y_pred_test)^2)
test_cor <- cor(y_pred_test, Y_test)
print(paste("Test MSE:", test_mse, "Test Pearson R:", test_cor))

train_data <- data.frame(
  y_pred = as.vector(y_pred_train),
  y = Y_train,
  type = "train"
)

test_data <- data.frame(
  y_pred = as.vector(y_pred_test),
  y = Y_test,
  type = "test"
)

# Combine train and test data
combined_data <- bind_rows(train_data, test_data)

# Create the plot
ggplot(combined_data, aes(x = y, y = y_pred, color = type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("train" = "blue", "test" = "red")) +
  labs(title = "Predicted vs Actual Values",
       x = "Predicted Values",
       y = "Actual Values",
       color = "Dataset") +
  theme_minimal()
