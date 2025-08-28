library(dplyr)
library(readxl)
library(glmnet)
library(dplyr)
library(caret)
library(sesame)
library(doParallel)
library(SummarizedExperiment)
library(glmnetUtils)
library(glmnet)
library(starnet)

registerDoParallel(cores = 6)

# Load data (assuming it's in a data frame called 'data')
se <- readRDS('~/Projects/hippoage/data/ADNI/20250314_COMET_SE.rds')
se <- se[,se$Case==1]

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
data$CoMeT <- coldata$CorticalLOAD_CoMeT
sex <- as.integer(coldata$Sex == 'Female')
covars <- coldata[,c('PlateNumber')]

# Split data into features (X) and target (Y)
X <- data[, !colnames(data) %in% c("CoMeT")]
Y <- data$CoMeT

# Create 80/20 train/test split
set.seed(42)
train_indices <- createDataPartition(sex, p = 0.8, list = FALSE)
X_train <- X[train_indices,]
Y_train <- Y[train_indices]
covars_train <- covars[train_indices,]

X_test <- X[-train_indices,]
Y_test <- Y[-train_indices]

# Select top CpGs based on correlation
select_top_cpgs_corr <- function(X, Y, corr) {
  cor_df <- cor(X,Y)
  top_corr_cpgs <- unlist(names(cor_df[abs(cor_df) >= corr,]))
  
  print(mean(abs(cor_df)))
  
  return(top_corr_cpgs)
}

# Function to perform EWAS and select top n CpGs
select_top_cpgs_var <- function(X, threshold = 0.9) {
  
  # Compute variance for each CpG site
  cpg_variances <- apply(X, 2, var)
  
  # Filter CpGs with variance above a threshold
  filter <- quantile(cpg_variances, threshold)  # Keep top 10% most variable
  high_var_cpgs <- names(cpg_variances[cpg_variances >= filter])
  
  # Return
  return(high_var_cpgs)
}

# Function to perform EWAS and select top n CpGs
select_top_cpgs <- function(X, Y, covars, n = 1000) {
  metadf <- as.data.frame(cbind(Y,covars))
  X <- t(X)
  ewas <- DML(X, ~ PlateNumber+Y, meta=metadf) # rows are probes and columns are samples
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

# Select top CpGs from univariate analysis
top_cpgs <- unlist(select_top_cpgs(X_train,Y_train,covars_train,1000))
top_n <- top_cpgs[1:375]

# Apply boruta
# library(Boruta)
# boruta <- Boruta(X_train,Y_train,pValue = 0.1,mcAdj = TRUE,doTrace = 2)
# features <- names(boruta$finalDecision[boruta$finalDecision == "Confirmed"])

x_num <- X_train[,top_n]
x <- as.matrix(x_num)
y <- Y_train

# Perform cross-validation
en <- cv.glmnet(x,y,type.measure='mse',family='gaussian',alpha=0.3,parallel=TRUE,trace.it=1)

# Assess cross validation CV metrics
mean(en$cvm)
var(y)

# Get the best parameters
best_lambda.1se <- en$lambda.1se
best_lambda.min <- en$lambda.min

# Evaluate on train set
y_pred_train <- predict(en, newx = x, s=best_lambda.min, type='response')
train_mse <- mean((Y_train - y_pred_train)^2)
train_r2 <- cor(y_pred_train, Y_train)**2
print(paste("Train MSE:", train_mse, "Train R2:", train_r2))

# Evaluate on test set
cgs <- rownames(coef(en))[-1]
x_test_num <- X_test[,cgs]
x_test <- as.matrix(x_test_num)
y_pred_test <- predict(en, newx = x_test,s=best_lambda.min,type='response')
test_mse <- mean((Y_test - y_pred_test)^2)
test_r2 <- cor(y_pred_test, Y_test)**2
print(paste("Test MSE:", test_mse, "Test R2:", test_r2))

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
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("train" = "blue", "test" = "red")) +
  labs(title = "Predicted vs Actual Values",
       x = "Actual Values",
       y = "Predicted Values",
       color = "Dataset") +
  theme_minimal()

# Feature importances
coefficients <- coef(en, s = "lambda.min")
importances <- as.matrix(coefficients)[-1, ] # Remove intercept
importances <- sort(importances, decreasing = TRUE)
imp_df <- data.frame(name = names(importances), value = importances)
rownames(imp_df) <- NULL
imp_df <- imp_df[abs(imp_df$value) > 0,]

write.csv(imp_df,'~/Projects/hippoage/data/ADNI/20250305_HV_EN_Features.csv',row.names = FALSE)
