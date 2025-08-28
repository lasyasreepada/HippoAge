library(glmnet)
library(dplyr)
library(caret)
library(sesame)
library(doParallel)
library(SummarizedExperiment)
library(glmnetUtils)
library(glmnet)
library(doParallel)
library(maxprobes)
library(readxl)

sesameDataCache()

# Load data (assuming it's in a data frame called 'data')
se <- readRDS('~/Projects/hippoage/data/ADNI/20250314_ADNI_SE.rds')

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

coldata$SPAREBAR <- coldata$SPARE_BA_slope / mean(coldata$SPARE_BA_slope)
data$SPAREBAR <- coldata$SPAREBAR

# Split data into features (X) and target (Y)
X <- data[, !colnames(data) %in% c("SPAREBAR")]
Y <- data$SPAREBAR
covars <- coldata[,c('Sex')]

# Remove unnecessary items
rm(se)
rm(data)
gc()

# Create 80/20 train/test split
set.seed(42)
train_indices <- createDataPartition(covars, p = 0.8, list = FALSE)

# Train
X_train <- X[train_indices,]
Y_train <- Y[train_indices]
covars_train <- covars[train_indices]

# Test
X_test <- X[-train_indices,]
Y_test <- Y[-train_indices]
covars_test <- covars[-train_indices]

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

# Select top CpGs based on correlation
select_top_cpgs_corr <- function(X, Y, corr) {
  cor_df <- cor(X,Y)
  top_corr_cpgs <- unlist(names(cor_df[abs(cor_df) >= corr,]))
  
  print(mean(abs(cor_df)))
  
  return(top_corr_cpgs)
}

# Include only reliable cpgs
sugden <- read_excel('~/Downloads/mmc2.xlsx',col_names = TRUE, skip = 2)
sugden_cpgs <- unlist(sugden[sugden$Reliability>0.2,'Illumina Probe ID'])

sugden_cpgs <- intersect(sugden_cpgs,rownames(betas))

# Select CpGs with most variance
# top_cpgs_var <- select_top_cpgs_var(X_train,0.8)

# Select CpGs with highest correlation with SPAREBA
top_cpgs_corr <- select_top_cpgs_corr(X_train,Y_train,0.3)

# Final subset
# top_cpgs <- intersect(top_cpgs_var,top_cpgs_corr)
# top_cpgs <- top_cpgs_var
# top_cpgs <- top_cpgs_corr
top_cpgs <- sugden_cpgs

# x_num <- X_train
x_num <- X_train[,top_cpgs]
x <- as.matrix(x_num)
y <- Y_train

library(caret)
# Define the train control with cross-validation
train_control <- trainControl(method = "cv", 
                              number = 5, # Number of folds
                              savePredictions = "final")

# Train the model with cross-validation
model <- train(x = x, 
               y = y, 
               method = "glmnet", 
               tuneGrid = expand.grid(alpha = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                                      lambda = exp(seq(log(0.0001), log(100), length.out = 100))),
               trControl = train_control)

final_model <- glmnet(x,y,alpha=0.9,lambda=0.1384886)

# Evaluate on train set
y_pred_train <- predict(final_model, newx = x)
train_mse <- mean((Y_train - y_pred_train)^2)
train_cor <- cor(y_pred_train, Y_train)
print(paste("Train MSE:", train_mse, "Train Pearson R:", train_cor))

# Evaluate on test set
x_num_test <- X_test[,top_cpgs]
x_test <- as.matrix(x_num_test)
y_pred_test <- predict(final_model, newx = x_test)
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
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("train" = "blue", "test" = "red")) +
  labs(title = "Predicted vs Actual Values",
       x = "Actual Values",
       y = "Predicted Values",
       color = "Dataset") +
  theme_minimal()

# Train in full dataset 
# cv_model_full <- cv.glmnet(as.matrix(X[,top_cpgs]), Y, type.measure = "mse", alpha = 0.1, nfolds = 10, parallel = TRUE, trace.it = 1)


