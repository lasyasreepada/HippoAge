library(starnet)
library(glmnet)
library(dplyr)
library(caret)
library(sesame)
library(doParallel)
library(SummarizedExperiment)
library(glmnetUtils)
library(glmnet)

sesameDataCache()

registerDoParallel(cores = 6)

# Load data (assuming it's in a data frame called 'data')
se <- readRDS('~/Projects/hippoage/data/ADNI/Methylation/20250228_HV_SE.rds')
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
data$Sex <- as.integer(coldata$Sex == 'Female')
data$Leuk <- as.integer(coldata$Leuk >= 0.5)
data$HV <- coldata$M_Hippo_VOL_ASHST1.adjusted

# Split data into features (X) and target (Y)
X <- data[, !colnames(data) %in% c("Sex", "Leuk", "HV")]
Y <- data$HV
covars <- data[,c('Sex', 'Leuk')]

# Remove unnecessary items
rm(data)
rm(se)
rm(betas)
gc()

# Create 80/20 train/test split
set.seed(42)
train_indices <- createDataPartition(covars[,'Sex'], p = 0.8, list = FALSE)
X_train <- X[train_indices,]
Y_train <- Y[train_indices]
covars_train <- covars[train_indices,]

X_test <- X[-train_indices,]
Y_test <- Y[-train_indices]
covars_test <- covars[-train_indices,]

# Function to perform EWAS and select top n CpGs
select_top_cpgs_ewas <- function(X, Y, covars_train, n = 20000) {
  metadf <- as.data.frame(cbind(Y,covars_train))
  X <- t(X)
  ewas <- DML(X, ~Sex+Leuk+Y, meta=metadf) # rows are probes and columns are samples
  smry <- summaryExtractTest(ewas)
  
  # Extract probe IDs for top n CpGs
  helper <- function(smry, feature="Y", top_n) {
    library(dplyr)
    pval = paste0("Pval_",feature)
    cgs <- smry %>% 
      arrange(!!sym(pval)) %>%
      slice_min(order_by = !!sym(pval), n = top_n) %>%
      select(Probe_ID)
    
    cgs <- unlist(cgs)
    return(cgs)
  }
  
  # Return top n cpgs
  top_ewas_cpgs <- helper(smry,top_n=n)
  return(top_ewas_cpgs)
}

top_ewas_cpgs <- select_top_cpgs_ewas(X_train, Y_train, covars_train)

x_num <- scale(X_train[,top_ewas_cpgs])
x <- as.matrix(cbind(x_num,covars_train))
y <- Y_train

# Fit starnet
starnet.hv <- starnet(y=y,X=x,family="gaussian",type.measure='mse',parallel=TRUE,trace.it=1)

# Get the coefficients
# TODO

# Evaluate on train set
y_pred_train_df <- predict(starnet.hv,newx=x,nzero = 500)

# Select top CpGs based on correlation
select_top_model <- function(starnet_df, y) {
  cor_df <- cor(starnet_df,y)
  top_corr_cpgs <- unlist(names(cor_df[max(cor_df)]))
  
  print(mean(abs(cor_df)))
  
  return(top_corr_cpgs)
}

y_pred_train <- y_pred_train_df$
train_mse <- mean((Y_train - y_pred_train)^2)
train_r2 <- cor(y_pred_train, Y_train)**2
print(paste("Train MSE:", train_mse, "Train R2:", train_r2))



# Evaluate on test set
x_test <- as.matrix(cbind(scale(X_test[,top_ewas_cpgs],center = attr(x_num, "scaled:center"), scale = attr(x_num, "scaled:scale")),covars_test))
y_pred_test_df <- predict(starnet.hv, newx = x_test)
y_pred_test <- y_pred_test_df$tune
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
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("train" = "blue", "test" = "red")) +
  labs(title = "Predicted vs Actual Values",
       x = "Actual Values",
       y = "Predicted Values",
       color = "Dataset") +
  theme_minimal()



