import pandas as pd
from sklearn.model_selection import GroupKFold, cross_val_score
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error

# Assuming 'df' is your DataFrame containing CpG values, hippocampal volume, and Sex
# Columns: CpG1, CpG2, ..., CpGn, hippocampal_volume, sex
df = pd.read_csv('~/Projects/hippoage/data/ADNI/20250218_EN_Data.csv')

# Load the data (assuming df has columns 'CpG1', 'CpG2', ..., 'hippocampal_volume', 'sex')
X = df.drop(columns=['Unnamed: 0','HV_res', 'Sex'])  # CpG values as features
y = df['HV_res']  # Target variable: hippocampal volume
sex = df['Sex']  # Used for grouping

# Initialize ElasticNetCV model
elastic_net_cv = ElasticNetCV(cv=5,  # Default 5-fold cross-validation
                             l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9, 1.0],  # Range of l1_ratio values to try
                             alphas=[0.001, 0.01, 0.1, 1, 10, 100],  # Range of alpha values to try
                             max_iter=10000)  # Set maximum iterations if needed

# Create a pipeline with StandardScaler and ElasticNetCV
pipeline = Pipeline([
    ('scaler', StandardScaler()),  # Standardize the CpG values
    ('elasticnet', elastic_net_cv)  # ElasticNetCV for automatic hyperparameter tuning
])

# Set up GroupKFold cross-validation, stratifying by 'sex'
gkf = GroupKFold(n_splits=2)  # Only 2 splits since we have 2 groups (Male, Female)
elastic_net.fit(X,y)
# Perform cross-validation and calculate mean squared error for each fold
mse_scores = cross_val_score(pipeline, X, y, cv=gkf.split(X, y, groups=sex), scoring='mean_absolute_error')

# Convert the negative MSE to positive
mse_scores = -mse_scores

# Print the cross-validation results
print(f"Mean squared error for each fold: {mse_scores}")
print(f"Average MSE: {mse_scores.mean()}")

# Access the best parameters
best_alpha = elastic_net_cv.alpha_
best_l1_ratio = elastic_net_cv.l1_ratio_

print(f"Best alpha: {best_alpha}")
print(f"Best l1_ratio: {best_l1_ratio}")