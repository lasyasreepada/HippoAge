import pandas as pd
import pyarrow as pq
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_absolute_percentage_error
from sklearn.feature_selection import SelectKBest, f_regression

# Load Data
# data = pd.read_parquet("/Users/lasyasreepada/Projects/hippoage/data/20250228_ADNI_EN.parquet")
data = pd.read_csv('/Users/lasyasreepada/Projects/hippoage/data/ADNI/20250301_HV_RF.csv')
X = data.drop(columns=['HV'])
y = data['HV']
sex = data['Sex']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,random_state=42,stratify=sex,shuffle=True)

# Identify Numeric & Categorical Columns
num_cols = [col for col in X.columns if col != "Sex"]
cat_cols = ['Sex']

# Create the SelectKBest transformer for numerical features
kbest = SelectKBest(score_func=f_regression, k=200)

# Create the ColumnTransformer
preprocessor = ColumnTransformer(
    transformers=[
        ('kbest', kbest, num_cols),
        ('num', StandardScaler(), num_cols),
        ("cat", OneHotEncoder(drop='if_binary'), cat_cols)] # One Hot Encode categorical features
) 

# Define ElasticNet model
elastic_net = ElasticNetCV(
    l1_ratio=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],  # Mix of L1 (Lasso) and L2 (Ridge)
    cv=10,  # Inner CV for alpha selection
    max_iter=10000
)

# Create Pipeline
en_pipeline = Pipeline([
    ("preprocessor", preprocessor),
    ("model", elastic_net) 
])

# Fit pipeline
en_pipeline.fit(X_train,y_train)

# Print results
model = en_pipeline.named_steps['model']
print("Optimal Alpha:", model.alpha_)
print("Optimal L1 ratio:", model.l1_ratio_)
print(f"Number of non-zero features: {np.sum(model.coef_ != 0)}")

plt.figure(figsize=(10, 6))
for i in range(model.mse_path_.shape[0]):
    plt.semilogx(model.alphas_, model.mse_path_[i], label=f'Fold {i+1}')
plt.semilogx(model.alphas_, model.mse_path_.mean(axis=0), 'k--', label='Mean MSE')
plt.xlabel('Alpha')
plt.ylabel('Mean Squared Error')
plt.title('Alpha vs. Mean Squared Error for each fold')
plt.legend()
plt.tight_layout()
plt.show()

y_hat_train = en_pipeline.predict(X_train)
y_hat_test = en_pipeline.predict(X_test)
print(pearsonr(y_hat_train,y_train))
print(pearsonr(y_hat_test,y_test))

# elastic_net.fit(X_train, y_train)
# y_hat_train = inverse_preprocess_age(elastic_net.predict(X_train))
# y_hat_test = inverse_preprocess_age(elastic_net.predict(X_test))

import matplotlib.pyplot as plt

# Scatter Plot for Train Set
plt.figure(figsize=(6, 6))
plt.scatter(y_train, y_hat_train, alpha=0.5, color="red", label="Train")
plt.scatter(y_test, y_hat_test, alpha=0.5, color="blue", label="Test")
plt.plot([y_hat_train.min(), y_hat_train.max()], [y_hat_train.min(), y_hat_train.max()], "k--", lw=2)  # Perfect prediction line
plt.xlabel("Observed Hippocampal Volume - Residual")
plt.ylabel("Predicted Hippocampal Volume - Residual")
plt.title("ElasticNet predicting Hippocampal Volume Residual")
plt.legend()
plt.show()

# Get feature names (after preprocessing)
feature_names = preprocessor.get_feature_names_out()

# Get nonzero coefficients
nonzero_mask = model.coef_ != 0
important_features = feature_names[nonzero_mask]
important_coefs = model.coef_[nonzero_mask]

# Sort by absolute importance
sorted_idx = np.argsort(np.abs(important_coefs))[::-1]  # Descending order
important_features = important_features[sorted_idx]
important_coefs = important_coefs[sorted_idx]

# Plot Top 20 Important Features
plt.figure(figsize=(10, 6))
plt.barh(important_features[:20], important_coefs[:20], color="blue")
plt.xlabel("Coefficient Value")
plt.ylabel("Feature")
plt.title("Top 20 Important Features (Elastic Net)")
plt.gca().invert_yaxis()
plt.show()


