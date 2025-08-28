# Imports
import pandas as pd
import numpy as np
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

# Read data
data = pd.read_csv('/Users/lasyasreepada/Projects/hippoage/data/ADNI/20250301_HV_RF.csv')

X = data.drop(columns=['HV'])
y = data['HV']
sex = data['Leuk']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=sex,random_state=42)

# Create the ColumnTransformer
num_cols = [col for col in X.columns if col not in ["Sex","Leuk"]]
cat_cols = ['Sex','Leuk']

preprocessor = ColumnTransformer(
    transformers=[
        ('num', StandardScaler(), num_cols),
        ], remainder='passthrough' # One Hot Encode categorical features
) 

# Initialize the Extra Trees Regressor
tree = ExtraTreesRegressor(random_state=42, warm_start=True)

# Create Pipeline 
tree_pipe = Pipeline([
    ("preprocessor", preprocessor),
    ("tree", tree)
])

# Set up hyperparameter grid for tuning
# Define the parameter grid
param_grid = {
    'tree__n_estimators': [50, 100],
    'tree__max_depth': [None, 5, 10],
    'tree__min_samples_split': [5, 10],
    'tree__min_samples_leaf': [2, 4],
    'tree__max_samples': [0.7,0.8],
    'tree__bootstrap': [True],
    'tree__criterion': ['squared_error', 'absolute_error', 'friedman_mse']
}

# Perform Grid Search with Cross-Validation
grid_search = GridSearchCV(estimator=tree_pipe, param_grid=param_grid, cv=5, scoring='r2', n_jobs=-1, verbose=1)
grid_search.fit(X_train, y_train)

# Best model from Grid Search
best_tree = grid_search.best_estimator_

# Make predictions on the test set
y_pred_train = best_tree.predict(X_train)
y_pred_test = best_tree.predict(X_test)

# Evaluate the model
mae_train = mean_squared_error(y_train, y_pred_train)
mae_test = mean_squared_error(y_test, y_pred_test)
r2_train = r2_score(y_train, y_pred_train)
r2_test = r2_score(y_test, y_pred_test)

# Output results
print("Best Parameters:", grid_search.best_params_)
print("Mean Absolute Error (MAE) Train:", mae_train)
print("R-squared (R²) Train:", r2_train)
print("Mean Absolute Error (MAE) Test:", mae_test)
print("R-squared (R²) Test:", r2_test)

# Feature importance
feature_names = X.columns
feature_importances = pd.DataFrame({'Feature': feature_names, 'Importance': best_tree['tree'].feature_importances_})
feature_importances = feature_importances.sort_values(by='Importance', ascending=False)
print("\nTop 10 Feature Importances:")
print(feature_importances.head(10))

import matplotlib.pyplot as plt

# Scatter Plot for Train Set
plt.figure(figsize=(6, 6))
plt.scatter(y_train, y_pred_train, alpha=0.5, color="red", label="Train")
plt.scatter(y_test, y_pred_test, alpha=0.5, color="blue", label="Test")
# plt.plot([y_pred_train.min(), y_pred_train.max()], [y_pred_train.min(), y_pred_train.max()], "k--", lw=2)  # Perfect prediction line
plt.xlabel("Observed Hippocampal Volume")
plt.ylabel("Predicted Hippocampal Volume")
plt.title("Random Forest predicting Hippocampal Volume")
plt.legend()
plt.show()
