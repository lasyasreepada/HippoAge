# Imports
import pandas as pd
import numpy as np
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

# Read data
data = pd.read_csv('/Users/lasyasreepada/Projects/hippoage/data/ADNI/20250301_HV_RF.csv')

X = data.drop(columns=['HV'])
y = data['HV']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

ols = LinearRegression()
ols.fit(X_train, y_train)

y_pred_train = ols.predict(X_train)
y_pred = ols.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(ols.coef_)      # Coefficients
print(ols.intercept_) # Intercept

print("Mean Squared Error (MAE) Test:", mse)
print("R-squared (R²) Test:", r2)

import matplotlib.pyplot as plt

# Scatter Plot for Train Set
plt.figure(figsize=(6, 6))
plt.scatter(y_train, y_pred_train, alpha=0.5, color="red", label="Train")
plt.scatter(y_test, y_pred, alpha=0.5, color="blue", label="Test")
# plt.plot([y_pred_train.min(), y_pred_train.max()], [y_pred_train.min(), y_pred_train.max()], "k--", lw=2)  # Perfect prediction line
plt.xlabel("Observed Hippocampal Volume")
plt.ylabel("Predicted Hippocampal Volume")
plt.title("Linear Regression Predicting Hippocampal Volume")
plt.legend()
plt.show()