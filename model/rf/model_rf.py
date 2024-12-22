import numpy as np
import pandas as pd
from sksurv.ensemble import RandomSurvivalForest
from sksurv.preprocessing import OneHotEncoder
from sksurv.util import Surv
import sys


args = sys.argv[1:]
project_name = args[0]

# Read in & Process Data
df = pd.read_csv('../../pathway_map/'+project_name+'/pathway_counts.csv')
df.dropna(inplace=True, how='any')

X = df.iloc[:,3:] # select columns with features (currently gene info)
surv_info = df[['status','months']]  # columns with survival info
y = np.array(list(zip(surv_info.status, surv_info.months)) , dtype =[('Status', '?'), ('Survival_in_months', '<f8')])
Xt = OneHotEncoder().fit_transform(X) # one hot encoder w/ categorical variables

# Initialize and fit the Random Survival Forest model
rsf = RandomSurvivalForest()

from sklearn.model_selection import RandomizedSearchCV, KFold
from sksurv.metrics import concordance_index_censored
from scipy.stats import randint

param_distributions = {
    'n_estimators': randint(50, 100),
    'min_samples_split': randint(2, 40),
    'min_samples_leaf': randint(1, 20),
    'max_depth': randint(5, 100),
    'max_features': ['sqrt', 'log2', None]
}
# Define the scoring function
def neg_cindex_scorer(estimator, X, y):
    y_pred = estimator.predict(X)
    return -concordance_index_censored(y['Status'], y['Survival_in_months'], y_pred)[0]



cv = KFold(n_splits=5, shuffle=True, random_state=42)
# Define the RandomizedSearchCV
random_search = RandomizedSearchCV(estimator=rsf, param_distributions=param_distributions, n_iter=100, cv=cv, n_jobs=-1, scoring=neg_cindex_scorer, random_state=42, verbose=1)

# Fit the GridSearchCV
random_search.fit(X, y)

# Get the best parameters
best_params = random_search.best_params_
print("Best parameters found: ", best_params)

# Get the best estimator
best_model = random_search.best_estimator_
print("Best model: ", best_model)

from sklearn.model_selection import train_test_split
Xtrain, Xtest, ytrain, ytest = train_test_split(Xt, y, test_size=0.1, shuffle = True, random_state=1)


Xtrain, Xtest, ytrain, ytest = train_test_split(Xt, y, test_size=0.1, shuffle = True, random_state=1)
best_model.fit(Xtrain, ytrain)

surv_fns = best_model.predict_survival_function(Xtrain)
print("Concordance: {:.2f}".format(best_model.score(Xtrain, ytrain)))
print("Concordance: {:.2f}".format(best_model.score(Xtest, ytest)))

import joblib
best_model.fit(X, y)
print("Concordance: {:.2f}".format(best_model.score(X, y)))
try:
    while True:
        save = input('Save?(y/n)')
        if save == 'y':
            joblib_file = project_name + "_rf.pkl" 
            joblib.dump(best_model, joblib_file)
            print('Model saved')
            break
        elif save == 'n':
            print('Model not saved')
            break
except:
    joblib_file = project_name + "_rf_auto_backup.pkl" 
    joblib.dump(best_model, joblib_file)