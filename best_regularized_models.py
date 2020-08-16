#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
Obtain the best regularized parameters
for common models using GridSearch and CV.
Present models: 
   -- Ridge 
   -- Lasso
   -- Elastic Net
   -- Kernel Ridge
   -- Gaussian process Regressor
   -- RandomForest Regressor
   -- Gradient Boosting Regressor
"""

import numpy as np
from sklearn.linear_model import LinearRegression, Lasso, Ridge, ElasticNet
from sklearn.kernel_ridge import KernelRidge
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.model_selection import GridSearchCV


class regularized_models():


   def __init__(self, X, y, cv=5, score='neg_mean_absolute_error', n_jobs=None): 
       self.X = X
       self.y = y
       self.cv = cv
       self.score = score
       self.n_jobs = n_jobs


   def best_Ridge(self, alphas):
  
       model = GridSearchCV(Ridge(), param_grid={"alpha": alphas}, 
                   cv=self.cv, scoring=self.score, n_jobs=self.n_jobs) 
       model = model.fit(self.X, self.y)
       return model.best_estimator_, model.best_score_


   def best_Lasso(self, alphas):
  
       model = GridSearchCV(Lasso(), param_grid={"alpha": alphas}, 
                   cv=self.cv, scoring=self.score, n_jobs=self.n_jobs) 
       model = model.fit(self.X, self.y)
       return model.best_estimator_, model.best_score_
    

   def best_ElasticNet(self, alphas):
 
#       l1_ratios = [.1, .5, .7, .9, .95, .99, 1]
       l1_ratios = [.1, .5, .9]
       model = GridSearchCV(ElasticNet(), param_grid={"alpha": alphas,
                                  	       "l1_ratio": l1_ratios},
                   cv=self.cv, scoring=self.score, n_jobs=self.n_jobs) 
       model = model.fit(self.X, self.y)
       return model.best_estimator_, model.best_score_


   def best_KernelRidge(self, alphas):
       
       gammas = np.logspace(-12, -9, 5) 
       model = GridSearchCV(KernelRidge(), param_grid={"alpha": alphas,
                                     		       "gamma": gammas, 
				      "kernel" : ['laplacian', 'rbf']},
                    cv=self.cv, scoring=self.score, n_jobs=self.n_jobs) 
       model = model.fit(self.X, self.y)
       return model.best_estimator_, model.best_score_


   def best_GPR(self, alphas):
 
       model = GridSearchCV(GaussianProcessRegressor(normalize_y=True),
                                          param_grid={"alpha": alphas},
                    cv=self.cv, scoring=self.score, n_jobs=self.n_jobs) 
       model = model.fit(self.X, self.y)
       return model.best_estimator_, model.best_score_


   def best_RFRegressor(self, estimators):
 
       model = GridSearchCV(RandomForestRegressor(), param_grid=
                            {"n_estimators": estimators.astype('int')},
                    cv=self.cv, scoring=self.score, n_jobs=self.n_jobs) 
       model = model.fit(self.X, self.y)
       return model.best_estimator_, model.best_score_


   def best_GBRegressor(self, estimators):
 
       model = GridSearchCV(GradientBoostingRegressor(), param_grid=
                            {"n_estimators": estimators.astype('int')},
                    cv=self.cv, scoring=self.score, n_jobs=self.n_jobs) 
       model = model.fit(self.X, self.y)
       return model.best_estimator_, model.best_score_


