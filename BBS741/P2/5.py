#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 16:55:29 2018

@author: shuhao
"""
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LassoCV
from sklearn.metrics import mean_squared_error

CTCFdata = np.loadtxt('TF-Matrix.chr2.txt', delimiter="\t", skiprows=1, usecols=1)
data = np.loadtxt('TF-Matrix.chr2.txt', delimiter="\t", skiprows=1, usecols=range(2,128))
newdata = np.loadtxt('TF-Matrix.chr2.txt', delimiter="\t", skiprows=1, usecols=I10)
OLSdata1 = sm.add_constant(data)
OLSdata2 = sm.add_constant(newdata)

MSE1 = mean_squared_error(CTCFdata, model1.predict(OLSdata1))
MSE2 = mean_squared_error(CTCFdata, model2.predict(data))
MSE3 = mean_squared_error(CTCFdata, model3.predict(OLSdata2))
MSE = [MSE1, MSE2, MSE3] # MSEs of the 3 models

# LASSO model has the best performance.