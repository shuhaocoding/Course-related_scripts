#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 01:03:14 2018

@author: shuhao
"""
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LassoCV
from sklearn.datasets import make_regression

data = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", skiprows=1, usecols=range(2,128))
CTCFdata = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", skiprows=1, usecols=1)
model2 = LassoCV(cv=10).fit(data, CTCFdata)
R2 = model2.score(data, CTCFdata) # R2 is 0.729
p = model2.coef_.tolist()[1:]
n = p.count(0) # 33 features were dropped
