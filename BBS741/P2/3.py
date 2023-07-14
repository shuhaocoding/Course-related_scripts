#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 01:38:23 2018

@author: shuhao
"""
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

data = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", skiprows=1, usecols=range(2,128))
CTCFdata = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", skiprows=1, usecols=1)
data = sm.add_constant(data)
model = sm.OLS(CTCFdata,data)
results = model.fit()

p = results.pvalues.tolist()[1:]
p_s = sorted(p)
I10 = [p.index(p_s[x])+2 for x in range(10)]
newdata = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", skiprows=1, usecols=I10)

newdata = sm.add_constant(newdata)
model = sm.OLS(CTCFdata,newdata)
model3 = model.fit()
print(model3.summary()) # R2 is 0.697