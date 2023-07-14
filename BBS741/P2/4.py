#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 14:29:17 2018

@author: shuhao
"""
import numpy as np
from sklearn.linear_model import LinearRegression
import random
import matplotlib.pyplot as plt

cols = [1] + I10
matrix = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", skiprows=1, usecols=cols)
R_sqrs = []
for i in range(10000):
    newmat = np.empty((0,11))
    for i in range(10407):
        newmat = np.vstack((newmat,random.choice(matrix)))
    data = newmat[:,1:]
    CTCFdata = newmat[:,0]
    reg = LinearRegression().fit(data, CTCFdata)
    R_sqrs.append(reg.score(data, CTCFdata))
lower = np.percentile(R_sqrs, 2.5) # output 2.5% percentiles
upper = np.percentile(R_sqrs, 97.5) # output 97.5% percentiles

plt.hist(R_sqrs, bins=30, facecolor='skyblue')
plt.xlabel('R2')
plt.ylabel('Frequency')
plt.axvline(x=lower, linewidth=1, color='r')
plt.axvline(x=upper, linewidth=1, color='r')
plt.show()
