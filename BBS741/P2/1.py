# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

data = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", skiprows=1, usecols=range(2,128))
CTCFdata = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", skiprows=1, usecols=1)
data = sm.add_constant(data)
model = sm.OLS(CTCFdata,data)
model1 = model.fit()
print(model1.summary())
# R2 is 0.733

p = model1.pvalues.tolist()[1:]
p_s = sorted(p)
TFnames = np.loadtxt('TF-Matrix.chr1.txt', delimiter="\t", dtype='str')[0,2:128]
features = [TFnames[p.index(p_s[x])] for x in range(5)]
pVals = p_s[:5]

#The top 5 important features are'TRIM22', 'RAD21', 'SMC3', 'EZH2', 'ZNF143'
#Among those TFs, 'RAD21', 'SMC3', 'EZH2', 'ZNF143' are all to some extent 
#associated with chromatin remodeling. Since CTCF is itself a component of insulator 
#complex and responsible for chromatin remodeling, those TFs are expected to predict 
#CTCF binding. However, 'TRIM22' is an E3 ubiquitin ligase, which seems to be 
#irrelevant to chromatin remodeling. Thus, the link between 'TRIM22' and CTCF
#requires further investigation.  
