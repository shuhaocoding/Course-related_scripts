import os
os.chdir("/Users/shuhao/Desktop/RNAseq_all/CNS")
import glob
fl = glob.glob("./*.csv")

import pandas as pd
dicl = []
gnls = []
gns = set()
for file in fl:
    dics = pd.read_csv(file, header=0).set_index('name').to_dict('index')
    dicl.append(dics)
    df = pd.read_csv(file, header=0)
    gnls.append(df.name.tolist())
    gns = gns|set(df.name.tolist())

ctd = {}
for gn in gns:
    count = 0
    for gnl in gnls:
        if gn in gnl:
            count+=1
    ctd[gn] = count

import numpy as np
from numpy import mean
expd = {}
chd = {}
for gn in gns:
    expl = []
    chl = []
    for dic in dicl:
        if gn in dic:
            expl.append(dic[gn]['baseMean'])
            chl.append(abs(dic[gn]['log2FoldChange']))
    aveexp = mean(expl)
    avech = mean(chl)
    expd[gn] = aveexp
    chd[gn] = avech

f = open('/Users/shuhao/Desktop/candidate gene list_CNS.txt', 'w')
header = 'name,occurrence,ave_exp,ave_log2ch'
f.write(header+'\n')
for gn in gns:
    f.write(gn+','+str(ctd[gn])+','+str(expd[gn])+','+str(chd[gn])+'\n')
f.close()
