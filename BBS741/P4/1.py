# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
os.chdir('/Users/shuhao/Desktop/Project-4')

import numpy as np
import matplotlib.pyplot as plt
from hmmlearn import hmm
import sys

cpg = {}
with open('cpg.fa') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            sequence_name = line[1:]
            if sequence_name not in cpg:
                cpg[sequence_name] = ''
            continue
        sequence = line
        cpg[sequence_name] += sequence
cpg_shuffled = {}
with open('cpg.shuffled.fa') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            sequence_name = line[1:]
            if sequence_name not in cpg_shuffled:
                cpg_shuffled[sequence_name] = ''
            continue
        sequence = line
        cpg_shuffled[sequence_name] += sequence

m_cp = hmm.MultinomialHMM(n_components=4)
m_cp.startprob_ = np.array([0.25, 0.25, 0.25, 0.25])
m_cp.transmat_ = np.array([[0.18, 0.274, 0.426, 0.12],[0.171, 0.368, 0.273, 0.188],[0.161, 0.339, 0.375, 0.125],[0.079, 0.355, 0.384, 0.182]])
m_cp.emissionprob_ = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

m = hmm.MultinomialHMM(n_components=4)
m.startprob_ = np.array([0.25, 0.25, 0.25, 0.25])
m.transmat_ = np.array([[0.3, 0.205, 0.285, 0.21],[0.322, 0.298, 0.078, 0.302],[0.248, 0.246, 0.298, 0.208],[0.177, 0.239, 0.292, 0.292]])
m.emissionprob_ = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

s = [value.lower() for key, value in cpg.items()][:500]
for i, seq in enumerate(s): 
    seq = seq.replace('c', 'b');
    seq = seq.replace('g', 'c');
    seq = seq.replace('t', 'd');
    s[i] = seq;
lr = [];
for seq in s: 
    x = []; 
    for nt in seq:
        x.append([ord(nt) - 97])
    lr.append((m_cp.score(x)-m.score(x))/len(seq))

s2 = [value.lower() for key, value in cpg_shuffled.items()][:500]
for i, seq in enumerate(s2):    
    seq = seq.replace('c', 'b');
    seq = seq.replace('g', 'c');
    seq = seq.replace('t', 'd');
    s2[i] = seq;
lr2 = [];
for seq in s2: 
    x = []; 
    for nt in seq:
        if nt == 'n':
            continue
        x.append([ord(nt) - 97])
    lr2.append((m_cp.score(x)-m.score(x))/len(seq))

import matplotlib.pyplot as plt
plt.hist(lr, bins='auto', label='cpg')
plt.hist(lr2, bins='auto', label='cpg_shuffled') 
plt.legend(loc='upper left')
plt.title("log likelihood ratios")
plt.show()

###Conclusion: The two sets separate to 2 distinct classes. -0.03 seems a good threshold 
#that can be used to identify CpG islands.





