#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 17:33:06 2019

@author: shuhao
"""
# Predict states for sequences from cpg.withbg
data = np.loadtxt('cpg.withbg.bed', delimiter="\t", usecols=(1,2,5,6), dtype='int')[:5]
data = np.array([row - row[0] for row in data])
st = [0] * 5; 
for i in range(5):
    st[i] = [0] * data[i,1]
    for j in range(data[i,2]):
        st[i][j]=0
    for j in range(data[i,2],data[i,3]):
        st[i][j]=1
    for j in range(data[i,3],data[i,1]):
        st[i][j]=0
cpg_st = {'island217': st[0], 'island218': st[1], 'island219': st[2], 'island220': st[3], 'island221': st[4]}

cpg_withbg = {}
with open('cpg.withbg.fa') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            sequence_name = line[1:]
            if sequence_name not in cpg_withbg:
                cpg_withbg[sequence_name] = ''
            continue
        sequence = line
        cpg_withbg[sequence_name] += sequence
cpg_sub = {'island217': cpg_withbg['island217'], 'island218': cpg_withbg['island218'], 'island219': cpg_withbg['island219'], 'island220': cpg_withbg['island220'], 'island221': cpg_withbg['island221']}

count = 0; l = 0; predst = [];
for key, value in sorted(cpg_sub.items()):
    seq = value.lower(); 
    seq = seq.replace('c', 'b');
    seq = seq.replace('g', 'c');
    seq = seq.replace('t', 'd');
    seq = seq.replace('n', '');
    x = []
    for nt in seq:
        x.append([ord(nt) - 97])
    Z = model.predict(x)
    predst.append(Z)
    for i in range(len(seq)):
        if Z[i] == cpg_st[key][i]:
            count += 1
    l += len(seq)
p = count/float(l)
print p
### Around 35.8% of the states predicted by HMM are correct for the subset I chose.

# Predict states for sequences from cpg.test
cpg_test = {}
with open('cpg.test.fa') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            sequence_name = line[1:]
            if sequence_name not in cpg_test:
                cpg_test[sequence_name] = ''
            continue
        sequence = line
        cpg_test[sequence_name] += sequence

cpg_coor2 = {}; l = 0; count = 0;
for key, value in sorted(cpg_test.items()):
    seq = value.lower(); 
    seq = seq.replace('c', 'b');
    seq = seq.replace('g', 'c');
    seq = seq.replace('t', 'd');
    seq = seq.replace('n', '');
    x = []
    for nt in seq:
        x.append([ord(nt) - 97])
    Z = model.predict(x)
    c = [];
    for i in range(len(Z)):
        if Z[i] == 0:
            continue
        if i != 0 and i != len(Z)-1 and Z[i]==Z[i-1] and Z[i]==Z[i+1]:
            continue
        c.append(i+1)
    c2 = [];
    for i in xrange(0, len(c), 2):
        c2.append((c[i],c[i+1]))
    cpg_coor2[key] = c2
    st = [0] * len(seq)
    for i in range(cpg_num[key]):
        lower = cpg_coor[key][i][0]-1
        upper = cpg_coor[key][i][1]
        for j in range(lower,upper):
            st[j]=1
    for i in range(len(seq)):
        if Z[i] == st[i]:
            count += 1
    l += len(seq)
p2 = count/float(l)
print cpg_coor2
print p2
### Coordiates predicted by HMM are stored in cpg_coor2. Around 80% of states are the same
#under the prediction of MM and HMM.
