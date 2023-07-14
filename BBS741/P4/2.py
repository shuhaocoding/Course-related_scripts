#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 03:22:37 2019

@author: shuhao
"""
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
        
#Test run to decide the window size
seq = cpg_test['test1'].lower(); 
seq = seq.replace('c', 'b');
seq = seq.replace('g', 'c');
seq = seq.replace('t', 'd');
lr1 = []; w = 50; c1 = range(len(seq)-w+2)[1:];
for i in range(len(seq)-w+1):
        sub = seq[i:i+w]
        x = []; 
        for nt in sub:
            x.append([ord(nt) - 97])
        lr1.append((m_cp.score(x)-m.score(x))/w)
lr2 = []; w = 100; c2 = range(len(seq)-w+2)[1:];
for i in range(len(seq)-w+1):
        sub = seq[i:i+w]
        x = []; 
        for nt in sub:
            x.append([ord(nt) - 97])
        lr2.append((m_cp.score(x)-m.score(x))/w)
lr3 = []; w = 150; c3 = range(len(seq)-w+2)[1:];
for i in range(len(seq)-w+1):
        sub = seq[i:i+w]
        x = []; 
        for nt in sub:
            x.append([ord(nt) - 97])
        lr3.append((m_cp.score(x)-m.score(x))/w)
lr4 = []; w = 300; c4 = range(len(seq)-w+2)[1:];
for i in range(len(seq)-w+1):
        sub = seq[i:i+w]
        x = []; 
        for nt in sub:
            x.append([ord(nt) - 97])
        lr4.append((m_cp.score(x)-m.score(x))/w)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True)
f.suptitle('Window test')
ax1.plot(c1, lr1)
ax2.plot(c2, lr2)
ax3.plot(c3, lr3)
ax4.plot(c4, lr4)
### Window size of 300bp seems good

# Determine the locations of CpG islands
cpg_coor = {}; cpg_num = {}
w = 300; thr = -0.03; 
for key, value in cpg_test.items():
    seq = value.lower(); 
    seq = seq.replace('c', 'b');
    seq = seq.replace('g', 'c');
    seq = seq.replace('t', 'd');
    c = []
    for i in range(len(seq)-w+1):
        sub = seq[i:i+w]
        x = []; 
        for nt in sub:
            x.append([ord(nt) - 97])
        if (m_cp.score(x)-m.score(x))/float(w) > thr:
            c.append((i+1, i+w))
    coor = []
    for begin,end in sorted(c):
        if coor and coor[-1][1] >= begin - 1:
            coor[-1] = (coor[-1][0], end)
        else:
            coor.append((begin, end))
    cpg_coor[key] = coor
    cpg_num[key] = len(coor)
print cpg_coor
print cpg_num

###Conclusion: Numbers of CpG islands in each of the test regions are stored in dict cpg_num.
#Coordinates of the CpG islands are stored in dict cpg_coor.


