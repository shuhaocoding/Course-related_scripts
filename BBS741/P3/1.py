# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
P1 = np.loadtxt('ENCFF146XEP.bed', delimiter="\t", dtype='str')
P2 = np.loadtxt('ENCFF553PWA.bed', delimiter="\t", dtype='str')
P3 = np.loadtxt('ENCFF596OJV.bed', delimiter="\t", dtype='str')
P4 = np.loadtxt('ENCFF809YFY.bed', delimiter="\t", dtype='str')

I = []; s = np.vstack((P2,P3,P4))
for i in range(len(P1)):
    for row in s:
        if P1[i,0] == row[0]:
            if min(int(P1[i,2]), int(row[2])) - max(int(P1[i,1]), int(row[1])) > 0:
                I.append(i)
                break
P1n = np.delete(P1, I, 0) # Eliminate peaks that are not unique
SV1 = P1n[:,6].astype(float).tolist(); SV1s = sorted(SV1, reverse=True)
P1n = P1n[[SV1.index(SV1s[x]) for x in range(500)],:] # Take the top 500 peaks ranked by signal values
beg = []; end = []
for i in range(500):
    beg.append((int(P1n[i,1])+int(P1n[i,2]))/2-120)
    end.append((int(P1n[i,1])+int(P1n[i,2]))/2+120)
beg = np.array(beg); end = np.array(end) # Resize each peak to 240bp
P1n[:,1] = beg; P1n[:,2] = end 
# P1n is a 2d array modified from the original ENCFF146XEP.bed file. 

I = []; s = np.vstack((P1,P3,P4))
for i in range(len(P2)):
    for row in s:
        if P2[i,0] == row[0]:
            if min(int(P2[i,2]), int(row[2])) - max(int(P2[i,1]), int(row[1])) > 0:
                I.append(i)
                break
P2n = np.delete(P2, I, 0)
SV1 = P2n[:,6].astype(float).tolist(); SV1s = sorted(SV1, reverse=True)
P2n = P2n[[SV1.index(SV1s[x]) for x in range(500)],:]
beg = []; end = []
for i in range(500):
    beg.append((int(P2n[i,1])+int(P2n[i,2]))/2-120)
    end.append((int(P2n[i,1])+int(P2n[i,2]))/2+120)
beg = np.array(beg); end = np.array(end)
P2n[:,1] = beg; P2n[:,2] = end
# P2n is a 2d array modified from the original ENCFF553PWA.bed file. 

I = []; s = np.vstack((P1,P2,P4))
for i in range(len(P3)):
    for row in s:
        if P3[i,0] == row[0]:
            if min(int(P3[i,2]), int(row[2])) - max(int(P3[i,1]), int(row[1])) > 0:
                I.append(i)
                break
P3n = np.delete(P3, I, 0)
SV1 = P3n[:,6].astype(float).tolist(); SV1s = sorted(SV1, reverse=True)
P3n = P3n[[SV1.index(SV1s[x]) for x in range(500)],:]
beg = []; end = []
for i in range(500):
    beg.append((int(P3n[i,1])+int(P3n[i,2]))/2-120)
    end.append((int(P3n[i,1])+int(P3n[i,2]))/2+120)
beg = np.array(beg); end = np.array(end)
P3n[:,1] = beg; P3n[:,2] = end
# P3n is a 2d array modified from the original ENCFF596OJV.bed file. 

I = []; s = np.vstack((P1,P2,P3))
for i in range(len(P4)):
    for row in s:
        if P4[i,0] == row[0]:
            if min(int(P4[i,2]), int(row[2])) - max(int(P4[i,1]), int(row[1])) > 0:
                I.append(i)
                break
P4n = np.delete(P4, I, 0)
SV1 = P4n[:,6].astype(float).tolist(); SV1s = sorted(SV1, reverse=True)
P4n = P4n[[SV1.index(SV1s[x]) for x in range(500)],:]
beg = []; end = []
for i in range(500):
    beg.append((int(P4n[i,1])+int(P4n[i,2]))/2-120)
    end.append((int(P4n[i,1])+int(P4n[i,2]))/2+120)
beg = np.array(beg); end = np.array(end)
P4n[:,1] = beg; P4n[:,2] = end
# P4n is a 2d array modified from the original ENCFF809YFY.bed file. 










