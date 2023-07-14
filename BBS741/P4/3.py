#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 23:09:32 2019

@author: shuhao
"""
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

s = [value.lower() for key, value in cpg.items()][:500]
s2 = [value.lower() for key, value in cpg_shuffled.items()][:500]
s = s + s2
for i, seq in enumerate(s): 
    seq = seq.replace('c', 'b');
    seq = seq.replace('g', 'c');
    seq = seq.replace('t', 'd');
    seq = seq.replace('n', '');
    s[i] = seq;
lengths = []; x = []; 
for seq in s: 
    for nt in seq:
        x.append([ord(nt) - 97])
    lengths.append(len(seq))

model = hmm.MultinomialHMM(n_components=2, init_params="", params="te")
model.emissionprob_ = np.array([[0.25, 0.25, 0.25, 0.25], [0.15, 0.35, 0.35, 0.15]])
model.transmat_ = np.array([[0.95,0.05],[0.05,0.95]])
model.startprob_ = np.array([0.5,0.5])
model.fit(x, lengths)
print model.transmat_
print model.emissionprob_

###Conclusion: I used HMM with 2 hidden states (1:CpG and 0:non-CpG): I initialized all the parameters including
#emission probabilities and transition probabilities. I fixed the initial probabilities and let both emission
#transition probabilities evolve. Final emission probabilities and transition probabilities are stored in
#model.emissionprob_ and model.transmat_.

#Test on the rest of cpg data
s = [value.lower() for key, value in cpg.items()][500:]
for i, seq in enumerate(s): 
    seq = seq.replace('c', 'b');
    seq = seq.replace('g', 'c');
    seq = seq.replace('t', 'd');
    seq = seq.replace('n', '');
    s[i] = seq;
lengths = []; x = [];
for seq in s: 
    for nt in seq:
        x.append([ord(nt) - 97])
    lengths.append(len(seq)) 
Z = model.predict(x, lengths)
p = sum(Z)/float(len(Z))
print p

