#BBS741
#Kmeans clustering skeleton script

import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


#Part 1 - Read in Enhancer-Brain-Matrix.txt and create numerical dataMatrix

#Need to fill in
dataMatrix = np.loadtxt('Enhancer-Brain-Matrix.txt', delimiter="\t", skiprows=1, usecols=range(1,9))

#Part 2 - Scaling data across timepoints 
standardized_dataMatrix = StandardScaler().fit_transform(dataMatrix.transpose()).transpose()

#Part 3 - Kmeans clustering of data Matrix with 4 clusters
kmeans = KMeans(n_clusters=4, random_state=0).fit(standardized_dataMatrix)

#Array with values assigning each enhancer to a cluster
dataLabels=kmeans.labels_

#Part 4 - Print out matrix to "Kmeans-Clustering.txt"
with open('Enhancer-Brain-Matrix.txt','r') as f:
    heads = f.readline().strip().split("\t")[1:] 
eh = np.loadtxt('Enhancer-Brain-Matrix.txt', delimiter="\t", dtype='str', skiprows=1, usecols=0)  
s = np.c_[eh, dataLabels, dataMatrix]
f = open('Kmeans-Clustering.txt','w')
f.write('Enhancer'+'\t'+'Cluster')
for i in heads:
    f.write('\t'+i)
f.write('\n')
for row in s:
    for i in range(9):
        f.write(row[i]+'\t')
    f.write(row[9]+'\n')
f.close()

#PLot heatmaps for each cluster
s2= np.c_[eh, dataLabels, standardized_dataMatrix]
cat1=[]; cat2=[]; cat3=[]; cat4=[]
for row in s2:
    if row[1] == '0':
        cat1.append(list(row))
    elif row[1] == '1':
        cat2.append(list(row))
    elif row[1] == '2':
        cat3.append(list(row))
    else:
        cat4.append(list(row))
cat1=np.array(cat1)[:,2:].astype(float); 
cat2=np.array(cat2)[:,2:].astype(float); 
cat3=np.array(cat3)[:,2:].astype(float); 
cat4=np.array(cat4)[:,2:].astype(float); 

import seaborn as sns
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(4, 8), dpi=200)
plt.subplot(221)
hp1=sns.heatmap(cat1, xticklabels=False)
plt.title('Class1')
plt.subplot(222)
hp2=sns.heatmap(cat2, xticklabels=False, yticklabels=False)
plt.title('Class2')
plt.subplot(223)
hp3=sns.heatmap(cat3, xticklabels=heads)
plt.title('Class3')
plt.subplot(224)
hp4=sns.heatmap(cat4, xticklabels=heads, yticklabels=False)
plt.title('Class4')
plt.tight_layout()
fig.savefig('heatmaps.png', dpi=200)

# Conclusion: Four classes of enhancers are activated at different stages during 
# embryonic state in forebrain. Class 1 enhancers are activated toward the end
# of the embryonic state. Class 2 enhancers are primarily activated from 13.5 to
# 16.5 days. Class 3 enhancers are activated at the beginning of the embryonic
# state. Class 4 enhancers are primarily activated from 11.5 to 13.5 days. 
