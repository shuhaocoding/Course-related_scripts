import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.preprocessing import StandardScaler

#Part 1 - Read in Enhancer-Mouse-Matrix.txt and create numerical dataMatrix
dataMatrix = np.loadtxt('Enhancer-Mouse-Matrix.txt', delimiter="\t", skiprows=1, usecols=range(1,73))
with open('Enhancer-Mouse-Matrix.txt','r') as f:
    heads = f.readline().strip().split("\t")[1:]
#Need to fill in

#Part 2 - Scaling data across enhancers and transpose 
dataMatrix=StandardScaler().fit_transform(dataMatrix)
dataMatrix=dataMatrix.transpose()

#Part 3 - Run Clustering 
data_dist = pdist(dataMatrix)
data_link = linkage(data_dist)
dendrogram(data_link, labels=heads)
plt.xlabel('Biosamples')  
plt.ylabel('Distance') 
plt.tight_layout()
plt.savefig('dendrogram.png',dpi=300)

# Conclusion: Biosamples from forebrain, midbrain, hindbrain, neurotube, limb, 
# and facial prominence are clustered together. Biosamples from stomach and 
# intestine are clustered together. Biosamples from other tissues (heart, liver, 
# lung, kidney) are clustered separately.

# Outliers: 'liver embryo 14.5 days', 'liver postnatal 0 days', 'facial prominence
# embryo 11.5 days', 'neural tube embyro 11.5 days', 'neural tube embyro 12.5 days', 
# and 'intestine postnatal 0 days'