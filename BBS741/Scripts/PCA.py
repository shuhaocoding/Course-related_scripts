from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn.decomposition import PCA

#Part 1 - Read in Enhancer-Mouse-Matrix.txt and create numerical dataMatrix
dataMatrix = np.loadtxt('Enhancer-Mouse-Matrix.txt', delimiter="\t", skiprows=1, usecols=range(1,73))
with open('Enhancer-Mouse-Matrix.txt','r') as f:
    heads = f.readline().strip().split("\t")[1:]

#Need to fill in
        
#Part 2 - Scaling data across enhancers and transpose 
dataMatrix=StandardScaler().fit_transform(dataMatrix)
dataMatrix=dataMatrix.transpose()
    
#Part 3 - Run PCA 
pca = PCA(n_components=5, svd_solver='full')

#Array of each PCA component for each biosample
transformedMatrix=pca.fit_transform(dataMatrix)

#Array of variance explained for each component
variance=pca.explained_variance_ratio_
with open('PCA_Variance.txt','w') as f:
    f.write('PC1'+'\t'+str(variance[0])+'\n'+'PC2'+'\t'+str(variance[1])+'\n'+'PC3'+'\t'+str(variance[2])+'\n'+'PC4'+'\t'+str(variance[3])+'\n'+'PC5'+'\t'+str(variance[4]))

Matrix = np.loadtxt('Enhancer-Mouse-Matrix.txt', delimiter="\t", dtype='str')
labels = Matrix.transpose()[1:,0]
s = np.c_[labels, transformedMatrix]
f = open('PCA_PC_Matrix.txt','w')
f.write(' '+'\t'+'PC1'+'\t'+'PC2'+'\t'+'PC3'+'\t'+'PC4'+'\t'+'PC5'+'\n')
for row in s:
    f.write(row[0]+'\t'+row[1]+'\t'+row[2]+'\t'+row[3]+'\t'+row[4]+'\t'+row[5]+'\n')
f.close()

#Part 4 - Make scatterplot of PC1 vs PC2 OR print out components for plotting in R

# PC1 vs PC2 colored by tissue type
import matplotlib.pyplot as plt
import matplotlib

fax=[]; fay=[]; fbx=[]; fby=[]; mbx=[]; mby=[]; hbx=[]; hby=[]; 
htx=[]; hty=[]; itx=[]; ity=[]; kdx=[]; kdy=[]; lbx=[]; lby=[];
lvx=[]; lvy=[]; lux=[]; luy=[]; ntx=[]; nty=[]; stx=[]; sty=[];

for row in s:
     PC1 = float(row[1]); PC2 = float(row[2]);
     if 'facial' in row[0]:
         fax.append(PC1); fay.append(PC2);     
     elif 'forebrain' in row[0]:
         fbx.append(PC1); fby.append(PC2)
     elif 'midbrain' in row[0]:
         mbx.append(PC1); mby.append(PC2)
     elif 'hindbrain' in row[0]:
         hbx.append(PC1); hby.append(PC2)
     elif 'heart' in row[0]:
         htx.append(PC1); hty.append(PC2)
     elif 'intestine' in row[0]:
         itx.append(PC1); ity.append(PC2)
     elif 'kidney' in row[0]:
         kdx.append(PC1); kdy.append(PC2)
     elif 'limb' in row[0]:
         lbx.append(PC1); lby.append(PC2)
     elif 'liver' in row[0]:
         lvx.append(PC1); lvy.append(PC2)
     elif 'lung' in row[0]:
         lux.append(PC1); luy.append(PC2)
     elif 'neural' in row[0]:
         ntx.append(PC1); nty.append(PC2)
     else:
         stx.append(PC1); sty.append(PC2)
         
fa = (np.array(fax), np.array(fay));
fb = (np.array(fbx), np.array(fby));
mb = (np.array(mbx), np.array(mby));
hb = (np.array(hbx), np.array(hby));
ht = (np.array(htx), np.array(hty));
it = (np.array(itx), np.array(ity));
kd = (np.array(kdx), np.array(kdy));
lb = (np.array(lbx), np.array(lby));
lv = (np.array(lvx), np.array(lvy));
lu = (np.array(lux), np.array(luy));
nt = (np.array(ntx), np.array(nty));
st = (np.array(stx), np.array(sty));
      
data = (fa, fb, mb, hb, ht, it, kd, lb, lv, lu, nt, st)
colors = matplotlib.cm.rainbow(np.linspace(0, 1, 12))
groups = ('facial', 'forebrain', 'midbrain', 'hindbrain', 'heart', 'intestine', 'kidney', 'limb', 'liver', 'lung', 'neural', 'stomach') 

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for data, color, group in zip(data, colors, groups):
    x, y = data
    ax.scatter(x, y, alpha=0.8, c=color, edgecolors='none', s=30, label=group)
plt.title('Matplot scatter plot')
plt.legend(loc='best', fontsize='small')
fig.savefig('PCA1.png', dpi=300) 


# PC1 vs PC2 colored by time point 
g0x=[]; g0y=[]; g10_5x=[]; g10_5y=[]; g11_5x=[]; g11_5y=[]; g12_5x=[]; g12_5y=[]; 
g13_5x=[]; g13_5y=[]; g14_5x=[]; g14_5y=[]; g15_5x=[]; g15_5y=[]; g16_5x=[]; g16_5y=[];
for row in s:
     PC1 = float(row[1]); PC2 = float(row[2]);
     if '_0' in row[0]:
         g0x.append(PC1); g0y.append(PC2);     
     elif '10.5' in row[0]:
         g10_5x.append(PC1); g10_5y.append(PC2)
     elif '11.5' in row[0]:
         g11_5x.append(PC1); g11_5y.append(PC2)
     elif '12.5' in row[0]:
         g12_5x.append(PC1); g12_5y.append(PC2)
     elif '13.5' in row[0]:
         g13_5x.append(PC1); g13_5y.append(PC2)
     elif '14.5' in row[0]:
         g14_5x.append(PC1); g14_5y.append(PC2)
     elif '15.5' in row[0]:
         g15_5x.append(PC1); g15_5y.append(PC2)
     else:
         g16_5x.append(PC1); g16_5y.append(PC2)
         
g0 = (np.array(g0x), np.array(g0y));
g10_5 = (np.array(g10_5x), np.array(g10_5y));
g11_5 = (np.array(g11_5x), np.array(g11_5y));
g12_5 = (np.array(g12_5x), np.array(g12_5y));
g13_5 = (np.array(g13_5x), np.array(g13_5y));
g14_5 = (np.array(g14_5x), np.array(g14_5y));
g15_5 = (np.array(g15_5x), np.array(g15_5y));
g16_5 = (np.array(g16_5x), np.array(g16_5y));

data = (g0, g10_5, g11_5, g12_5, g13_5, g14_5, g15_5, g16_5)
colors = matplotlib.cm.rainbow(np.linspace(0, 1, 8))
groups = ("0d", "10.5d", "11.5d", "12.5d", "13.5d", "14.5d", "15.5d", "16.5d") 

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for data, color, group in zip(data, colors, groups):
    x, y = data
    ax.scatter(x, y, alpha=0.8, c=color, edgecolors='none', s=30, label=group)
plt.title('Matplot scatter plot')
plt.legend(loc='best', fontsize='small')
fig.savefig('PCA2.png', dpi=300) 
      
# Conclusion: In PCA1, biosamples from liver are clustered together. Biosamples
# from limb and facial prominence are clustered together. Biosamples from forebrain,
# midbrain, hindbrain, and neural tubes are clustered together. All other biosamples
# are clustered together.

# The clusters tend to separate by tissue types instead of time points.