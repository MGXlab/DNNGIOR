# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 17:38:28 2023

@author: danie
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import scipy.stats as sts


import cobra
from cobra import Reaction
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
from sklearn.linear_model import LinearRegression as LR
plt.style.use('seaborn-bright')




def getExpData(modelGenes, expDict):
    reactions = []
    expr = []
    
    for gene in modelGenes:
        reactions.append(modelGenes[gene])
        expr.append(np.array(list(map(floatTry, expDict[gene]))))
                
    expr = np.array(expr).T
    
    trExpr = np.array([sts.yeojohnson(i)[0] for i in expr])
    
    
    return trExpr, reactions



def getGroupComparison(idx1, idx2, reactions, trExp):
    
    
    m1 = np.mean(trExp[idx1], axis=0)
    m2 = np.mean(trExp[idx2], axis=0)
    
    lr = LR()
    lr.fit(m1.reshape(-1,1), m2.reshape(-1,1))
    
    pred = lr.predict(m1.reshape(-1,1)).flatten()
    
    
    
    upDownNeutral = []
    up = []
    down = []
    sig = []
    d1 = m2-pred
    
    
    
    z = (d1-np.mean(d1))/np.std(d1)
    
    for i,v in enumerate(z):
        if abs(v)>1.96:
            sig.append(d1[i])
            if d1[i]<0: 
                upDownNeutral.append(-1)
                up.append(reactions[i])
            elif d1[i]>0:
                upDownNeutral.append(1)
                down.append(reactions[i])
        else:
            upDownNeutral.append(0)
    sig = np.array(sig)
    top = min(sig[sig>0])
    bot = max(sig[sig<0])
    x = np.linspace(0, max(m1)+0.5, 1000)
    y = lr.coef_[0]*x + lr.intercept_
    return m1,m2, x,y, d1, up, down, np.array(upDownNeutral), top, bot

in_modelPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bh', 'gsmms')

in_genExp = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bh','genes')

out_modelPath = os.path.join(in_modelPath, 'efms')

model = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'bh_final.xml'))


ncbi2patric = {}

with open(os.path.join(in_genExp, 'bh_BVBRC_feature.txt') ) as f:
    f.readline()
    for line in f:
        a= line.strip().split('\t')
        
        ncbi2patric[a[4].replace('"', '')] = a[3].replace('"', '')
        


geneNormedExp = {}

with open(os.path.join(in_genExp, 'bh_tpm.txt') ) as f:
    f.readline()
    for line in f:
        a= line.strip().split('\t')
        if a[0] in ncbi2patric:
            geneNormedExp[ncbi2patric[a[0]]] = a[4::]


patric2agora = {}
agora2patric = {}
    

with open(os.path.join(in_modelPath, 'BH_AGORAtoPATRIC.tsv')) as f:
    for line in f:
        a=line.strip().split('\t')
        agora = a[0]
        patric = a[1].split('.')
        patric = patric[0] + '.' + patric[1] + '.' + patric[2] + '.' + patric[3]
        
        agora2patric[agora] = 'fig|' + patric
        patric2agora[patric] = agora



modelGenes = {}
modelReactions = {}


for reac in model.reactions:
    
    
    genes = []
    
    for gene in reac.genes:
        if gene.id in agora2patric:
            genes.append(agora2patric[gene.id])
        else:
            genes.append('fig|' + gene.id)
    
    modelReactions[reac.id] = genes[:]



for reac in modelReactions:
    for gene in modelReactions[reac]:
        if gene not in modelGenes:
            modelGenes[gene] = []
        modelGenes[gene].append(reac)

for gene in modelGenes:
    s = list(set(modelGenes[gene]))
    modelGenes[gene] = ' / '.join(s)




    
trExp, reactions = getExpData(modelGenes, geneNormedExp)

# for i,v in enumerate(reactions):
#     if len(v)>15:
#         reactions[i] = v[0:15] + '...'

tre,glc,x,y, d1, up, down, upDownNeutral, top, bot = getGroupComparison([0,1,2], [3,4,5], reactions, trExp)

cols = []

for i in upDownNeutral:
    if i==1:
        cols.append('#7700a6')
    elif i==-1:
        cols.append('#1afe49')
    else:
        cols.append('#8386f5')
        
cols=np.array(cols)

fig, ax = plt.subplots()


ax.scatter(tre, glc, s=20, facecolor=cols, edgecolor='black', linewidth=.2)
ax.plot(x, y+top,'--', color='k', lw=0.7)
ax.plot(x, y+bot,'--', color='k', lw=0.7)

#axins = ax.inset_axes([0.1, 0.5, 0.5, 0.5])
# texts = [plt.text(tre[i], glc[i], reactions[i], fontsize=6) for i,v in enumerate(upDownNeutral) if v!=0]
# adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))


# ax.annotate('TREpts', xy = (15.60931674678146,7.835653234399946), xytext = (17,4.5), arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))

# ax.annotate('FDH', xy = (9.12208769,3.72577229), xytext = (12,0), arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))

# ax.annotate('POR4', xy = (10.36042541,3.93223249), xytext = (14,2.5), arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))

# ax.annotate('pbiosynthesis', xy = (7.259327759729246,2.5100173075039134), xytext = (23,5), arrowprops=dict(arrowstyle="->",connectionstyle="arc3"), alpha=0)

# ax.annotate('pbiosynthesis', xy = (20.010220086942805,9.869446823158482), xytext = (23,5), arrowprops=dict(arrowstyle="->",connectionstyle="arc3"), alpha=0)



#axins.scatter(tre[upDownNeutral==1], glc[upDownNeutral==1], s=20, facecolor=cols[upDownNeutral==1], edgecolor='black', linewidth=.2)
#axins.plot(x, y+top,'--', color='k', lw=0.7)
#axins.plot(x, y+bot,'--', color='k', lw=0.7)
#axins.set_ylim(0,sum(upDownNeutral==1)*2)
#axins.set_xlim(-10.5,20)
#axins.set_axis_off()

#reactions = np.array(reactions)
#ts = np.linspace(sum(upDownNeutral==1)*2, 0, sum(upDownNeutral==1))



#adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))

