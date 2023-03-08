# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:10:05 2023

@author: danie
"""

import os
from pathlib import Path
import numpy as np
import scipy.stats as sts
import cobra
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')

from parseGenExpData import *


geneFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', 'genes')






geneFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', 'genes')
modelFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', 'gsmms') 

figuresFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'bhExperiments')


model = cobra.io.read_sbml_model(os.path.join(modelFolder, 'bh_final.xml'))

bhGE = GeneExpr(geneFolder = geneFolder, 
                tpmFile = 'bh_tpm.txt', 
                groupLabels = ['t14', 't32', 't72'], 
                groupIDX = [[4,5,6],
                            [7,8,9],
                            [10,11,12]
                            ], 
                groupComparison = {('t14', 't32'):'t14vst32_deseq.txt',
                                   ('t14', 't72'): 't14vst72_deseq.txt',
                                   ('t32', 't72'): 't32vst72_deseq.txt'}, 
                featuresFile =  'bh_BVBRC_features.txt', 
                sbmlModel = model
                )   



#########################t14 vs 32 #########################
reactionList = []

with open(os.path.join(geneFolder, 'wcReactions.txt')) as f:
    for line in f:
        reactionList.append(line.strip())


reactionList = np.array(reactionList)

fc = bhGE.powerM['t32'] - bhGE.powerM['t14']

group = ('t14', 't32')

x,p,g = extracReactions(bhGE, reactionList, fc, group)

sorter = np.argsort(x)

x = x[sorter]
g = g[sorter]
reactionList = reactionList[sorter]
p = p[sorter]


c= np.array(['#6B6ACF']*len(x))



c[(p < 0.05) & (x > 0)] = '#001eff'

c[(p < 0.05) & (x < 0)] = '#44A043'





fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-8, 8, color = 'k', lw=.5, ls='--')    

labels = reactionList.copy()

labels[labels=='CODH_ACS'] = 'CODH'
labels[labels=='OOR2r'] = 'OOR2'
labels[labels=='OOR2r'] = 'OOR2'
labels[labels=='PFK(ppi)'] = 'PFK'
labels[labels=='HYDFDNrfdx'] = 'HydABC'


ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)



ax.set_xlim(-8.0, 8.0)
ax.set_xlabel('gene expression FC')

ax.text(-4, len(x) + .2, 't14')
ax.text(4, len(x) + .2, 't32')

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_t14vst32.png'), dpi = 300)



#########################t14 vs 72 #########################

fc = bhGE.powerM['t72'] - bhGE.powerM['t14']

group = ('t14', 't72')

x,p,g = extracReactions(bhGE, reactionList, fc, group)

c= np.array(['#6B6ACF']*len(x))



c[(p < 0.05) & (x > 0)] = '#bd00ff'

c[(p < 0.05) & (x < 0)] = '#44A043'





fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-8, 8, color = 'k', lw=.5, ls='--')    

ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)
ax.set_xlim(-8.0, 8.0)
ax.set_xlabel('gene expression FC')

ax.text(-4, len(x) + .2, 't14')
ax.text(4, len(x) + .2, 't72')

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_t14vst72.png'), dpi = 300)


#########################t32 vs 72 #########################

fc = bhGE.powerM['t72'] - bhGE.powerM['t32']

group = ('t32', 't72')

x,p,g = extracReactions(bhGE, reactionList, fc, group)

c= np.array(['#6B6ACF']*len(x))



c[(p < 0.05) & (x > 0)] = '#bd00ff'

c[(p < 0.05) & (x < 0)] = '#001eff'





fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-8, 8, color = 'k', lw=.5, ls='--')    

ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)
ax.set_xlim(-8.0, 8.0)
ax.set_xlabel('gene expression FC')

ax.text(-4, len(x) + .2, 't32')
ax.text(4, len(x) + .2, 't72')

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_t32vst72.png'), dpi = 300)