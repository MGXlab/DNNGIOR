# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 11:06:59 2023

@author: danie
"""
import os
from pathlib import Path
import numpy as np
import scipy.stats as sts
import cobra
from sklearn.linear_model import LinearRegression as LR
import scipy.stats as sts
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')


class GeneExpr:
    
    def __init__(self, 
                 geneFolder,
                 tpmFile,
                 groupLabels,
                 groupIDX,
                 groupComparison,
                 featuresFile,
                 sbmlModel):
        
        self.geneFolder = geneFolder
        self.groupLabels = groupLabels
        self.groupIDX = groupIDX #in tpm file
        self.groupIDX_ = []
        counter = -1
        for idxs in self.groupIDX:
            a = []
            for c in idxs:
                counter+=1
                a.append(counter)
            self.groupIDX_.append(a[:])
        self.groupComparison = groupComparison
        self.genes, self.tpm = self.__parseTPM(tpmFile)
        self.__parsePatricFeatures(featuresFile)
        self.means = self.__getMeans(self.tpm)
        self.powerM, self.powerSTD = self.__getPowerT(self.tpm)
        self.pvalsDeseq = self.__getPvalsDeseq()
        self.genes2reactions, self.reactions2genes = self.__parseModelGenes(sbmlModel)
        
        
    def __parseTPM(self, tpmFile):
        with open(os.path.join(self.geneFolder, tpmFile)) as f:
            d = {}
            genes = []
            f.readline()
            for line in f:
                a = line.strip().split('\t')
                genes.append(a[1])
                d[a[1]] = {}
                for i, group in enumerate(self.groupLabels):
                    #the patric id is assumed at position 1 of the tpm file
                    d[a[1]][group] = np.array(list(map(self.floatTry, a[self.groupIDX[i][0]: self.groupIDX[i][-1] + 1])))
        return np.array(genes), d
    
    def __parsePatricFeatures(self, featuresFile):
        
        ncbi2patric = {}
        patric2function = {}
        
        with open(os.path.join(self.geneFolder, featuresFile)) as f:
            f.readline()
            for line in f:
                a = line.strip().split('\t')
                
                ncbi2patric[a[6].replace('"', '')] = a[5].replace('"', '')
                patric2function[a[5].replace('"', '')] = a[12]
        self.ncbi2patric = ncbi2patric
        self.patric2function = patric2function
            
    def __parseModelGenes(self, model):
        g2r = {i:[] for i in self.genes}
        reacs = {}
        
        for reaction in model.reactions:
            reacs[reaction.id] = []
            if reaction.genes != frozenset():
                for gene in reaction.genes:
                    geneID = 'fig|' + gene.id
                    g2r[geneID].append(reaction.id)
                    reacs[reaction.id].append(geneID)
        return g2r, reacs
        
        
        
    def __getGeneMatrix(self, geneD):
        m = []
        
        for gene in self.genes:
            m.append(np.concatenate([geneD[gene][group] for group in self.groupLabels]))
            
        return np.array(m)
    
    def __getGeneDict(self, geneM):
        d = {}
        
        for gi, gene in enumerate(self.genes):
            d[gene] = {}
            for gri, group in enumerate(self.groupLabels):
                d[gene][group] = np.array(geneM[gi][self.groupIDX_[gri][0]: self.groupIDX_[gri][-1] + 1])
                
        return d
    
                
    
    def __getMeans(self, geneD):
        
        return {group:np.array([np.mean(geneD[gene][group]) for gene in self.genes]) for group in self.groupLabels}
    
    def __getStds(self, geneD):
        
        return {group:np.array([np.std(geneD[gene][group]) for gene in self.genes]) for group in self.groupLabels}
    
    
    def __getPowerT(self, geneD):
        expr = self.__getGeneMatrix(geneD).T
        trExpr = np.array([sts.yeojohnson(i)[0] for i in expr]).T
        
        trD = self.__getGeneDict(trExpr)
        
        return self.__getMeans(trD), self.__getStds(trD)
        
    def __getPvalsDeseq(self):
        pvals = {gene:{comp:np.nan for comp in self.groupComparison} for gene in self.genes}
        self.notInPatric = []
        
        for comp in self.groupComparison:
            with open(os.path.join(self.geneFolder, self.groupComparison[comp])) as f:
                f.readline()
                for line in f:
                    a = line.strip().split('\t')
                    geneID = a[0].replace('"', '')
                    if geneID not in self.ncbi2patric:
                        self.notInPatric.append(geneID)
                    else:
                        if abs(self.floatTry(a[2]))>1:
                            pvals[self.ncbi2patric[geneID]][comp] = self.floatTry(a[-1], 1.0)
                        else:
                            pvals[self.ncbi2patric[geneID]][comp] = 1
                    
                    
        
        return pvals
    
    
    
        
    
        
        
                    
    
    @staticmethod
    def floatTry(string, r=0.0):
        try:
            return float(string)
        
        except:
            return r






def extracReactions(exprObj, reactionList, fcL, group):
    x, p, g = [],[],[]
    
    for reac in reactionList:
        genes = np.array(exprObj.reactions2genes[reac])
        
        if len(genes)>0:
            
            change=[]
            pv = []
            gids = []
            
            
            for gene in genes:
                gidx = np.arange(len(exprObj.genes))[exprObj.genes==gene]
                change.append(fcL[gidx])
                pv.append(exprObj.pvalsDeseq[gene][group])
                gids.append(gene)
            
            change = np.array(change).flatten()
            pv = np.array(pv).flatten()
            gids = np.array(gids).flatten()
            
            x.append(change[np.abs(pv)==min(np.abs(pv))][0])
            p.append(pv[np.abs(pv)==min(np.abs(pv))][0])
            g.append(gids[np.abs(pv)==min(np.abs(pv))][0])
        
        else:
            x.append(0)
            p.append(1)
            g.append('NF')

    return np.array(x), np.array(p), np.array(g)
    
    
    

            
#example usage
# geneFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', 'genes')
# modelFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', 'gsmms') 
# model = cobra.io.read_sbml_model(os.path.join(modelFolder, 'bh_final.xml'))

# bhGE = GeneExpr(geneFolder = geneFolder, 
#                 tpmFile = 'bh_tpm.txt', 
#                 groupLabels = ['t14', 't32', 't72'], 
#                 groupIDX = [[4,5,6],
#                             [7,8,9],
#                             [10,11,12]
#                             ], 
#                 groupComparison = {('t14', 't32'):'t14vst32_deseq.txt',
#                                    ('t14', 't72'): 't14vst72_deseq.txt',
#                                    ('t32', 't72'): 't32vst72_deseq.txt'}, 
#                 featuresFile =  'bh_BVBRC_features.txt', 
#                 sbmlModel = model
#                 )   




# reactionList = []

# with open(os.path.join(geneFolder, 'wcReactions.txt')) as f:
#     for line in f:
#         reactionList.append(line.strip())


# reactionList = np.array(reactionList)

# fc = bhGE.powerM['t32'] - bhGE.powerM['t14']

# group = ('t14', 't32')

# x,p,g = extracReactions(bhGE, reactionList, fc, group)

# sorter = np.argsort(x)

# x = x[sorter]
# g = g[sorter]
# reactionList = reactionList[sorter]
# p = p[sorter]


# c= np.array(['#6B6ACF']*len(x))

# c[(p < 0.05) & (x > 0)] = '#3A68AE'

# c[(p < 0.05) & (x < 0)] = '#44A043'


# fig, ax = plt.subplots()

# ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

# ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

# ax.hlines(np.arange(len(x)),-8, 8, color = 'k', lw=.5, ls='--')    


# ax.set_yticks(np.arange(len(x)), labels=reactionList, fontsize=10)
# ax.set_xlim(-max(abs(x))-0.5, max(abs(x))+0.5)
# plt.rcParams["figure.figsize"] = (3 ,10)
# plt.tight_layout()