# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 17:38:28 2023

@author: danie
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt



import cobra
from cobra import Reaction
import numpy as np


in_modelPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bh', 'gsmms')

in_genExp = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bh','genes')

out_modelPath = os.path.join(in_modelPath, 'efms')

model = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'Blautia_hydrogenotrophica_DSM_10507_DGarza.xml'))


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
modelFreq = {}



for reac in model.reactions:
    modelReactions[reac.id] = [reac.id, reac.name, reac.build_reaction_string(use_metabolite_names=1)]
    
    genes = []
    
    for gene in reac.genes:
        if gene.id in agora2patric:
            genes.append(agora2patric[gene.id])
        else:
            genes.append('fig|' + gene.id)
    
    modelGenes[reac.id] = genes[:]

with open(os.path.join(in_modelPath, 'efms', 'bh_trehalose.tsv')) as f:
    f.readline()
    for line in f:
        a=line.strip().split('\t')
        modelFreq[a[0]] = [a[1]]


with open(os.path.join(in_modelPath, 'efms', 'bh_pyruvate.tsv')) as f:
    f.readline()
    for line in f:
        a=line.strip().split('\t')
        modelFreq[a[0]].append(a[1])
        
with open(os.path.join(in_modelPath, 'efms', 'bh_glucose.tsv')) as f:
    f.readline()
    for line in f:
        a=line.strip().split('\t')
        modelFreq[a[0]].append(a[1])
        

ncbi2patric = {}

with open(os.path.join(in_genExp, 'bh_BVBRC_feature.txt') ) as f:
    f.readline()
    for line in f:
        a= line.strip().split('\t')
        
        ncbi2patric[a[4].replace('"', '')] = a[3].replace('"', '')
        


geneNormedExp = {}

with open(os.path.join(in_genExp, 'bh_normalized_counts.txt') ) as f:
    f.readline()
    for line in f:
        a= line.strip().split('\t')
        if a[0] in ncbi2patric:
            geneNormedExp[ncbi2patric[a[0]]] = a[1::]

    

genePvals = {}

with open(os.path.join(in_genExp, 'bh_tpm.txt') ) as f:
    f.readline()
    for line in f:
        a= line.strip().split('\t')
        if a[0] in ncbi2patric:
            genePvals[ncbi2patric[a[0]]] = a[3::]

with open(os.path.join(in_genExp, 'metabolicReactions.tsv'), 'w') as f:
    f.write('reactionID\tname\tstring\ttrehalose\tpyruvate\tglucose\tgene\tfunction\tbh_t14_A\tbh_t14_B\tbh_t14_C\tbh_t32_A\tbh_t32_B\tbh_t32_C\tbh_t72_A\tbh_t72_B\tbh_t72_C\tbh_t14_A\tbh_t14_B\tbh_t14_C\tbh_t32_A\tbh_t32_B\tbh_t32_C\tbh_t72_A\tbh_t72_B\tbh_t72_C\tl2fc_32_12\tadjP_32_12\tl2fc_32_72\tadjP_32_72\n')
    
    for reac in modelReactions:
        f.write('\t'.join(modelReactions[reac]))
        f.write('\t')
        f.write('\t'.join(modelFreq[reac]))
        f.write('\n')
        
        for gene in modelGenes[reac]:
            f.write('\t\t\t\t\t\t')
            f.write(gene + '\t' + genePvals[gene][0] + '\t' + '\t'.join(geneNormedExp[gene]) + '\t' + '\t'.join(genePvals[gene][1::]) + '\n')
        f.write('\n')
            

