# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 15:19:41 2023

@author: danie
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt



import cobra
from cobra import Reaction
import numpy as np

in_modelPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bh', 'gsmms')
out_modelPath = os.path.join(in_modelPath, 'efms')
model = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'Blautia_hydrogenotrophica_DSM_10507_DGarza.xml'))

modelOut = 'bh_final.xml'


patric2agora = {}
agora2patric = {}
    

with open(os.path.join(in_modelPath, 'BH_AGORAtoPATRIC.tsv')) as f:
    for line in f:
        a=line.strip().split('\t')
        agora = a[0]
        patric = a[1].split('.')
        patric = patric[0] + '.' + patric[1] + '.' + patric[2] + '.' + patric[3]
        
        agora2patric[agora] = patric
        patric2agora[patric] = agora
agora2patric['476272.21.peg.1245'] = '476272.21.peg.1245'
agora2patric['476272.21.peg.3366'] = '476272.21.peg.3366'
agora2patric['476272.21.peg.3367'] = '476272.21.peg.3367'
agora2patric['476272.21.peg.3368'] = '476272.21.peg.3368'
agora2patric['476272.21.peg.308'] = '476272.21.peg.308'

#######################fixpyruvate###############
model.reactions.EX_cpd00020_e.remove_from_model()


pyruvate_e = model.metabolites.get_by_id('pyr[e]')
pyruvate_c = model.metabolites.get_by_id('pyr[c]')
proton_e = model.metabolites.get_by_id('h[e]')
proton_c = model.metabolites.get_by_id('h[c]')

reaction = Reaction('EX_cpd00020_e')

reaction.name = 'pyruvate exchange'


reaction.add_metabolites({pyruvate_e:-1})
reaction.lower_bound=-1000
reaction.upper_bound=1000

model.add_reactions([reaction])


model.repair()
model.optimize()

reaction = Reaction('PYRt2r')

reaction.name = 'pyruvate reversible transport via proton symport'


reaction.add_metabolites({pyruvate_e:-1, proton_e:-1, pyruvate_c:1, proton_c:1})
reaction.lower_bound=-1000
reaction.upper_bound=1000

model.add_reactions([reaction])


model.repair()
model.optimize()

###############################################################


#################add the FDH###############################

co2_c = model.metabolites.get_by_id('co2[c]')
formate_c = model.metabolites.get_by_id('for[c]')
nad_c =  model.metabolites.get_by_id('nad[c]')
nadh_c =  model.metabolites.get_by_id('nadh[c]')

reaction = Reaction('FDH')

reaction.name = 'formate:NAD+ oxidoreductase'


reaction.add_metabolites({co2_c:-1, proton_c:-1, nadh_c:-1, formate_c:1, nad_c:1})
reaction.lower_bound=-1000
reaction.upper_bound=1000
reaction.gene_reaction_rule = '( 476272.21.peg.308 )'

model.add_reactions([reaction])

model.repair()
model.optimize()
########################################################

#################FDHfdx################################

model.reactions.FDHfdx.gene_reaction_rule = '( 476272.21.peg.3366 and 476272.21.peg.3367 and 476272.21.peg.3368)'

########################################################



#################ChangeGeneNamesToPatricIDS

for reac in model.reactions:
    if reac.genes != frozenset():
        
        genes = []
        
        for gene in reac.genes:
            genes.append(gene.id)
        genes.sort(key=len, reverse=True)
        rule = reac.gene_reaction_rule[:]
        
        
        
        for gene in genes:
            rule = rule.replace(gene, agora2patric[gene])[:]
        
        reac.gene_reaction_rule = rule[:]

cobra.io.write_sbml_model(model, os.path.join(in_modelPath, modelOut))