# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 14:42:12 2022

@author: u0139894
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt



import cobra
from cobra import Reaction
import numpy as np


def removeReaction(model, reactionID, optimal, threshold = 0.00001, low=1, up=0):
    ub = model.reactions.get_by_id(reactionID).upper_bound
    
    lb = model.reactions.get_by_id(reactionID).lower_bound
    
    if up:
        model.reactions.get_by_id(reactionID).upper_bound = 0
    if low:
        model.reactions.get_by_id(reactionID).lower_bound = 0
    
    if model.slim_optimize() < threshold*optimal:
        model.reactions.get_by_id(reactionID).upper_bound = ub
        model.reactions.get_by_id(reactionID).lower_bound = lb
        return 1
    return 0


def get_media_vector(toDel, model):
    for i in toDel:
        model.reactions.get_by_id(i).lower_bound = -1000
    sol = model.optimize()
    b = np.zeros(len(toDel))
    for i,v in enumerate(toDel):
        b[i] = removeReaction(model, v, optimal= sol.objective_value, threshold=0.00001, low=1, up=0)
    
    return b


in_modelPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bh', 'gsmms')
out_modelPath = os.path.join(in_modelPath, 'efms')
model = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'Blautia_hydrogenotrophica_DSM_10507_DGarza.xml'))

modelOut = 'bh_pyruvate.xml'





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




# ############################
# #use as glutamate synthase
# model.reactions.GLUOX.lower_bound=0
# model.reactions.GLUOX.upper_bound=0
# ##################################

###########block succinate secretion#########

model.reactions.get_by_id('EX_succ(e)').upper_bound = 0

model.reactions.get_by_id('EX_for(e)').upper_bound = 0

####################################################

################Fix LDH direction##################

model.reactions.LDH_D.upper_bound = 0

model.reactions.LDH_L.upper_bound = 0


############constrain carbon sources#################


model.reactions.get_by_id('EX_tre(e)').lower_bound = -10

model.reactions.get_by_id('EX_glc_D(e)').lower_bound = -10

model.reactions.get_by_id('EX_cpd00020_e').lower_bound = -10

#################################################





# #add the mock objective reaction
prot_c = model.metabolites.get_by_id('h[c]').copy()
nad_c = model.metabolites.get_by_id('nad[c]').copy()
nadh_c = model.metabolites.get_by_id('nadh[c]').copy()
atp_c = model.metabolites.get_by_id('atp[c]').copy()
adp_c = model.metabolites.get_by_id('adp[c]').copy()
accoA_c = model.metabolites.get_by_id('accoa[c]').copy()
coA_c = model.metabolites.get_by_id('coa[c]').copy()
biomass_c = model.metabolites.get_by_id('biomass[c]').copy()
biomass = Reaction('mock_biomass')
biomass.name='Mock biomass function'
biomass.lower_bound=0
biomass.upper_bound=1000
biomass.add_metabolites({atp_c:-3, 
                          accoA_c:-2,
                          nadh_c:-2, 
                          prot_c:-2, 
                          adp_c:3,
                          nad_c:2,
                          coA_c:2
                          })

model.add_reactions([biomass])
model.reactions.biomass525.objective_coefficient = 0

model.reactions.mock_biomass.objective_coefficient = 1

optM = model.optimize()

print(optM.objective_value)

###################


met2keep = ['EX_glc_D(e)', 'EX_glu_L(e)', 'EX_co2(e)']

simulN = 1000





        
        
l = np.array([i for i in model.medium if i not in met2keep])


simuls = np.zeros((simulN,len(l)))

for i in range(simulN):
    print(i)
    a = np.arange(len(l))
    np.random.shuffle(a)
    sorter = np.argsort(a)
    simuls[i] = get_media_vector(l[a], model)[sorter]


s = np.sum(simuls, axis=0)


for i,v in enumerate(l):
    if s[i] == 0:
        model.reactions.get_by_id(v).lower_bound = 0


cobra.io.write_sbml_model(model, os.path.join(out_modelPath, modelOut))
