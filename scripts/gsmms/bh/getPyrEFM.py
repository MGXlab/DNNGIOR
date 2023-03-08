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


def applyEFM(model, rlist, binaryV):
    for i,v in enumerate(rlist):
        if binaryV[i]==0:
            model.reactions.get_by_id(v).lower_bound=0.0
        else:
            model.reactions.get_by_id(v).lower_bound=-10.0
            

def getBinaryReactionSet(model, reactionList):
    model.optimize()
    
    binaryV = np.zeros(len(reactionList))
    for i,v in enumerate(reactionList):
        if np.round(model.reactions.get_by_id(v).flux,3)!=0:
            binaryV[i] = 1.0
    return binaryV
    
def buildConsensusEFMs(model, trList, rList, conditions, binaryM):
    model.optimize()
    accepted = []
    for i,v in enumerate(binaryM):
        applyEFM(model, trList, v)
        if all([eval(z) for z in conditions]):
            accepted.append(getBinaryReactionSet(model, rList))
    return accepted
            
    

in_modelPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bh', 'gsmms')
out_modelPath = os.path.join(in_modelPath, 'efms')
model = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'Blautia_hydrogenotrophica_DSM_10507_DGarza.xml'))

modelOut = 'bh_pyruvate.tsv'








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

co2_c = model.metabolites.get_by_id('co2[c]')
formate_c = model.metabolites.get_by_id('for[c]')
nad_c =  model.metabolites.get_by_id('nad[c]')
nadh_c =  model.metabolites.get_by_id('nadh[c]')

reaction = Reaction('FDH')

reaction.name = 'formate:NAD+ oxidoreductase'


reaction.add_metabolites({co2_c:-1, proton_c:-1, nadh_c:-1, formate_c:1, nad_c:1})
reaction.lower_bound=-1000
reaction.upper_bound=1000

model.add_reactions([reaction])

model.repair()
model.optimize()
########################################################

# ############################
# #use as glutamate synthase
# model.reactions.GLUOX.lower_bound=0
# model.reactions.GLUOX.upper_bound=0
# ##################################

###########block succinate secretion#########

model.reactions.get_by_id('EX_succ(e)').upper_bound = 0

###############################################

################Fix LDH direction##################

model.reactions.LDH_D.upper_bound = 0

model.reactions.LDH_L.upper_bound = 0


###########################################
model.reactions.ForTrsp.upper_bound=0
model.reactions.ForTrsp.lower_bound=0


############constrain carbon sources#################


# model.reactions.get_by_id('EX_tre(e)').lower_bound = -10

# model.reactions.get_by_id('EX_glc_D(e)').lower_bound = -10

# model.reactions.get_by_id('EX_cpd00020_e').lower_bound = -10

#################################################


#model.reactions.LDH_D.upper_bound=-3.5
#model.reactions.ACKr.upper_bound=-25
#model.reactions.get_by_id('EX_ac(e)').lower_bound = 25
#model.reactions.get_by_id('EX_lac_D(e)').lower_bound = 3.5

met2keep = ['EX_cpd00020_e']

simulN = 1000





        
        
l = np.array([i for i in model.medium if i not in met2keep])


simuls = np.zeros((simulN,len(l)))

for i in range(simulN):
    print(i)
    a = np.arange(len(l))
    np.random.shuffle(a)
    sorter = np.argsort(a)
    simuls[i] = get_media_vector(l[a], model)[sorter]


uniqSimuls = np.unique(simuls, axis=0)
modelReactions = [i.id for i in model.reactions]
conditions = ["model.reactions.get_by_id('EX_cpd00020_e').flux<0"]

accepted = buildConsensusEFMs(model, l, modelReactions,conditions,uniqSimuls)
accepted = np.array(accepted)

s1 = np.sum(accepted,axis=0)

with open(os.path.join(in_modelPath, 'efms', modelOut), 'w') as f:
    f.write('reaction\tpyruvteMM\n')
    for i,v in enumerate(modelReactions):
        
        f.write(v + '\t' + str(np.round(s1[i]/len(accepted),3)) + '\n')
    