# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 18:08:19 2023

@author: danie
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt



import cobra
from cobra import Reaction
import numpy as np






def compareModels(model1, model2):
    
    obj1 = model1.optimize()
    obj2 = model2.optimize()
    
    d = {}
    for i in model1.reactions:
        
        a = i.flux/obj1.objective_value
        b = model2.reactions.get_by_id(i.id).flux/obj2.objective_value
        
        diff = abs(a) - abs(b)
        
        if np.abs(np.round(diff,3))>0:
            d[i.id] = diff
    
    return d
            

        
in_modelPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', 'gsmms', 'efms')

model1 = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'bh_trehalose.xml'))
model2 = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'bh_glucose.xml'))

d = compareModels(model1, model2)