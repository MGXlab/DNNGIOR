# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 18:40:11 2021

@author: u0139894
"""

import logging
logging.getLogger("cobra").setLevel(logging.ERROR)

from reaction_class import Reaction
import gapfill_function

from build_model import *

import os
from pathlib import Path
path = Path.cwd()

import numpy as np


# Example 1. Gapfilling a model

# Load the draft model

model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-draft.xml')
draft_model = Reaction(model=model_location) 


draft_reaction_ids = set(draft_model.reactions) #Set of reaction ids in draft model.


# Load the database with reactions used for gapfilling
biochem = os.path.join(path.parent, 'files',  'biochemistry', 'reactions.tsv')

db_reactions = Reaction(biochem_input=biochem)

#Combine all reactions into one object

all_reactions = Reaction()
all_reactions.reactions = all_reactions.add_dict(draft_model.reactions, db_reactions.reactions) 
#add_dict returns all items from both input dictionaries

# Gapfill

gapfilled_model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-draft_gpfilled.xml')

model, obj, new_reacs = gapfill_function.gapfill(all_reactions, draft_reaction_ids, {}, 'bio1', result_selection = 'min_reactions')


model = refine_model(model)
# # Add a different media

# media_file = os.path.join(path.parent, 'files', 'biochemistry', 'Nitrogen-Nitrite_media.tsv')

# Nit_media = {}

# with open(media_file) as f:
#     f.readline()
#     for line in f:
#         a = line.strip().split('\t')
#         Nit_media['EX_' + a[0] + '_e0'] = {'lower_bound':-1, 'upper_bound':1, 'metabolites':{a[0]+'_e0':-1.0}}

# model2, obj2, new_reacs2 = gapfill_function.gapfill(all_reactions, draft_reaction_ids, {}, 'bio1', result_selection = 'min_reactions', medium=Nit_media)


# # Add costs to reactions

# #select random reaction and add random costs between 0 and 1
# candidate_reactions = {}
# for i in  all_reactions.reactions:
#     if np.random.uniform()<0.02:
#         candidate_reactions[i] = np.random.uniform()


# model3, obj3, new_reacs3 = gapfill_function.gapfill(all_reactions, draft_reaction_ids, candidate_reactions, 'bio1', result_selection = 'min_reactions')

