# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 18:40:11 2021

@author: u0139894
"""
import os
import logging
logging.getLogger("cobra").setLevel(logging.ERROR)
import numpy as np
import sys

from pathlib import Path
path = Path.cwd()

base_path  = "/".join(os.path.abspath(__file__).split("/")[:-2])
sys.path.insert(0, base_path)


from dnngior.reaction_class import Reaction

from dnngior.build_model import *

from pathlib import Path
path = Path.cwd()

import dnngior
import cobra
import os

import gurobipy as gp
from gurobipy import GRB


# Example 1. Gapfilling a model using a complete medium
# -----------------------------------------------------

base_path  = "/".join(os.path.abspath(__file__).split("/")[:-2])
draftModel = os.path.join(base_path, "docs/models/bh_ungapfilled_model.sbml")

gapfill_compl      = dnngior.Gapfill(draftModel, medium = None, objectiveName = 'bio1')
gf_model_compl_med = gapfill_compl.gapfilledModel.copy()

print("NN gapfilling added {} new readctions".format(len(gapfill_compl.added_reactions)))
print("The NN gapfilled model, comes with {} reactions and {} metabolites".format(len(gf_model_compl_med.metabolites), len(gf_model_compl_med.reactions)))

for reaction in gapfill_compl.added_reactions:
    print(reaction)

# Example 2. Gapfilling a model using a defined medium
# ------------------------------------------------------

media_file = os.path.join(base_path, 'docs/biochemistry/Nitrogen-Nitrite_media.tsv')

Nit_media = {}
with open(media_file) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        Nit_media['EX_' + a[0] + '_e0'] = {'lower_bound':-1, 'upper_bound':1, 'metabolites':{a[0]+'_e0':-1.0}}

gapfill_nitr     = dnngior.Gapfill(draftModel, medium = Nit_media, objectiveName = 'bio1')
gf_model_Nit_med = gapfill_nitr.gapfilledModel.copy()

print("NN gapfilling added {} new readctions".format(len(gapfill_nitr.added_reactions)))
print("The NN gapfilled model, comes with {} reactions and {} metabolites".format(len(gf_model_Nit_med.metabolites), len(gf_model_Nit_med.reactions)))

# 
# for reaction in gapfill_nitr.added_reactions:
#     print(reaction)
