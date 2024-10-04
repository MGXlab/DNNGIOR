# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 18:40:11 2021

@author: u0139894
"""
import os
import sys
import cobra
import logging
logging.getLogger("cobra").setLevel(logging.ERROR)
import numpy as np
import pandas as pd

from pathlib import Path
path = Path.cwd()

base_path  = "/".join(os.path.abspath(__file__).split("/")[:-2])
sys.path.insert(0, base_path)

from dnngior import gapfill_function
from dnngior.gapfill_class  import Gapfill
from dnngior.reaction_class import Reaction
from dnngior.build_model import *
from dnngior import NN_Trainer
import gurobipy as gp
from gurobipy import GRB


# Example 1. Gapfilling a model using a complete medium
# -----------------------------------------------------


draftModelMS = os.path.join(base_path, "docs/models/E_coli_KTE31_388739.3_draft.sbml")
draftModelBiGG = os.path.join(base_path, "docs/models/bigg_example.xml")

grey_list = ['rxn11062_c0','rxn42178_c0','rxn05017_c0','rxn40445_c0','rxn42091','rxn47890','rxn39398_c0', 'rxn21619_c0', 'rxn21618_c0', 'rxn31418_c0', 'rxn03190_c0', 'rxn45845_c0', 'rxn21663', 'rxn41716_c0','rxn45646_c0']


gapfill_compl_ms      = Gapfill(draftModelMS, medium = None, objectiveName = 'bio1', grey_list=grey_list)
gapfill_compl_bg      = Gapfill(draftModelBiGG, objectiveName='Growth', dbType = 'BiGG')


# Example 2. Gapfilling a model using a defined medium
# ------------------------------------------------------
Nit_media_file = os.path.join(base_path, 'docs/biochemistry/Nitrogen-Nitrite_media.tsv')
gapfill_nitr     = Gapfill(draftModelMS, medium_file = Nit_media_file)

#Example 3. training a network

file_path = os.path.join(path.parent,'docs', 'NN')
data = pd.read_csv(os.path.join(file_path, 'Sample_reaction_presence.csv'), index_col=0)
network = NN_Trainer.train(data=data, modeltype='ModelSEED',output_path=os.path.join(file_path,'custom_networks','test.npz'), save=True)

#
# for reaction in gapfill_nitr.added_reactions:
#     print(reaction)
