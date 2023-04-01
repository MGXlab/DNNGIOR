# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 08:35:26 2023

@author: danie
"""

import cobra
from reaction_class import Reaction
import os
import sys

from pathlib import Path
path = Path.cwd()
sys.path.append(path)

from MSEED_compounds import Compounds
from MSEED_reactions import Reactions

compounds_helper = Compounds()
compounds_dict   = compounds_helper.loadCompounds()



blank_model = cobra.Model("exchangeReactions")


biochem = os.path.join(path.parent, 'files',  'biochemistry', 'reactions.tsv')
db_reactions = Reaction(biochem_input=biochem)

reactionsToadd = []

addedMets = []

for reaction in db_reactions.reactions:
    for met in db_reactions.reactions[reaction]["metabolites"]:
        if "e0" in met: 

            if met not in addedMets:

                mObj = compounds_dict[met.split('_')[0]]

                r = cobra.Reaction("EX_" + met)

                m = cobra.Metabolite(id = met, formula = mObj['formula'], name = mObj['name'], compartment="e0")

                r.add_metabolites({m:-1})

                r.upper_bound = 1000
                r.lower_bound = -1000
                
                reactionsToadd.append(r.copy())

                addedMets.append(met)


blank_model.add_reactions(reactionsToadd)

models = os.path.join(path.parent, 'files',  'models', 'exchangeReactions.sbml')
cobra.io.write_sbml_model(blank_model, filename=models)
