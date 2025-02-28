# -*- coding: utf-8 -*-
"""
adapted from addExchanges on Nov 23 2023

@author: meine
"""

import os
import sys
import cobra
from pandas import read_csv
from pathlib import Path
from reaction_class import Reaction
from variables import *


blank_model = cobra.Model("BiGG_exchanges")
compounds = read_csv(BIGG_COMPOUNDS, index_col = 0, sep='\t', dtype='str')
compounds['name'].fillna('NoName', inplace=True)
addedMets = []
for met in compounds.index:
    if "_e" in met:
        if met not in addedMets:
            met_name = compounds['name'][met].replace(' ','_')
            exchange = cobra.Reaction("EX_" + met)

            m = cobra.Metabolite(id = met, name = met_name, compartment="e")
            blank_model.add_metabolites(m)
            blank_model.add_boundary(blank_model.metabolites.get_by_id(m.id), type="exchange")
            addedMets.append(met)

model_path = os.path.join(MODELS_PATH, 'BiGG_exchanges.sbml')
cobra.io.write_sbml_model(blank_model, filename=model_path)
