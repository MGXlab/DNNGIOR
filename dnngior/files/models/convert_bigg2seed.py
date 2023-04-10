#!/usr/bin/env python

import cobra
import random 
import os 

cobra_config = cobra.Configuration()
cobra_config.solver = 'glpk'

models_path  = os.getcwd()
model_file   = os.path.join(models_path, "e_coli_core.xml")
bigg_model   = cobra.io.read_sbml_model(model_file)

bigg_mets2seed = {}
for i in bigg_model.metabolites:
    bigg_mets2seed[i.id] = i.annotation["seed.compound"]


for i in bigg_model.metabolites:
    # print(i.id)
        
    if isinstance(bigg_mets2seed[i.id],str):
        if i.id[-1] == "c":
            i.id="M_"+bigg_mets2seed[i.id]+"_c0"
        else:
            i.id="M_"+bigg_mets2seed[i.id]+"_e0"
    else:
            seed_ids=bigg_mets2seed[i.id]
            # print(seed_ids)
            if seed_ids[0] in bigg_model.metabolites:
                if i.id[-1] == "c":
                    i.id="M_"+seed_ids[1]+"_c0"
                else:
                    i.id="M_"+seed_ids[1]+"_e0"
            else:
                if i.id[-1] == "c":
                    i.id="M_"+seed_ids[1]+"_c0"
                else:
                    i.id="M_"+seed_ids[1]+"_e0"

bigg_reactions2seed = {}
for i in bigg_model.reactions:
    try:
        bigg_reactions2seed[i.id]=i.annotation['seed.reaction']
    except:
        bigg_reactions2seed[i.id]=None

bigg_reactions2seed['PFL']="rxn00157"
bigg_reactions2seed['FBP']="rxn33589"
bigg_reactions2seed['TKT2']="rxn31366"
bigg_reactions2seed['TKT1']="rxn00785"
bigg_reactions2seed['TALA']="rxn29919"
bigg_reactions2seed['BIOMASS_Ecoli_core_w_GAM']="bio1"
bigg_reactions2seed['RPI']="rxn28141"
bigg_reactions2seed['PGI']="rxn33838"
bigg_reactions2seed['PFK']="rxn33754"
bigg_reactions2seed['FRUpts2']="rxn37578"


for i in bigg_model.reactions:

    if i.id == "SUCDi": 
        bigg_model.remove_reactions([i])

    if isinstance(bigg_reactions2seed[i.id], str):
        i.id = bigg_reactions2seed[i.id]

    else:

        check = True

        while check: 

            r_id = random.choice(bigg_reactions2seed[i.id])

            if r_id not in bigg_model.reactions:
                i.id = r_id
                check = False

print(bigg_model.summary())

cobra.io.write_sbml_model(bigg_model, filename="e_coli_core_Seed.sbml")
