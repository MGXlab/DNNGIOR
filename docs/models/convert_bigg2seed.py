#!/usr/bin/env python

import cobra
import random 
import os, sys 

cobra_config = cobra.Configuration()
cobra_config.solver = 'glpk'

models_path  = os.getcwd()
model_file   = os.path.join(models_path, sys.argv[1])
bigg_model   = cobra.io.read_sbml_model(model_file)


# COMPOUNDS PARTS

bigg_mets2seed = {}
for i in bigg_model.metabolites:
    try:
        bigg_mets2seed[i.id] = i.annotation["seed.compound"]
    except:
        pass

bigg_mets2seed["pheme_e"]    = "cpd00028"
bigg_mets2seed["air_c"]      = "cpd02140"
bigg_mets2seed["2ahhmp_c"]   = "cpd00954"
bigg_mets2seed["acmama_c"]   = "cpd34754"
bigg_mets2seed["pheme_c"]    = "cpd00028"
bigg_mets2seed["25dhpp_c"]   = "cpd00957"
bigg_mets2seed["myrsACP_c"]  = "cpd11466"
bigg_mets2seed["palmACP_c"]  = "cpd11476"
bigg_mets2seed["ocdcaACP_c"] = "cpd15268"
bigg_mets2seed["octeACP_c"]  = "cpd29680"
bigg_mets2seed["o2s_c"]      = ""
bigg_mets2seed["12dgr_HP_c"] = "cpd11427"
bigg_mets2seed["h2co3_e"]    = "cpd00242"
bigg_mets2seed["no_e"]       = "cpd00418"
bigg_mets2seed["no_c"]       = "cpd00418"
bigg_mets2seed["mlthf_c"]    = "cpd00125"
bigg_mets2seed["h2co3_c"]    = "cpd00242"
bigg_mets2seed["fpram_c"]    = "cpd02826"
bigg_mets2seed["kdo_c"]      = "cpd00875"
bigg_mets2seed["ugmd_c"]     = "cpd02964"
bigg_mets2seed["u3aga_HP_c"] = "cpd02886"



for i in bigg_model.metabolites:
    if i.id == "o2s_c" or i.id == "h2co3_c":
        bigg_model.remove_metabolites([i])

    if isinstance(bigg_mets2seed[i.id],str):
        if i.id[-1] == "c":
            i.id = "M_" + bigg_mets2seed[i.id] +"_c0"
        else:
            i.id = "M_" + bigg_mets2seed[i.id] + "_e0"

    else:
            seed_ids = bigg_mets2seed[i.id]

            if seed_ids[0] in bigg_model.metabolites:
                if i.id[-1] == "c":
                    i.id = "M_" + seed_ids[1] + "_c0"
                else:
                    i.id = "M_" + seed_ids[1 ] +"_e0"
            else:
                if i.id[-1] == "c":
                    i.id = "M_" + seed_ids[1] + "_c0"
                else:
                    i.id = "M_" + seed_ids[1] + "_e0"


# REACTIONS PART 

bigg_reactions2seed = {}
bigg_reactions2seed['PFL']="rxn00157"
bigg_reactions2seed['FBP']="rxn33589"
bigg_reactions2seed['TKT2']="rxn31366"
bigg_reactions2seed['TKT1']="rxn00785"
bigg_reactions2seed['TALA']="rxn29919"
bigg_reactions2seed['BIOMASS_HP_published']="bio1"
bigg_reactions2seed['RPI']="rxn28141"
bigg_reactions2seed['PGI']="rxn33838"
bigg_reactions2seed['PFK']="rxn33754"
bigg_reactions2seed['FRUpts2']="rxn37578"
bigg_reactions2seed["CYSS"] = "rxn00649"
bigg_reactions2seed["ADSK"] = "rxn00361"

counter = 0 
for i in bigg_model.reactions:
    try:
        bigg_reactions2seed[i.id]=i.annotation['seed.reaction']
    except:
        # bigg_reactions2seed[i.id] = None
        # bigg_model.remove_reactions([i])
        counter += 1
print(counter)


print("~~~~~")
print(bigg_reactions2seed['BIOMASS_HP_published'])
print("~~~~~")

for i in bigg_model.reactions:

    try: 
        isinstance(bigg_reactions2seed[i.id], str)
    except:
        bigg_model.remove_reactions([i])
        continue

    if isinstance(bigg_reactions2seed[i.id], str):

        try:
            i.id = bigg_reactions2seed[i.id]

        except:
            bigg_model.remove_reactions([i])


    else:

        check = True

        while check: 

            r_id = random.choice(bigg_reactions2seed[i.id])

            if r_id not in bigg_model.reactions:
                i.id = r_id
                check = False


for i in bigg_model.reactions:
    if i.id[:3] != "rxn" or i.id != "bio1":
        bigg_model.remove_reactions([i])



print(bigg_model.summary())

fname = sys.argv[1].split(".")[0] + "_mseed.xml"

cobra.io.write_sbml_model(bigg_model, filename=fname)
