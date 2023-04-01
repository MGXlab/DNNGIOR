# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 11:00:09 2023

@author: danie
"""

import os, sys
import numpy as np
import cobra

from pathlib import Path
path = Path.cwd()

from modelseedpy import MSBuilder, MSGenome
from modelseedpy.core import msmedia
from modelseedpy.core.rast_client import RastClient

sys.path.append(os.path.join(path.parents[1], 'scripts'))


import build_model, gapfill_function, NN_Predictor
from   reaction_class import Reaction
""" 
To install modelseedpy
Follow instructions from: https://github.com/ModelSEED/ModelSEEDpy 
using the git clone option
"""
from modelseedpy import MSBuilder, MSGenome
from modelseedpy.core import msmedia


"""
To annotate a genome with RAST (the following commands are to be executed on terminal): 
    1. Installation 
        Following the instructions here: https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html we get RAST locally:
        curl -O -L https://github.com/BV-BRC/BV-BRC-CLI/releases/download/1.040/bvbrc-cli-1.040.deb
        sudo dpkg -i bvbrc-cli-1.040.deb 
        sudo apt-get -f install
    
    2. Build a new genome object
    The genetic code used in translating to protein sequences. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information on genetic codes.
    rast-create-genome --scientific-name "Salmonella infantis" --genetic-code 11 --domain Bacteria --contigs *.fna > genome.gto
    3. Run the annotation step
    
    rast-process-genome < *.gto > annotated_genome.gto2
    4. Export as a .faa file 
    
    rast-export-genome protein_fasta < annotated_genome.gto2 > annotated_genome.faa
This script will use the annotated_genome.faa output of the RAST annotation as its input
"""

# Set the path to your genome
print("Build MSGenome object")
patric_genome = MSGenome.from_fasta(os.path.join(path.parents[1], 'files', 'genomes', 'annotated_genome.faa'), split = ' ')

rast = RastClient()
rast.annotate_genome(patric_genome)

print("Start building base model")
base_model = MSBuilder.build_metabolic_model(model_id = "Blautia hydrogenotrophica DSM 10507", 
                                             genome   = patric_genome, 
                                             index    = "0",
                                             classic_biomass = True, 
                                             gapfill_model   = False, 
                                             gapfill_media   = None, 
                                             annotate_with_rast = False,
                                             allow_all_non_grp_reactions = True
                                            )


model_name = "bh_ungapfilled_model.sbml"
cobra.io.write_sbml_model(cobra_model = base_model, filename = os.path.join(path.parents[0], 'models', model_name))

print("An ungapfilled modelSEED reconstructed model in now available.")

# Set paths to 
files_path  = os.path.join(path.parents[0])
models_path = os.path.join(files_path, 'models')
path_to_NN  = os.path.join(files_path, 'NN')
path_to_biochem  = os.path.join(files_path,  'biochemistry', 'reactions.tsv')

# Load NN predictor
trainedNN_MS = NN_Predictor.load_NN( path = os.path.join(path_to_NN, 'NN_MS.h5') )







# Here goes your medium      
def_med = {}

#H2O
def_med["EX_cpd00001_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00001_e0":-1}}
#Phosphate
def_med["EX_cpd00009_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00009_e0":-1}}
#O2
def_med["EX_cpd00007_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00007_e0":-1}}
#CO2
def_med["EX_cpd00011_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00011_e0":-1}}
#NH3
def_med["EX_cpd00013_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00013_e0":-1}}
#D-Glucose
def_med["EX_cpd00027_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00027_e0":-1}}
#Mn2+
def_med["EX_cpd00030_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00030_e0":-1}}
#Zn2+
def_med["EX_cpd00034_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00034_e0":-1}}
#Sulfate
def_med["EX_cpd00048_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00048_e0":-1}}
#Cu2+
def_med["EX_cpd00058_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00058_e0":-1}}
#Ca2+
def_med["EX_cpd00063_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00063_e0":-1}}
#H+
def_med["EX_cpd00067_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00067_e0":-1}}
#Cl-
def_med["EX_cpd00099_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00099_e0":-1}}
#Co2+
def_med["EX_cpd00149_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00149_e0":-1}}
#K+
def_med["EX_cpd00205_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00205_e0":-1}}
#Mg
def_med["EX_cpd00254_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00254_e0":-1}}
#Na+
def_med["EX_cpd00971_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00971_e0":-1}}
#Fe2+
def_med["EX_cpd10515_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd10515_e0":-1}}
#fe3	
def_med["EX_cpd10516_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd10516_e0":-1}}
#core oligosaccharide lipid A_e0
def_med["EX_cpd15432_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd15432_e0":-1}}
#two linked disacharide pentapeptide murein units (uncrosslinked, middle of chain)_e0
def_med["EX_cpd15511_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd15511_e0":-1}}
#Farnesylfarnesylgeraniol_e0
def_med["EX_cpd02557_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd02557_e0":-1}}
#diphosphate_e0
def_med["EX_cpd02229_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd02229_e0":-1}}



draft_model = cobra.io.read_sbml_model(os.path.join(models_path, 'bh_ungapfilled_model.sbml')) 

# Build a Reaction object for the reactions present on your model
# draft_reaction = Reaction( model = os.path.join(models_path, 's_infantis_base_model.sbml') )
draft_reaction = Reaction( model = os.path.join(models_path, 'bh_ungapfilled_model.sbml') )

# Build a Reaction object for the exchange reactions; if you have a defined medium, set the fixed_bounds argument accordingly
exchange_reacs = Reaction( model = os.path.join(models_path, 'exchangeReactions.sbml'),
                          fixed_bounds = def_med)


# Build a Reaction object for the reactions under the biochemistry folder 
db_reactions = Reaction(biochem_input = path_to_biochem)
db_reactions.reactions = db_reactions.add_dict(exchange_reacs.reactions, db_reactions.reactions)

#Combine all reactions into one object
all_reactions = Reaction(fixed_bounds = def_med)

all_reactions.reactions = all_reactions.add_dict(draft_reaction.reactions, db_reactions.reactions)





draft_reaction_ids = set(draft_reaction.reactions)

#####remove the predefined exchange reactions#####
for react in exchange_reacs.reactions:
    if react in draft_reaction_ids:
        draft_reaction_ids.remove(react)
##################################################


p = NN_Predictor.predict( draft_reaction_ids, trainedNN = path_to_NN) 
weights = {}
for i in p:
    weights[i]  = np.round(1-p[i], 10)


#######apply the defined media condition############
for react in exchange_reacs.reactions:
    if react not in def_med:
        all_reactions.reactions[react]["lower_bound"] = 0
####################################################


model_NN_gf, obj, new_reacs = gapfill_function.gapfill( all_reactions, 
                                                        draft_reaction_ids, 
                                                        weights, 
                                                        'bio1', 
                                                        result_selection = 'min_reactions',
                                                        
                                                    )

if model_NN_gf is None:
    #try to make use the complete media with high costs, then check if there are essential exchange
    #reactions that should have been added to your media.
    for react in exchange_reacs.reactions:
        if react not in def_med:
            all_reactions.reactions[react]["lower_bound"] = -1
            weights[react] = 1000
    
    model_NN_gf, obj, new_reacs = gapfill_function.gapfill( all_reactions, 
                                                            draft_reaction_ids, 
                                                            weights, 
                                                            'bio1', 
                                                            result_selection = 'min_reactions',
                                                            
                                                        )
    
    
    

    

ref_model = build_model.refine_model(model_NN_gf, draft_model, unscalled = list(def_med.keys()))

for reaction in ref_model.reactions:
     if not draft_model.reactions.has_id(reaction.id):
            print(reaction.id, "~~", reaction.build_reaction_string(use_metabolite_names = 1), '\n\n')
            
cobra.io.write_sbml_model(ref_model, filename = os.path.join(models_path, "bh_gapfilled_model.sbml"))