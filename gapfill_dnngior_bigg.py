from dnngior import gapfill_class
import os
import sys
from pathlib import Path
path = Path.cwd()
sys.path.append(path)

file_path = os.path.join(path.parent, 'files')

import cobra
import pandas as pd
import numpy as np

input_folder =  '/home/meine/NOBINFBACKUP/models/plantsphere/carveme/Geneless_removed'
output_folder = '/home/meine/NOBINFBACKUP/models/plantsphere/carveme_dnngior_mm/'

db_path = '/home/meine/Documents/git/DNNGIOR/dnngior/files/biochemistry/bigg_reactions.tsv'
db = pd.read_csv(db_path, sep='\t', index_col=0)
file_path = '/home/meine/Documents/git/transGit/files/'
min_med_ids = pd.read_csv(os.path.join(file_path, 'min_med_plant.csv'))['BiGG_IDs'].values
carbon_ids = pd.read_csv(os.path.join(file_path, 'compound_ids_planttest.csv'))['BiGG_IDs'].values

med = 'EX_'+min_med_ids+'_e'
med_new = set(med)
med_new.add('EX_glc__D_e')
med_new_dic = {}
for i in med_new:
    med_new_dic[i] = {'lower_bound':-100.0, 'upper_bound':100.0, 'metabolites':{i[3:]:-1.0}}
for model_name in os.listdir(input_folder)[:1]:
    gf_model = gapfill_class.Gapfill(os.path.join(input_folder, model_name), medium=med_new_dic, gapfill=False, dbType='BiGG', objectiveName='Growth', default_cost=5)
    #print([gf_model.exchange_reacs.reactions[i] for i in gf_model.exchange_reacs.reactions if gf_model.exchange_reacs.reactions[i]['lower_bound'] < 0])
    gf_model.gapfill(result_selection='min_cost')
    cobra.io.write_sbml_model(cobra_model = gf_model.gapfilledModel, filename = os.path.join(output_folder, 'gf_nn_'+model_name.split('.')[-2][2:]+'.xml'))

gf_model.all_reactions.reactions['Sd']
