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

# test4 = ['rxn09157_c0','rxn11946_c0','rxn31732_c0','rxn05818_c0','rxn45357_c0','rxn05296_c0', 'rxn11598_c0', 'rxn06646_c0','rxn45894_c0','rxn01179_c0','rxn31598_c0','rxn31601_c0','rxn31598_c0','rxn47459_c0','rxn07061_c0','rxn06998_c0','rxn32035_c0','rxn45209_c0','rxn44224_c0', 'rxn07589_c0','rxn16020_c0','rxn43319_c0','rxn05029_c0']
# ess_mass_imb = ["rxn03386_c0"]
#rxn46694_c0
#
#,
#,"rxn03386_c0",,
#
#,
#set_test = set(test4)

input_folder =  '/home/meine/NOBINFBACKUP/models/plantsphere/ungapfilled/'
output_folder = '/home/meine/NOBINFBACKUP/models/plantsphere/gapfilled/nn_mm_dir'
db_path = '/home/meine/Documents/git/DNNGIOR/dnngior/files/biochemistry/reactions.tsv'
db = pd.read_csv(db_path, sep='\t', index_col=0)
file_path = '/home/meine/Documents/git/transGit/files/'
min_med_ids = pd.read_csv(os.path.join(file_path, 'min_med_plant_test.csv'))['ModelSEED_IDs'].values
carbon_ids = pd.read_csv(os.path.join(file_path, 'compound_ids_planttest.csv'))['ModelSEED_IDs'].values
med = 'EX_'+min_med_ids+'_e0'
ess_qm = ['EX_cpd15432_e0', 'EX_cpd00027_e0']
carbon_ex = 'EX_'+carbon_ids+'_e0'
# car_aa_all = ['EX_cpd00035_e0','EX_cpd00039_e0','EX_cpd00051_e0','EX_cpd00054_e0','EX_cpd00060_e0','EX_cpd00064_e0','EX_cpd00065_e0','EX_cpd00069_e0','EX_cpd00084_e0','EX_cpd00107_e0','EX_cpd00119_e0','EX_cpd00116_e0','EX_cpd00156_e0','EX_cpd00161_e0','EX_cpd00322_e0']
# # car_aa = ['EX_cpd00039_e0','EX_cpd00107_e0','EX_cpd00156_e0','EX_cpd00322_e0','EX_cpd00116_e0']
# car_woaa = set(carbon_ex) - set(car_aa_all)
# len(db) - len(db[db['status'].str.contains('MI')])
bad_reactions = [i+'_c0' for i in db.index[~db['status'].str.contains("OK")]]
mi_reactions = [i+'_c0' for i in db.index[db['status'].str.contains("MI")]]
black_list = ['rxn11062_c0','rxn42178_c0','rxn05017_c0','rxn40445_c0','rxn42091','rxn47890','rxn39398_c0', 'rxn21619_c0', 'rxn21618_c0', 'rxn31418_c0', 'rxn03190_c0', 'rxn45845_c0', 'rxn21663', 'rxn41716_c0','rxn45646_c0']
med_new = set(med)|set(ess_qm)
med_new_dic = {}
for i in med_new:
    med_new_dic[i] = {'lower_bound':-1.0, 'upper_bound':1.0, 'metabolites':{i[3:]:-1.0}}

lom = os.listdir(input_folder)
len_lom = len(lom)
track_bad = {}
gapfill_data = {}
for model_name in lom:
    print('gapfilling: {}, {}/{}'.format(model_name, lom.index(model_name),len_lom))
    gf_model = gapfill_class.Gapfill(os.path.join(input_folder, model_name), medium=med_new_dic, gapfill=False, default_cost=5)
    p_r = set([k for k in gf_model.predicted_reactions if gf_model.predicted_reactions[k]>0.5])
    gf_model.draft_reaction_ids = gf_model.draft_reaction_ids.union(p_r)
    #print([gf_model.exchange_reacs.reactions[i] for i in gf_model.exchange_reacs.reactions if gf_model.exchange_reacs.reactions[i]['lower_bound'] < 0])
    #gf_model.set_weights({i:5.0 for i in gf_model.all_reactions.reactions.keys()})
    for i in bad_reactions:
        if i not in gf_model.draft_reaction_ids:
            gf_model.weights[i] = 1000
    for i in mi_reactions:
        if i not in gf_model.draft_reaction_ids:
            gf_model.weights[i] = 10000

    for i in black_list:
        if i in gf_model.all_reactions.reactions.keys():
            del gf_model.all_reactions.reactions[i]
        if i in gf_model.weights.keys():
            del gf_model.weights[i]

    print("{} in all reactions".format(len(gf_model.all_reactions.reactions)))
    gf_model.gapfilledModel = 42
    gf_model.gapfill(result_selection='min_cost')
    if not gf_model.gapfilledModel == 42:
        gapfill_data[model_name[3:]] = {"added_reactions": {k:gf_model.weights[k] for k in gf_model.added_reactions},"number": len(gf_model.added_reactions)}
        cobra.io.write_sbml_model(cobra_model = gf_model.gapfilledModel, filename = os.path.join(output_folder, 'gf_nn_'+model_name.split('.')[-2][3:]+'.xml'))
        set_of_bad = {}
        set_of_ndbad = set()
        set_of_mi  = {}
        set_of_ndmi= set()
        for i in gf_model.gapfilledModel.reactions.list_attr('id'):
            if i in bad_reactions:
                if i in gf_model.draft_reaction_ids:
                    set_of_bad[i] = '42'
                else:
                    set_of_ndbad.add(i)
                    set_of_bad[i] = gf_model.weights[i]
            if i in mi_reactions:
                if i in gf_model.draft_reaction_ids:
                    set_of_mi[i] = '42'
                else:
                    set_of_ndmi.add(i)
                    set_of_mi[i] = gf_model.weights[i]
        print("mass imbalanced: {}".format(set_of_mi))
        track_bad[model_name[3:]] = {"mi_reactions": set_of_mi,"non draft mi_reactions": set_of_ndmi,"bad reactions":set_of_bad, "number":len(set_of_bad)}
    else:
        track_bad[model_name[3:]] = {"mi_reactions": set(["gapfilling","failed"]),"non draft_mi_reactions": set(), "bad reactions":{}, "number":42}
mi_df = pd.DataFrame(track_bad).T
mi_df.to_csv(file_path+'bad_nn_0812.csv')

gf_df = pd.DataFrame(gapfill_data).T
gf_df.to_csv(file_path+'gf_data_2501.csv')
gf_model.gapfilledModel
