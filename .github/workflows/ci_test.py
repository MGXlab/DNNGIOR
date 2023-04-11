import dnngior
import cobra
import os

import gurobipy as gp
from gurobipy import GRB

params = {}
params["WLSACCESSID"] = "8e3edcd8-f820-4fe3-bef2-12c072984d6a"
params["WLSSECRET"]   = "d23976cf-9364-4193-8c5a-c02109353796"
params["LICENSEID"]   = 962454

gp.Env(params=params) 

base_path  = "/".join(os.path.abspath(__file__).split("/")[:-3])
draftModel = os.path.join(base_path, "dnngior/files/models/e_coli_core_Seed.sbml")

gapfill            = dnngior.Gapfill(draftModel, medium = None, objectiveName = 'bio1')
gf_model_compl_med = gapfill.gapfilledModel.copy()

print("NN gapfilling added {} new readctions".format(len(gapfill.added_reactions)))
print("The NN gapfilled model, comes with {} reactions and {} metabolites".format(len(gf_model_compl_med.metabolites), len(gf_model_compl_med.reactions)))
