import dnngior
import cobra
import os


import gurobipy as gp
from gurobipy import GRB

params = {}
params["WLSACCESSID"]= "8e3edcd8-f820-4fe3-bef2-12c072984d6a"
params["WLSSECRET"]  = "d23976cf-9364-4193-8c5a-c02109353796"
params["LICENSEID"]  = 962454

with gp.Env(params=params) as env:
    with gp.Model(env=env) as m:


base_path  = os.getcwd()
draftModel = os.path.join(base_path, 'S_infantis_mvl3A.sbml')

gapfill            = dnngior.Gapfill(draftModel, medium = None, objectiveName = 'bio1')
gf_model_compl_med = gapfill.gapfilledModel.copy()

modelname = os.path.join(base_path, "S_infantis_mvl3A.sbml")

cobra.io.write_sbml_model(cobra_model = gf_model_compl_med, filename = modelname)
