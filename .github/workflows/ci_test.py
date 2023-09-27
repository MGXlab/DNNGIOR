import dnngior
import gurobipy as gp
import os

# Gurobipy WLS credentials; for how to get one check: 
# https://support.gurobi.com/hc/en-us/articles/13232844297489-How-do-I-set-up-a-Web-License-Service-WLS-client-license
params = {}
params["WLSACCESSID"] = "cb73b00e-795f-440e-b0b5-32d7f704cb13"
params["WLSSECRET"]   = "a86dbc31-628d-4f7c-a5b5-a8de7016f29a"
params["LICENSEID"]   = 964844

gp.Env(params = params)


base_path = "/".join(os.path.abspath(__file__).split("/")[:-3])

# BiGG case
print("gapfill using a BiGG model")

# Draft model
draftModel = os.path.join(base_path, "docs/models/e_coli_core.xml")

# Gapfill it
gf_object = dnngior.Gapfill(draftModel, dbType="BiGG", objectiveName="BIOMASS_Ecoli_core_w_GAM")

# ModelSEED case
print("gapfill using a ModelSEED model")

# Draft model
draftModel = os.path.join(base_path, "docs/models/e_coli_core_Seed.sbml")

# Gapfill it
gf_object = dnngior.Gapfill(draftModel, objectiveName="bio1")
