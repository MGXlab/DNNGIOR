import dnngior

import os


"""

import gurobipy as gp

# Gurobipy WLS credentials; for how to get one check: 
# https://support.gurobi.com/hc/en-us/articles/13232844297489-How-do-I-set-up-a-Web-License-Service-WLS-client-license
params = {}
params["WLSACCESSID"] = "d5419c87-0d36-4a93-9385-773f5483b3c1"
params["WLSSECRET"] = "afa5d95f-ad0b-4a38-9550-a8913aacb7c0"
params["LICENSEID"] = 964844

gp.Env(params = params)

base_path = "/".join(os.path.abspath(__file__).split("/")[:-3])

# BiGG case
print("gapfill using a BiGG model")

# Draft model
draftModel = os.path.join(base_path, "docs/models/e_coli_core.xml")

# Gapfill it
# gf_object = dnngior.Gapfill(draftModel, dbType="BiGG", objectiveName="BIOMASS_Ecoli_core_w_GAM")

# ModelSEED case
print("gapfill using a ModelSEED model")

# Draft model
draftModel = os.path.join(base_path, "docs/models/e_coli_core_Seed.sbml")

# Gapfill it
# gf_object = dnngior.Gapfill(draftModel, objectiveName="bio1")
"""

print("This workflow only tests that dnngior can be installed properly and that's imported in Python 3.10\
      The above workflow will be able to run only once a full Gurobi license is available.\
      Until the, please, make sure you perform the tests/test.py check before making a PR without any errors.")


