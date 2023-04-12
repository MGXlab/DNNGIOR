import dnngior
from dnngior.reaction_class import Reaction
from dnngior.NN_Predictor import NN

import gurobipy as gp
from gurobipy import GRB

import cobra
import numpy as np
import os, random


params = {}
params["WLSACCESSID"] = "8e3edcd8-f820-4fe3-bef2-12c072984d6a"
params["WLSSECRET"]   = "d23976cf-9364-4193-8c5a-c02109353796"
params["LICENSEID"]   = 962454

gp.Env(params = params)

base_path = "/".join(os.path.abspath(__file__).split("/")[:-3])

# Draft model
draftModel = os.path.join(base_path, "docs/models/bh_ungapfilled_model.sbml")

print(NN)
"""
Add a test with no use of the NN Predictor.
"""
