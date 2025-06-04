import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
path = Path.cwd()

base_path  = "/".join(os.path.abspath(__file__).split("/")[:-2])
sys.path.insert(0, base_path)
draftModelMS = os.path.join(base_path, "docs/models/E_coli_KTE31_388739.3_draft.sbml")


from dnngior import NN_Trainer 
from dnngior.NN_Predictor import NN
from dnngior.gapfill_class  import Gapfill
#Example 3. training a network

NN_path = os.path.join(path.parent,'docs', 'NN')
data = pd.read_csv(os.path.join(NN_path, 'Sample_reaction_presence.csv'), index_col=0)
network = NN_Trainer.train(data=data, modeltype='ModelSEED',output_path=os.path.join(NN_path,'custom_networks','test.npz'), save=True)
tensor_network = NN_Trainer.train(data=data, modeltype='ModelSEED',return_full_network=True, save=False)
#Custom network
Gapfill(draftModelMS, trainedNNPath=os.path.join(NN_path, "custom_networks","test.npz"))
model = read_sbml_model(draftModelMS)
tensor_network.predict(model)
network.predict(model)

