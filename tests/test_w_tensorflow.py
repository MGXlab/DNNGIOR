from dnngior import NN_Trainer
from dnngior.gapfill_class  import Gapfill
#Example 3. training a network

NN_path = os.path.join(path.parent,'docs', 'NN')
data = pd.read_csv(os.path.join(NN_path, 'Sample_reaction_presence.csv'), index_col=0)
network = NN_Trainer.train(data=data, modeltype='ModelSEED',output_path=os.path.join(NN_path,'custom_networks','test.npz'), save=True)
tensor_network = NN_Trainer.train(data=data, modeltype='ModelSEED',return_full_network=True, save=False)
#Custom network
Gapfill(draftModelMS, trainedNNPath=os.path.join(NN_path, "custom_networks","test.npz"))
model = cobra.io.read_sbml_model(draftModelMS)
tensor_network.predict(model)
network.predict(model)
