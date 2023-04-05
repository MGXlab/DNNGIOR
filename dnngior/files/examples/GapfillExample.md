# Gapfill your metabolic model using weights learned with a deep neural network

for this example, two things are needed:

1. **A draft model reconstructed using the ModelSEED pipeline**

Can be obtained either from [ModelSeed's web interface](https://modelseed.org/), from [Kbase](https://www.kbase.us/), or following the instructions in [this notebook](https://github.com/MGXlab/DNNGIOR/blob/main/files/examples/Reconstructing_a_ModelSEED_model_from_a_fasta_genome.md).

In this example, we use a [draft model](https://github.com/MGXlab/DNNGIOR/blob/main/files/models/bh_ungapfilled_model.sbml) that was obtained from an ensembl genome using the last approach. 

2. **A trained neural network**

Can be obtained following the instruction on [this notebook](). 

In this example, we use the previously trained network that is shipped with the package.



We start by gapfilling a model in a complete media:

```python
import dnngior
import cobra

from pathlib import Path
import os

files_path     = Path(os.getcwd()).parents[0]
draftModelPath = os.path.join(files_path, 'models', 'bh_ungapfilled_model.sbml')
trainedNNPath  = os.path.join(files_path, 'NN', 'NN_MS.h5')

gf = dnngior.Gapfill(draftModelPath, trainedNNPath, medium = None, objectiveName = 'bio1')

model_completeMedium = gf.gapfilledModel.copy()
```



The "model_completeMedium" object is a gapfilled cobrapy model. A list with the freshly added reactions is available in the "gf" object.



```python
addedReacts = gf.added_reactions

for reaction in addedReacts:
    print(reaction, "~~", model_completeMedium.reactions.get_by_id(reaction).build_reaction_string(use_metabolite_names = 1), '\n\n')
```


Now, if we have a defined medium where the model should grow. We repeat the the same steps as above, but provide a "medium" argument to the "dnngior.Gapfill" class as the following example:



```python
def_med = {}

#H2O
def_med["EX_cpd00001_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00001_e0":-1}}
#Phosphate
def_med["EX_cpd00009_e0"] = {'lower_bound': -100, 'upper_bound': 100, 'metabolites': {"cpd00009_e0":-1}}
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

draftModelPath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'models', 'bh_ungapfilled_model.sbml')

trainedNNPath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'NN', 'NN_MS.h5')

gf = dnngior.Gapfill(draftModelPath, trainedNNPath, medium = def_med, objectiveName = 'bio1')

#gapfilled model
model_defMedium = gf.gapfilledModel.copy()
model_defMedium.optimize()

```



