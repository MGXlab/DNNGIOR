Load dependencies


```python
import cobra
#Avoid the warnings about the sbml formating
import logging
logging.getLogger("cobra").setLevel(logging.ERROR)

import os
import sys
from pathlib import Path
path = Path.cwd()
sys.path.append(path)

import numpy as np
```


```python
from reaction_class import Reaction
import gapfill_function
from build_model import *
```

### Cobrapy crash tutorial for browsing GSMMs

I find cobrapy one the best ways to interact with GSMMs in a programatic way. It's quite well documented (https://cobrapy.readthedocs.io/en/latest/)

I am only going to illustrate some basic functionalities:

1) load a model from an XML file


```python
model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-draft.xml')
gpseq_model = cobra.io.read_sbml_model(model_location)
```

As expected from GSMMs, the model obj has genes, reactions, metabolites, medium, demand reactions, and compartments.

They are all iterable python objects that have additional information as id, name, charge, formula, etc


```python

print('number of reactions : ', len(gpseq_model.reactions))


print('number of metabolites : ', len(gpseq_model.metabolites))


print('number of genes : ', len(gpseq_model.genes))


print('Compartments : ', gpseq_model.compartments)

#1 reaction info
print('id : ', gpseq_model.reactions.rxn00004_c0.id)
print('name : ', gpseq_model.reactions.rxn00004_c0.name)

print(gpseq_model.reactions.rxn00004_c0.reaction)
print(gpseq_model.reactions.rxn00004_c0.build_reaction_string(use_metabolite_names=1))

print('\n\n', 'reaction genes : ',  gpseq_model.reactions.rxn00004_c0.genes , '\n\n')
print('\n\n', 'reaction gene rules : ', gpseq_model.reactions.rxn00004_c0.gene_reaction_rule, '\n\n',)


#1 metabolite info
print('id : ', gpseq_model.metabolites.cpd00020_c0.id)
print('name : ', gpseq_model.metabolites.cpd00020_c0.name)
print('formula : ', gpseq_model.metabolites.cpd00020_c0.formula)
print('charge : ', gpseq_model.metabolites.cpd00020_c0.charge)
print('elements : ', gpseq_model.metabolites.cpd00020_c0.elements)
print('weight : ', gpseq_model.metabolites.cpd00020_c0.formula_weight)
print('reactions that use this metabolite : ', [i.id for i in gpseq_model.metabolites.cpd00020_c0.reactions])

```

    number of reactions :  3846
    number of metabolites :  3113
    number of genes :  17464
    Compartments :  {'c0': '', 'e0': '', 'p0': ''}
    id :  rxn00004_c0
    name :  4-hydroxy-4-methyl-2-oxoglutarate pyruvate-lyase (pyruvate-forming)
    cpd02570_c0 <=> 2.0 cpd00020_c0
    Parapyruvate-c0 <=> 2.0 Pyruvate-c0
    
    
     reaction genes :  frozenset({<Gene gp_NODE___THREE____EIGHT____ZERO___length___FIVE____FIVE____SIX____SEVEN____SIX___cov___SIX____THREE____DOT____TWO____EIGHT____NINE____TWO____NINE____SEVEN_____THREE____ONE____THREE____SEVEN____ONE_____THREE____ZERO____SIX____SEVEN____NINE__ at 0x202df870610>, <Gene gp_NODE___THREE____EIGHT____THREE___length___FIVE____FIVE____ZERO____FOUR____TWO___cov___FOUR____ZERO____DOT____ONE____THREE____FIVE____FIVE____SEVEN____EIGHT_____FOUR____SEVEN____NINE____SIX____TWO_____FOUR____EIGHT____SIX____ONE____FIVE__ at 0x202df870640>, <Gene gp_NODE___SEVEN____FIVE____SIX___length___THREE____THREE____TWO____FOUR____SEVEN___cov___ONE____EIGHT____DOT____ONE____FIVE____FOUR____NINE____FOUR____SEVEN_____THREE____ONE____ZERO____TWO____SEVEN_____THREE____ZERO____THREE____SIX____TWO__ at 0x202df870670>, <Gene gp_NODE___EIGHT____ONE____NINE____THREE___length___FOUR____ZERO____SIX____FIVE___cov___ONE____SIX____DOT____ONE____TWO____ONE____SIX____NINE____SIX_____THREE____SEVEN____EIGHT____FIVE_____THREE____THREE____ZERO____THREE__ at 0x202df8706a0>, <Gene gp_NODE___EIGHT____FOUR___length___ONE____FOUR____THREE____ZERO____TWO____NINE___cov___SIX____FIVE____DOT____FIVE____TWO____SIX____NINE____ONE____FOUR_____TWO____TWO____SEVEN____SIX____EIGHT_____TWO____THREE____FOUR____THREE____SIX__ at 0x202df8706d0>, <Gene gp_NODE___ONE____THREE____SIX___length___ONE____ZERO____FOUR____SEVEN____EIGHT____FOUR___cov___THREE____ONE____DOT____SEVEN____TWO____EIGHT____SIX____EIGHT____ONE_____ONE____EIGHT____EIGHT____THREE____FOUR_____ONE____EIGHT____THREE____SEVEN____THREE__ at 0x202df8704f0>, <Gene gp_NODE___NINE____NINE____ZERO___length___TWO____SIX____EIGHT____THREE____FIVE___cov___ONE____SIX____DOT____SEVEN____FIVE____ONE____SEVEN____ONE____EIGHT_____ONE____THREE____NINE____NINE____TWO_____ONE____FOUR____SEVEN____ZERO____FIVE__ at 0x202df870700>, <Gene gp_NODE___ONE____NINE____SEVEN___length___EIGHT____ONE____SIX____ZERO____FIVE___cov___EIGHT____NINE____DOT____EIGHT____NINE____ONE____THREE____SIX____SEVEN_____FIVE____TWO____THREE____EIGHT____NINE_____FIVE____THREE____ZERO____FIVE____SEVEN__ at 0x202df870520>, <Gene gp_NODE___TWO____ONE____SIX___length___SEVEN____SEVEN____FOUR____THREE____ZERO___cov___THREE____ZERO____DOT____FIVE____SIX____SIX____FIVE____SIX____FIVE_____ONE____FOUR____TWO____TWO____FIVE_____ONE____THREE____SEVEN____FOUR____ZERO__ at 0x202df870550>, <Gene gp_NODE___TWO____FIVE____THREE____SIX___length___ONE____ONE____FIVE____SIX____THREE___cov___TWO____TWO____DOT____EIGHT____THREE____ZERO____NINE____EIGHT____SEVEN_____TWO____ONE____EIGHT____ZERO_____TWO____SIX____SIX____FIVE__ at 0x202df870580>, <Gene gp_NODE___TWO____EIGHT____ONE____ONE___length___ONE____ZERO____FIVE____FIVE____EIGHT___cov___ONE____THREE____DOT____TWO____FOUR____ONE____SIX____FOUR____FIVE_____FOUR____ONE____SIX_____ONE____ZERO____EIGHT____FOUR__ at 0x202df8705b0>, <Gene gp_NODE___THREE____ZERO____FIVE____TWO___length___NINE____EIGHT____TWO____TWO___cov___ONE____THREE____DOT____FIVE____NINE____ONE____NINE____NINE____THREE_____EIGHT____THREE____EIGHT____ONE_____SEVEN____EIGHT____NINE____SIX__ at 0x202df8705e0>}) 
    
    
    
    
     reaction gene rules :  gp_NODE___ONE____THREE____SIX___length___ONE____ZERO____FOUR____SEVEN____EIGHT____FOUR___cov___THREE____ONE____DOT____SEVEN____TWO____EIGHT____SIX____EIGHT____ONE_____ONE____EIGHT____EIGHT____THREE____FOUR_____ONE____EIGHT____THREE____SEVEN____THREE__ or gp_NODE___ONE____NINE____SEVEN___length___EIGHT____ONE____SIX____ZERO____FIVE___cov___EIGHT____NINE____DOT____EIGHT____NINE____ONE____THREE____SIX____SEVEN_____FIVE____TWO____THREE____EIGHT____NINE_____FIVE____THREE____ZERO____FIVE____SEVEN__ or gp_NODE___TWO____ONE____SIX___length___SEVEN____SEVEN____FOUR____THREE____ZERO___cov___THREE____ZERO____DOT____FIVE____SIX____SIX____FIVE____SIX____FIVE_____ONE____FOUR____TWO____TWO____FIVE_____ONE____THREE____SEVEN____FOUR____ZERO__ or gp_NODE___TWO____FIVE____THREE____SIX___length___ONE____ONE____FIVE____SIX____THREE___cov___TWO____TWO____DOT____EIGHT____THREE____ZERO____NINE____EIGHT____SEVEN_____TWO____ONE____EIGHT____ZERO_____TWO____SIX____SIX____FIVE__ or gp_NODE___TWO____EIGHT____ONE____ONE___length___ONE____ZERO____FIVE____FIVE____EIGHT___cov___ONE____THREE____DOT____TWO____FOUR____ONE____SIX____FOUR____FIVE_____FOUR____ONE____SIX_____ONE____ZERO____EIGHT____FOUR__ or gp_NODE___THREE____ZERO____FIVE____TWO___length___NINE____EIGHT____TWO____TWO___cov___ONE____THREE____DOT____FIVE____NINE____ONE____NINE____NINE____THREE_____EIGHT____THREE____EIGHT____ONE_____SEVEN____EIGHT____NINE____SIX__ or gp_NODE___THREE____EIGHT____ZERO___length___FIVE____FIVE____SIX____SEVEN____SIX___cov___SIX____THREE____DOT____TWO____EIGHT____NINE____TWO____NINE____SEVEN_____THREE____ONE____THREE____SEVEN____ONE_____THREE____ZERO____SIX____SEVEN____NINE__ or gp_NODE___THREE____EIGHT____THREE___length___FIVE____FIVE____ZERO____FOUR____TWO___cov___FOUR____ZERO____DOT____ONE____THREE____FIVE____FIVE____SEVEN____EIGHT_____FOUR____SEVEN____NINE____SIX____TWO_____FOUR____EIGHT____SIX____ONE____FIVE__ or gp_NODE___SEVEN____FIVE____SIX___length___THREE____THREE____TWO____FOUR____SEVEN___cov___ONE____EIGHT____DOT____ONE____FIVE____FOUR____NINE____FOUR____SEVEN_____THREE____ONE____ZERO____TWO____SEVEN_____THREE____ZERO____THREE____SIX____TWO__ or gp_NODE___EIGHT____ONE____NINE____THREE___length___FOUR____ZERO____SIX____FIVE___cov___ONE____SIX____DOT____ONE____TWO____ONE____SIX____NINE____SIX_____THREE____SEVEN____EIGHT____FIVE_____THREE____THREE____ZERO____THREE__ or gp_NODE___EIGHT____FOUR___length___ONE____FOUR____THREE____ZERO____TWO____NINE___cov___SIX____FIVE____DOT____FIVE____TWO____SIX____NINE____ONE____FOUR_____TWO____TWO____SEVEN____SIX____EIGHT_____TWO____THREE____FOUR____THREE____SIX__ or gp_NODE___NINE____NINE____ZERO___length___TWO____SIX____EIGHT____THREE____FIVE___cov___ONE____SIX____DOT____SEVEN____FIVE____ONE____SEVEN____ONE____EIGHT_____ONE____THREE____NINE____NINE____TWO_____ONE____FOUR____SEVEN____ZERO____FIVE__ 
    
    
    id :  cpd00020_c0
    name :  Pyruvate-c0
    formula :  C3H3O3
    charge :  -1
    elements :  {'C': 3, 'H': 3, 'O': 3}
    weight :  87.05412
    reactions that use this metabolite :  ['rxn02900_c0', 'rxn90067_c0', 'rxn00145_c0', 'rxn00966_c0', 'rxn08043_c0', 'rxn03610_c0', 'rxn00422_c0', 'rxn12707_c0', 'rxn01667_c0', 'rxn00168_c0', 'rxn00517_c0', 'rxn02346_c0', 'rxn00727_c0', 'rxn05485_c0', 'rxn00165_c0', 'rxn15021_c0', 'rxn02230_c0', 'rxn00148_c0', 'rxn00162_c0', 'rxn15266_c0', 'rxn00146_c0', 'rxn00149_c0', 'rxn14068_c0', 'rxn03909_c0', 'rxn00255_c0', 'rxn00161_c0', 'rxn01203_c0', 'rxn00289_c0', 'rxn00160_c0', 'rxn05469_c0', 'rxn07450_c0', 'rxn00278_c0', 'rxn00177_c0', 'rxn00166_c0', 'rxn02210_c0', 'rxn05560_c0', 'rxn14245_c0', 'rxn00157_c0', 'rxn01644_c0', 'rxn04675_c0', 'rxn01322_c0', 'rxn08535_c0', 'rxn03884_c0', 'rxn10951_c0', 'rxn04092_c0', 'rxn00801_c0', 'rxn00258_c0', 'rxn00304_c0', 'rxn00566_c0', 'rxn05938_c0', 'rxn00783_c0', 'rxn09167_c0', 'rxn05647_c0', 'rxn00164_c0', 'rxn01354_c0', 'rxn00460_c0', 'rxn00499_c0', 'rxn05617_c0', 'rxn15166_c0', 'rxn00678_c0', 'rxn05607_c0', 'rxn00147_c0', 'rxn00656_c0', 'rxn06282_c0', 'rxn15491_c0', 'rxn00272_c0', 'rxn09402_c0', 'rxn05109_c0', 'rxn05988_c0', 'rxn00004_c0', 'rxn00411_c0', 'rxn00159_c0', 'rxn00491_c0', 'rxn12670_c0', 'rxn16114_c0', 'rxn00726_c0', 'rxn00840_c0', 'rxn11731_c0', 'rxn02185_c0', 'rxn05226_c0', 'rxn03841_c0', 'rxn15044_c0', 'rxn15267_c0', 'rxn00950_c0', 'rxn00328_c0', 'rxn01735_c0', 'rxn01990_c0', 'rxn03379_c0', 'rxn00003_c0', 'rxn00154_c0', 'rxn00904_c0', 'rxn00011_c0', 'rxn00473_c0', 'rxn01242_c0', 'rxn14077_c0', 'rxn00191_c0', 'rxn00250_c0', 'rxn01202_c0', 'rxn11703_c0', 'rxn00540_c0', 'rxn03924_c0', 'rxn00500_c0', 'rxn01619_c0', 'rxn06404_c0', 'rxn01307_c0', 'rxn02925_c0', 'rxn00987_c0']
    

All these fields are editable. 

Additionally, there is a suite of functions (check the docs).

To run FBA (the draft model is not functional)


```python
solution = gpseq_model.optimize()
print(solution.objective_value)

print(solution.fluxes)
```

    0.0
    rxn00001_c0       0.0
    rxn00002_c0       0.0
    rxn00003_c0       0.0
    rxn00004_c0       0.0
    rxn00006_c0       0.0
                     ... 
    EX_cpd03662_e0    0.0
    EX_cpd00178_e0    0.0
    EX_cpd01618_e0    0.0
    bio1              0.0
    DM_cpd01042_c0    0.0
    Name: fluxes, Length: 3846, dtype: float64
    

As a quick illustration. I am going to create a reaction for ATP hydrolysis, add it to the model and use this as an objective instead of the biomass.


```python
reaction = cobra.Reaction('ATHyd')
reaction.name = 'ATP hydrolysis'
atp = gpseq_model.metabolites.cpd00002_c0
water = gpseq_model.metabolites.cpd00001_c0
adp = gpseq_model.metabolites.cpd00008_c0
pi = gpseq_model.metabolites.cpd00009_c0
proton = gpseq_model.metabolites.cpd00067_c0

reaction.add_metabolites({atp:-1, water:-1, adp:1, pi:1, proton:1})

reaction.lower_bound=0
reaction.upper_bound=1000

gpseq_model.add_reactions([reaction])


gpseq_model.reactions.ATHyd.objective_coefficient=1
gpseq_model.reactions.bio1.objective_coefficient=0
```


```python
solution = gpseq_model.optimize()
print(solution.objective_value)

print(solution.fluxes)
```

    0.0
    rxn00001_c0       0.0
    rxn00002_c0       0.0
    rxn00003_c0       0.0
    rxn00004_c0       0.0
    rxn00006_c0       0.0
                     ... 
    EX_cpd00178_e0    0.0
    EX_cpd01618_e0    0.0
    bio1              0.0
    DM_cpd01042_c0    0.0
    ATHyd             0.0
    Name: fluxes, Length: 3847, dtype: float64
    

Finally, to save the modified model as an XML file


```python
new_model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-draft_ATP_Obj.xml')
cobra.io.write_sbml_model(gpseq_model, new_model_location)
```

### Gapfilling a model on complete media with equally penalized database reactions

Use the function gapfill_function.gapfill

as inputs use:

1) The reaction ids of your draft model as a set that we call 'draft_reaction_ids'




```python
model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-draft_ATP_Obj.xml')
draft_model = Reaction(model=model_location)
draft_model_cobra_obj = cobra.io.read_sbml_model(model_location) #this will later be used to retrieve the genes of the draft model
draft_reaction_ids = set(draft_model.reactions) 
```

2) A Reaction class object that contains all the reactions of the ModelSeed database, called 'all_reactions' 


```python
# Load the database reactions.tsv file
biochem = os.path.join(path.parent, 'files',  'biochemistry', 'reactions.tsv')
db_reactions = Reaction(biochem_input=biochem)

#Combine all reactions into one object

all_reactions = Reaction()
all_reactions.reactions = all_reactions.add_dict(draft_model.reactions, db_reactions.reactions) 
```

3) the id of an objective function, in our case 'bio1' or 'ATHyd'

4) A method to select the best gapfilled sets. In this case we use the minimum number of added reactions 'min_reactions'


```python
model, obj, new_reacs = gapfill_function.gapfill(all_reactions, draft_reaction_ids, {}, 'bio1', medium = 'complete', result_selection = 'min_reactions')

```

    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Delta is 70671.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 35335.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 17667.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 8833.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 4416.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 2208.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 1104.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 552.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 276.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 138.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 69.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 34.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 17.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 8.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 0.60693822
    Delta is 4.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through objective reaction is 0.00000000
    Delta is 6.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through objective reaction is 0.00000000
    Delta is 7.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 0.36469980
    Objective value is 1.000000.
    

The added reactions as listed in the 'new_reacs' list and a cobrapy model object is generated called 'model'


```python
print('number of added reactions : ', len(new_reacs))
print('number of reactions in model : ', len(model.reactions))
print('number of metabolites in model : ', len(model.metabolites))

#run FBA
solution = model.optimize()
print('objective value :' , solution.objective_value)
```

    number of added reactions :  79
    number of reactions in model :  3742
    number of metabolites in model :  3130
    objective value : 1.0
    

The returned model only has reaction and metabolite ids and the upper and lower bounds are fixed between -1 and 1.

To obtain a model with all the annotations from the modelSEED database and internal fluxes constrained between -1000 and 1000, plus the gene annotations from  the draft model run:


```python
refined_model = refine_model(model, draft_model_cobra_obj)
```

    Read LP format model from file C:\Users\u0139894\AppData\Local\Temp\tmpsosy0t93.lp
    Reading time = 0.04 seconds
    : 3130 rows, 7484 columns, 35582 nonzeros
    


```python
solution = refined_model.optimize()
print(solution.objective_value) #the exchange reactions are still bounded between -1 and 1
```

    1.0
    

save the model


```python
gapfilled_model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-full_media_gpfilled.xml')
cobra.io.write_sbml_model(refined_model, gapfilled_model_location)
```

### Gapfilling a model on a defined media

We need the same inputs as above, but now we also need a media




```python
media_file = os.path.join(path.parent, 'files', 'biochemistry', 'Nitrogen-Nitrite_media.tsv')

Nit_media = {}

with open(media_file) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        Nit_media['EX_' + a[0] + '_e0'] = {'lower_bound':-1, 'upper_bound':1, 'metabolites':{a[0]+'_e0':-1.0}}
```


```python
model2, obj2, new_reacs2 = gapfill_function.gapfill(all_reactions, draft_reaction_ids, {}, 'ATHyd', result_selection = 'min_reactions', medium=Nit_media)

```

    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Delta is 65919.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 32959.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 16479.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 8239.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 4119.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 2059.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 1029.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 514.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 257.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 128.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 64.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 32.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 16.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 8.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 4.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 2.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 1.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Objective value is 1.000000.
    


```python
print('number of added reactions : ', len(new_reacs2))
solution = model2.optimize()
print('objective value : ', solution.objective_value)
```

    number of added reactions :  4
    objective value :  1.0
    


```python
refined_model2 = refine_model(model2, draft_model_cobra_obj)
gapfilled_model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-defined_media_gpfilled.xml')
cobra.io.write_sbml_model(refined_model2, gapfilled_model_location)

```

    Read LP format model from file C:\Users\u0139894\AppData\Local\Temp\tmp_0yv8p3p.lp
    Reading time = 0.04 seconds
    : 3113 rows, 7334 columns, 34746 nonzeros
    

### Gapfilig a model with different costs for reactions

Now we need to provide a python dictionary with reaction costs. The function will automatically give a cost of zero to reactions that are already in the draft model and also to the exchange reactions of the defined media. It will give a default cost (set to 1) to reactions that are not in the draft model nor in the dictionary with costs.


```python
#select random reaction and add random costs between 0 and 1
candidate_reactions = {}
for i in  all_reactions.reactions:
    if np.random.uniform()<0.02:
        candidate_reactions[i] = np.random.uniform()
```

We will also use a different criteria to choose the gapfilled reaction. Instead of 'min_reactions' we use 'min_cost', selecting the set with a minimum sum of costs.


```python
model3, obj3, new_reacs3 = gapfill_function.gapfill(all_reactions, draft_reaction_ids, candidate_reactions, 'bio1', result_selection = 'min_costs', medium=Nit_media)

```

    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Delta is 65919.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 32959.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 16479.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 8239.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 4119.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 2059.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 1029.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 514.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 257.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 128.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 64.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 32.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 16.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 8.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 1.00000000
    Delta is 4.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 0.08789180
    Delta is 2.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through objective reaction is 0.04747566
    Delta is 3.00
    Warning for adding constraints: zero or small (< 1e-13) coefficients, ignored
    Flux through biomass reaction is 0.08339543
    Objective value is 0.083395.
    


```python
print('number of added reactions : ', len(new_reacs3))
solution = model3.optimize()
print(solution.objective_value)
```

    number of added reactions :  37
    0.08339543116137185
    


```python
refined_model3 = refine_model(model3, draft_model_cobra_obj)
gapfilled_model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-random_weights_gpfilled.xml')
cobra.io.write_sbml_model(refined_model3, gapfilled_model_location)
```

    Read LP format model from file C:\Users\u0139894\AppData\Local\Temp\tmpwl6otvlt.lp
    Reading time = 0.04 seconds
    : 3116 rows, 7400 columns, 34956 nonzeros
    

### Quick idea of what to do when no flux is observed on a particular media


1) Run a gapfill on the with complete media option (as in [5])

2) Run flux variability analysis (for big models this takes a long time to run):


```python
model_exchanges = [i for i in model.reactions if ('EX_' in i.id) and ('_e0' in i.id)]
```


```python
fva = cobra.flux_analysis.flux_variability_analysis(model, model_exchanges)
print(fva)
```

The exchange reaction where the minimum and maximum are negative are likely required in the medium. They are essential for the objective and likely the database reactions have no way to make them, so gapdfilling does not wotrk.
