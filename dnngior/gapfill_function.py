
import gurobipy as gu
from dnngior.reaction_class import Reaction
import cobra
from copy import deepcopy
import numpy as np


def create_EX_reactions(metab_dict, direction = 'both'):
    '''
    Create an exchange reaction  for external metabolites in a metab_dict.
    
    Direction can be either 'import', 'export', 'both'
    '''
    EX_reactions = {}
    for metab in metab_dict:
        if '_e' in metab:
            EX_reaction_name = 'EX_' + metab[0:8] + '_e0'
            if direction=='import':
                EX_reactions[EX_reaction_name] = {'lower_bound': -1., 'upper_bound': 0., 'metabolites': {metab:-1.}} 
            elif direction=='export':
                EX_reactions[EX_reaction_name] = {'lower_bound': 0, 'upper_bound': 1., 'metabolites': {metab:-1.}}
            else:
                EX_reactions[EX_reaction_name] = {'lower_bound': -1., 'upper_bound': 1., 'metabolites': {metab:-1.}}
        
    return EX_reactions

def build_gurobi_model(reaction_dict, metabolite_dict, objective_name, model_reactions, candidate_reactions, delta):
    """
    Builds a gurobi model with a variable for each reaction in reaction dict,
    constrains variables using the metabolite dict, and sets the objective of the model.
    
    The objective is based on the candidate dict, method, objective name and delta.
    
    Candidate dict includes all reactions that are considered for gap filling, mapping to their cost.
    
    The objective of the model is defined as:
        MAXIMIZE: delta * objective_name (biomass reaction) - sum(cost * variable in candidate reactions)
    
    Parameters
    ----------
    reaction_dict : dict
    reaction_dict formated for gurobi. All reactions (draft + database)
    This is generated by the Reaction class with the .get_gurobi_reaction_dict() function.
    [reaction_id] : (lower_bound, upper_bound)
    
    metabolite_dict :dict
    metabolite_dict formated for gurobi. All metabolites (draft + database)
    This is generated by the Reaction class with the .get_gurobi_metabolite_dict() function.
    [metabolite_id] : {[reactio_id] : stoichiometry}
    
    objective_name : str
    name of the reaction to use as an objective.
    
    model_reactions : set
    reactions that are in the draft model
        
    candidate_reactions : dict
    candidate reactions to add (reactions that are not in the draft model) mapped to a cost
    reactions with low costs are more likely to be included
    [reaction_id] : cost,
    
    
    #export_reactions : 
    #media reactions
    
    delta :float
    parameter of the objective function
    """
    #define the model
    m= gu.Model('mipl')
    #define variables
    for i in reaction_dict:
        b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, lb = reaction_dict[i][0], ub = reaction_dict[i][1])
            
    m.update()
    #set the objective function
    var = m.getVars()
    coef=[]
    
    
    if delta >1:
        for i in var:
            if i.VarName==objective_name:
                coef.append(delta)
                
            elif i.VarName in candidate_reactions:
                cost = candidate_reactions[i.VarName]
                coef.append(-cost) #Costs are assumed to be positive numbers in input.
                
            elif i.VarName in model_reactions: #reactions in the draft model have cost zero
                coef.append(0)
            
        # elif i.VarName in export_reactions: #export reactions have cost zero
        #     coef.append(0)
            
            else:
                print ('Cannot set objective for %s, not objective_name, model, export or candidate reaction!' %i.VarName)
        
    
        

        
    
    else:
        for i in var:
            if i.VarName==objective_name:
                coef.append(1)
                
            elif i.VarName in candidate_reactions:
                cost = candidate_reactions[i.VarName]
                coef.append(0) #Costs are assumed to be positive numbers in input.
                
            elif i.VarName in model_reactions: #reactions in the draft model have cost zero
                coef.append(0)
            
        # elif i.VarName in export_reactions: #export reactions have cost zero
        #     coef.append(0)
            
            else:
                print ('Cannot set objective for %s, not objective_name, model, export or candidate reaction!' %i.VarName)
        
        
        
        
    m.setObjective(gu.LinExpr(coef, var), gu.GRB.MAXIMIZE) #set the objective expression
    m.update()    
    #set the stoichiometric constraints for metabolites    
    for i in metabolite_dict:
        var = metabolite_dict[i].keys()
        var = [m.getVarByName(z) for z in var]
        coef = metabolite_dict[i].values()
        m.addLConstr(gu.LinExpr(coef, var), 'E', 0, i)  
        
    m.update()
    m.setParam('OutputFlag', False) #Keep gurobi silent
    m.optimize()
    return m

def get_minimum(R, R_flux, R_cost, criteria=None):
    
    '''
    Each iteraction of the latendress_gapfill generates a viable gapfilling set.
    Select one based on different criteria: 
        - min_cost: minimun sum of costs in the added reactions
        
        - min_reactions: smallest number of added reactions
        
        - min_flux: smallest sum of fluxes
        
        - else the last iteration is chosen    
    
    Returns
    ------
    The selected reaction set and the value of delta
        
    '''
    
    if criteria == 'min_cost':
        sum_cost = [sum(cost) for cost in R_cost]
        if sum_cost != []:
            min_cost = min(sum_cost)
            min_cost_index = sum_cost.index(min_cost)
            return R[min_cost_index]
    
    if criteria == 'min_reactions':
        len_reactions = [len(e) for e in R]
        if len_reactions != []:
            min_reactions = min(len_reactions)
            min_reactions_index = len_reactions.index(min_reactions)
            return R[min_reactions_index]
    
    if criteria == 'min_flux':
        sum_flux = [sum(flux) for flux in R_flux]
        if sum_flux != []:
            min_flux = min(sum_flux)
            min_flux_index = sum_flux.index(min_flux)
            return R[min_flux_index]
    
    else:
        return R[-1]#, 0
    
def binarySearch(all_reactions_split, N, M, B, result_selection):
    """
    Function based on Latendresse BMC Bioinformatics 2014.
    Input:
        all_reactions_split is Reaction class object containing data for all reactions.
        N is reactions in draft model.
        M is dictionary of all candidate reactions mapping to their cost.
        split_EX_reactions are all environment defining reactions that the model can use, 
            these are not associated to a cost, hence not found in M.
        B is the name of the objective function.
    Output:
        List of the minimum set of candidate reactions to add to gapfill the model.
    """
    R = [] #reaction sets
    R_flux = [] #reaction fluxes
    R_cost = [] #reaction costs
    D = [] #deltas
    proposed_model = []
    
    alpha = 0
    beta = 2 * len(M)
    
    x = [i for i in M if i not in N]
    
    
    #format reactions for gurobi model
    reaction_dict = all_reactions_split.get_gurobi_reaction_dict(all_reactions_split.reactions.keys())
    
    #format metabolites for gurobi_model
    metabolite_dict = all_reactions_split.get_gurobi_metabolite_dict(all_reactions_split.reactions.keys())
    
    #check if the model is gapfillable :)
    gu_model = build_gurobi_model(reaction_dict, metabolite_dict, B, N, M, delta=1)
    
    
    R.append([var.VarName for var in gu_model.getVars() if (var.VarName not in N) and (var.X != 0)])

    max_obj = gu_model.getVarByName(B).X
    print ('Flux through biomass reaction is {:.8f}'.format(max_obj))
    
    
    
    
    #main loop
    
    
    #there is no solution for the model,
    #probably the medium is too restrictive
    if np.round(gu_model.getVarByName(B).X,6) ==0:
        
        return None
    
    
    print ('Flux through biomass reaction is {:.8f}'.format(gu_model.getVarByName(B).X))
    
    
    while abs(alpha - beta) > 1:#converting the "loop until" into a while loop
    
        
        sizeR = len(R[-1])
        
        
        print("current R is: ", sizeR)
        
        delta = int((alpha + beta) / 2.0)
        print ('Delta is {:.2f}'.format(delta))
        
        gu_model = build_gurobi_model(reaction_dict, metabolite_dict, B, N, M, delta)
        
        if np.round(gu_model.getVarByName(B).X,6) > 0:
            print ('Flux through biomass reaction is {:.8f}'.format(gu_model.getVarByName(B).X))
            
            R.append([var.VarName for var in gu_model.getVars() if (var.VarName not in N) and (var.X != 0)])
            #proposed_model.append([var.VarName for var in gu_model.getVars() if var.VarName in N or gu_model.getVarByName(var.VarName).X > 0])
            
            
            # reaction fluxes
            R_flux.append([gu_model.getVarByName(e).X for e in M if np.round(gu_model.getVarByName(e).X,6) > 0])
            
            #costs of candidate reactions retained in the model
            R_cost.append([M[e] for e in M if gu_model.getVarByName(e).X > 0])
            
            #delta
            D.append(delta)
            
            beta = delta
            
            #if len(R[-1])<sizeR:
                
                
            #else:
                
           
        else:
            alpha = delta
            #R[-1] = R[-1]
            #proposed_model[-1] = proposed_model[-2]
            print ('Flux through objective reaction is {:.8f}'.format(gu_model.getVarByName(B).X))
            #alpha = delta
            pass
        print('\n\n', 'condition is currently: ', abs(alpha - beta), '\n\n')
    minimum_set = get_minimum(R, R_flux, R_cost, result_selection) #List that has minimum nr of reactions, sum of cost or sum of flux, dependend on output.
    return  minimum_set

def make_cobra_metabolites(metab_dict):
    '''
    Make a cobra metabolite object for each metabolite in metab_dict.
    Return a dict mapping the metab_id (str) : object.
    '''
    cobra_metabs = {}
    for metab in metab_dict:
        cobra_metabs[metab] = cobra.Metabolite(metab)
        if '_c0' in metab:
            cobra_metabs[metab].compartment = 'c0'
        
        elif '_e0' in metab:
            cobra_metabs[metab].compartment = 'e0'
            
    return cobra_metabs
        
def make_cobra_reaction(reaction_dict, cobra_metabs, rxn):
    '''
    Make a cobra reaction object from a rxn, setting the bounds, and adding metabolites.
    Return the cobra reaction object.
    '''
    reaction = cobra.Reaction(rxn)
    reaction.name = rxn
    reaction.lower_bound = reaction_dict[rxn]['lower_bound']
    reaction.upper_bound = reaction_dict[rxn]['upper_bound']
    metabolites = {cobra_metabs[k]:v for k,v in reaction_dict[rxn]['metabolites'].items()}
    reaction.add_metabolites(metabolites)
    return reaction

def make_cobra_model(reaction_dict, metab_dict, reactions_in_model, objective_name):
    '''
    Make a cobra model of all reactions in reactions_in_model, 
    using a reaction class instance and a metab dict containing all metabs used in the reactions.
    '''
    #Make cobra metabolites and reaction objects.
    cobra_metabs = make_cobra_metabolites(metab_dict)
    cobra_reactions = [make_cobra_reaction(reaction_dict, cobra_metabs, e) for e in reactions_in_model]
    #Make cobra model.
    cobra_model = cobra.Model('tempmodel')
    cobra_model.add_reactions(cobra_reactions)
    cobra_model.objective = objective_name
    return cobra_model

def gapfill(all_reactions, draft_reaction_ids, candidate_reactions, obj_id, default_cost = 1, result_selection = 'min_cost'):
    '''
    Gapfill an incomplete model
    
    Parameters
    ----------
    all_reactions, class_obj
    Reaction class object containing the chemical information of all reactions 
            (bounds, stoichiometry and metabolites).
    
    draft_reaction_ids, set
    Set containing all the reaction ids in the input model for gap filling.
            
    candidate_reactions, dict
    This is a dictionary mapping reaction_ids to their cost during gap filling. 
    When reactions are not present in candidate_reactions, their cost will be default_cost.
         
    obj_id, str
    the reaction id correspoding to the objective.
            
    medium, dict
    The  exchange reactions that the model can use. 
    As default a complete medium is assumed, where the model can use all extracellular metabolites 
    that are used by some reactions present in all_reactions. 
    (This means the metabolites are available extracellularly, not that a transport reaction is added) 
    A custom media can be provided as a dict
    [reaction_id]:{['lower_bound']:-1, ['upper_bound']:1, [metabolites] : [metabolite]: stoichiometry}
            
    default_cost, float
    Reactions in all_reactions that are not in candidate_reactions or in the
    draft model will get this cost for gap filling.
    
    Returns
    -------
    cobra model obj
    the value of the obj
    list of added reactions
    '''
    
    
    
    
    all_reacs = deepcopy(all_reactions.reactions)
    
    
    
    
    
    all_reacs_obj = Reaction()
    all_reacs_obj.reactions = all_reacs.copy()
    cand_reacs = candidate_reactions.copy()
    
    
    
    

    
    
    #Add reactions from all_reactions to candidate_reactions, with cost = default_cost.
    for reaction in all_reacs_obj.reactions:
        if (reaction not in draft_reaction_ids) and (reaction not in cand_reacs):
            cand_reacs[reaction] = default_cost
    
    
    for reaction in candidate_reactions:
        if reaction in draft_reaction_ids:
            del cand_reacs[reaction] #Delete reaction from candidate_reactions if it is present in the starting model.
        
    #Split bidirectional reactions into a forward and reverse reaction.
    all_reactions_split = Reaction()
    
    all_reactions_split.reactions = all_reacs_obj.split_all_bidirectional_reactions(all_reacs_obj.reactions)
    
    #Add reverse reactions to draft_reaction_ids_split and candidate_reactions.
    draft_reaction_ids_split = set()
    
    for reaction in all_reactions_split.reactions:
        forward_version = reaction.replace('_r', '')
        if forward_version in draft_reaction_ids:
            draft_reaction_ids_split.add(reaction)
        else:#If forward version of a reverse reaction is in candidate_reactions.
            if '_r' in reaction:
                cand_reacs[reaction] = cand_reacs[forward_version] #Give reverse reaction same cost as forward version.
    
        
    
    #Run gapfilling algorithm
    split_gapfill_result = binarySearch(all_reactions_split, draft_reaction_ids_split, cand_reacs, obj_id, result_selection)
    
    
    
    if split_gapfill_result is None:
        print("\n\n\n", "media is too restrictive. No growing model can be found :(", "\n\n\n")
        return None, None, None
    
    gapfill_result = set([r.replace('_r', '') for r in split_gapfill_result])
    
    added_reactions = list(gapfill_result) #All reactions that are added to the model during gapfilling.
    
    gapfill_result.update(draft_reaction_ids)
    
    #Create cobra model
    
   
   
    metab_dict      = all_reacs_obj.get_gurobi_metabolite_dict()
    cobra_model     = make_cobra_model(all_reacs, metab_dict, gapfill_result, obj_id) 
    objective_value = cobra_model.optimize().objective_value    
    print ('Objective value is %f.' %objective_value)

   
    return cobra_model, objective_value, added_reactions#, gapfill_result, all_reacs_obj
