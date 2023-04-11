# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 22:42:00 2023

@author: danie
"""
from pathlib import Path
import os
import gurobipy as gu
import cobra
from copy import deepcopy
import numpy as np

from dnngior.variables import *
from dnngior.reaction_class import Reaction
from dnngior.NN_Predictor import NN
from dnngior import build_model


class Gapfill:
    def __init__(self, 
                draftModel, 
                trainedNNPath = None, 
                medium        = None, 
                objectiveName = 'bio1', 
                dbType        = 'ModelSEED'):

        self.objectiveName    = objectiveName
        self.trainedNNPath    = trainedNNPath
        self.dbType           = dbType
        self.draftModel       = cobra.io.read_sbml_model(draftModel)
        self.draft_reaction   = Reaction( model = draftModel )
        self.medium           = medium
        self.result_selection = "min_reactions"

        if dbType == "ModelSEED":
            self.path_to_biochem  = MODELSEED_REACTIONS
            if trainedNNPath is None: 
                self.trainedNNPath = TRAINED_NN_MSEED
        else:
            if trainedNNPath is None:
                self.trainedNNPath = TRAINED_NN_BIGG
            return "dbType %s is currently not supported" % dbType


        # Build a Reaction object for the exchange reactions; 
        # if you have a defined medium, set the fixed_bounds argument accordingly
        self.exchange_reacs = Reaction( model = os.path.join(MODELS_PATH, 'exchangeReactions.sbml'),
                                       fixed_bounds = self.medium)
        
        self.db_reactions           = Reaction(biochem_input = self.path_to_biochem)
        self.db_reactions.reactions = self.db_reactions.add_dict(self.exchange_reacs.reactions, self.db_reactions.reactions)
        
        self.all_reactions           = Reaction(fixed_bounds = self.medium) 
        self.all_reactions.reactions = self.all_reactions.add_dict(self.draft_reaction.reactions, self.db_reactions.reactions)
        
        self.draft_reaction_ids = set(self.draft_reaction.reactions)
        
        if self.medium is not None:
            #####remove the predefined exchange reactions#####
            for react in self.exchange_reacs.reactions:
                if react in self.draft_reaction_ids:
                    self.draft_reaction_ids.remove(react)
            ##################################################
            
            
            #######apply the defined media condition############
            for react in self.exchange_reacs.reactions:
                if react not in self.medium:
                    self.all_reactions.reactions[react]["lower_bound"] = 0
            ####################################################
        
        self.weights = {}
        
        if self.trainedNNPath is not None:
            
            # Predict weights
            p = NN(path = self.trainedNNPath).predict( self.draft_reaction_ids) 
            
            for i in p:
                self.weights[i]  = np.round(1-p[i], 10)

            model_NN_gf = self.gapfill(self.all_reactions,
                                       self.draft_reaction_ids,
                                       self.weights,
                                       self.objectiveName,
                                       self.result_selection
                                       )

            if self.medium is not None:
                self.gapfilledModel = build_model.refine_model(model_NN_gf, self.draftModel, unscalled = list(self.medium.keys()))

            else:

                if self.dbType=='BiGG':
                    self.gapfilledModel = model_NN_gf

                else:
                    self.gapfilledModel = build_model.refine_model(model_NN_gf, self.draftModel)
                
    def gapfill(self, 
                all_reactions, 
                draft_reaction_ids, 
                candidate_reactions, 
                obj_id, 
                result_selection,
                default_cost = 1, 
                ):
        
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
        
        # Add reactions from all_reactions to candidate_reactions, with cost = default_cost.
        for reaction in all_reacs_obj.reactions:
            if (reaction not in draft_reaction_ids) and (reaction not in cand_reacs):
                cand_reacs[reaction] = default_cost
        
        # Delete reaction from candidate_reactions if it is present in the starting model.
        for reaction in candidate_reactions:
            if reaction in draft_reaction_ids:
                del cand_reacs[reaction]    
        
        #Split bidirectional reactions into a forward and reverse reaction.
        all_reactions_split = Reaction()
        
        all_reactions_split.reactions = all_reacs_obj.split_all_bidirectional_reactions(all_reacs_obj.reactions)
        
        # Add reverse reactions to draft_reaction_ids_split and candidate_reactions.
        draft_reaction_ids_split = set()
        
        for reaction in all_reactions_split.reactions:
            forward_version = reaction.replace('_r', '')
            if forward_version in draft_reaction_ids:
                draft_reaction_ids_split.add(reaction)
            
            # If forward version of a reverse reaction is in candidate_reactions.
            else:
                if '_r' in reaction:
                    # Give reverse reaction same cost as forward version.
                    cand_reacs[reaction] = cand_reacs[forward_version] 


        # Run gapfilling algorithm
        split_gapfill_result = self.binarySearch(all_reactions_split, 
                                                 draft_reaction_ids_split, 
                                                 cand_reacs, 
                                                 obj_id, 
                                                 )
        
        
        # If the database and media did not result in a functional model
        if split_gapfill_result is None:
            print("\n\n\n", "media is too restrictive. No growing model can be found :(", "\n\n\n")
            return None, None, None
        
        gapfill_result = set([r.replace('_r', '') for r in split_gapfill_result])
        
        self.added_reactions = list(gapfill_result) #All reactions that are added to the model during gapfilling.
        
        gapfill_result.update(draft_reaction_ids)
        
        # Create cobra model
        metab_dict      = all_reacs_obj.get_gurobi_metabolite_dict()
        cobra_model     = self.make_cobra_model(all_reacs, 
                                                metab_dict, 
                                                gapfill_result, 
                                                obj_id) 
        self.objective_value = cobra_model.optimize().objective_value    
        print ('Objective value is %f.' %self.objective_value)
       
        return cobra_model
    
    def binarySearch(self, 
                     all_reactions_split, 
                     N, 
                     M, 
                     B
                     ):
        """
        Function modified from Latendresse BMC Bioinformatics 2014.
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
        R      = [] #reaction sets
        R_flux = [] #reaction fluxes
        R_cost = [] #reaction costs
        D      = [] #deltas
        proposed_model = []
        
        alpha = 0
        beta = 2 * len(M)
        
        x = [i for i in M if i not in N]
        
        
        # Format reactions for gurobi model
        reaction_dict = all_reactions_split.get_gurobi_reaction_dict(all_reactions_split.reactions.keys())
        
        # Format metabolites for gurobi_model
        metabolite_dict = all_reactions_split.get_gurobi_metabolite_dict(all_reactions_split.reactions.keys())
        
        # Check if the model is gapfillable :)
        gu_model = self.build_gurobi_model(reaction_dict, metabolite_dict, B, N, M, delta = 1)
        
        
        R.append([var.VarName for var in gu_model.getVars() if (var.VarName not in N) and (var.X != 0)])
    
        max_obj = gu_model.getVarByName(B).X
        print ('Flux through biomass reaction is {:.8f}'.format(max_obj))
                
        
        #there is no solution for the model,
        #probably the medium is too restrictive
        if np.round(gu_model.getVarByName(B).X,6) ==0:
            
            return None
        
        print ('Flux through biomass reaction is {:.8f}'.format(gu_model.getVarByName(B).X))
        
        
        while abs(alpha - beta) > 1:
            
            sizeR = len(R[-1])
            
            
            print("current R is: ", sizeR)
            
            delta = int((alpha + beta) / 2.0)
            print ('Delta is {:.2f}'.format(delta))
            
            gu_model = self.build_gurobi_model(reaction_dict, metabolite_dict, B, N, M, delta)
            
            if np.round(gu_model.getVarByName(B).X,6) > 0:
                print ('Flux through biomass reaction is {:.8f}'.format(gu_model.getVarByName(B).X))
                
                R.append([var.VarName for var in gu_model.getVars() if (var.VarName not in N) and (var.X != 0)])
                
                # Reaction fluxes
                R_flux.append([gu_model.getVarByName(e).X for e in M if np.round(gu_model.getVarByName(e).X,6) > 0])
                
                # Costs of candidate reactions retained in the model
                R_cost.append([M[e] for e in M if gu_model.getVarByName(e).X > 0])
                
                # Delta
                D.append(delta)
                
                beta = delta
               
            else:
                alpha = delta
                #R[-1] = R[-1]
                #proposed_model[-1] = proposed_model[-2]
                print ('Flux through objective reaction is {:.8f}'.format(gu_model.getVarByName(B).X))
                #alpha = delta
                pass
            print('\n\n', 'condition is currently: ', abs(alpha - beta), '\n\n')
        
        # List that has minimum nr of reactions, sum of cost or sum of flux, dependend on output.
        minimum_set = self.get_minimum(R, R_flux, R_cost, criteria = self.result_selection) 

        return  minimum_set
    
    def build_gurobi_model(self, 
                           reaction_dict, 
                           metabolite_dict, 
                           objective_name, 
                           model_reactions, 
                           candidate_reactions, 
                           delta):
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
        
        m = gu.Model('mipl')

        # Define variables
        for i in reaction_dict:
            b = m.addVar(vtype = gu.GRB.CONTINUOUS, name = i, lb = reaction_dict[i][0], ub = reaction_dict[i][1])

        m.update()

        # Set the objective function
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
                
                else:
                    print ('Cannot set objective for %s, not objective_name, model, export or candidate reaction!' %i.VarName)
        
        else:
            for i in var:
                if i.VarName == objective_name:
                    coef.append(1)

                elif i.VarName in candidate_reactions:
                    cost = candidate_reactions[i.VarName]
                    coef.append(0) #Costs are assumed to be positive numbers in input.

                elif i.VarName in model_reactions: #reactions in the draft model have cost zero
                    coef.append(0)

                else:
                    print ('Cannot set objective for %s, not objective_name, model, export or candidate reaction!' %i.VarName)
            
 
        # Set the objective expression
        m.setObjective(gu.LinExpr(coef, var), gu.GRB.MAXIMIZE) 
        m.update()    
        
        # Set the stoichiometric constraints for metabolites    
        for i in metabolite_dict:
            var = metabolite_dict[i].keys()
            var = [m.getVarByName(z) for z in var]
            coef = metabolite_dict[i].values()
            m.addLConstr(gu.LinExpr(coef, var), 'E', 0, i)  
            
        m.update()
        # Keep gurobi silent
        m.setParam('OutputFlag', False) 
        m.optimize()
 
        return m
    
    def get_minimum(self, 
                    R, 
                    R_flux, 
                    R_cost, 
                    criteria=None):
        
        '''
        Each iteraction of the binarySearch generates a viable gapfilling set.
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
            return R[-1]
    
    def make_cobra_metabolites(self, 
                               metab_dict):
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
    
    def make_cobra_reaction(self, 
                            reaction_dict, 
                            cobra_metabs, 
                            rxn):
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
    
    def make_cobra_model(self, 
                         reaction_dict, 
                         metab_dict, 
                         reactions_in_model, 
                         objective_name):
        '''
        Make a cobra model of all reactions in reactions_in_model, 
        using a reaction class instance and a metab dict containing all metabs used in the reactions.
        '''
        #Make cobra metabolites and reaction objects.
        cobra_metabs = self.make_cobra_metabolites(metab_dict)
        cobra_reactions = [self.make_cobra_reaction(reaction_dict, cobra_metabs, e) for e in reactions_in_model]
        #Make cobra model.
        cobra_model = cobra.Model('tempmodel')
        cobra_model.add_reactions(cobra_reactions)
        cobra_model.objective = objective_name
        return cobra_model

