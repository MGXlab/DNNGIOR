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
from pandas import read_csv

from dnngior.variables import *
from dnngior.reaction_class import Reaction
from dnngior.NN_Predictor import NN
from dnngior import build_model


class Gapfill:
    def __init__(self,
                draftModel,                     #Model to be gap-filled, required
                trainedNNPath = None,           #Path to Neural network, if None is provided will default to dbType
                medium        = None,           #User defined medium
                medium_file   = None,           #Path to used defined medium file
                black_list    = None,           #reactions to remove from candidates
                grey_list     = None,           #reactions to punish with high cost
                punish_cost   = 1000.0,         #high cost to punish with
                default_cost  = 1.0,            #default_cost of non-predicted_reactions
                objectiveName = 'bio1',         #Name of objective function
                dbType        = 'ModelSEED',    #Type of database for reactions
                gapfill       = True):          #Boolean of whether to auto-continue gapfilling

        self.objectiveName    = objectiveName
        self.trainedNNPath    = trainedNNPath
        self.dbType           = dbType
        self.draftModel       = cobra.io.read_sbml_model(draftModel)
        self.draft_reaction   = Reaction( model = draftModel )
        self.black_list       = black_list
        self.grey_list        = grey_list
        self.punish_cost      = punish_cost
        self.medium           = medium
        self.medium_file      = medium_file
        self.default_cost     = default_cost

        print("Gap-filling database = ", self.dbType)

        if dbType == "ModelSEED":
            self.path_to_biochem   = MODELSEED_REACTIONS
            self.path_to_exchanges = MODELSEED_EXCHANGES
            if trainedNNPath is None:
                self.trainedNNPath = TRAINED_NN_MSEED
        elif dbType == "BiGG":
            self.path_to_biochem = BIGG_REACTIONS
            self.path_to_exchanges = BIGG_EXCHANGES
            if trainedNNPath is None:
                self.trainedNNPath = TRAINED_NN_BIGG
        else:
            return "dbType %s is not supported" % dbType

        self.draft_reaction_ids = set(self.draft_reaction.reactions)
        if not self.objectiveName in self.draft_reaction_ids:
            raise Exception("Objective not found in draftModel: {}".format(self.objectiveName))

        # Build a Reaction object for the exchange reactions;
        # if you have a defined medium, set the fixed_bounds argument accordingly
        self.exchange_reacs =   Reaction(model = self.path_to_exchanges, fixed_bounds = self.medium)
        self.db_reactions   =   Reaction(biochem_input = self.path_to_biochem, dbType=self.dbType)

        # Merge reactions from db with those of the draft model
        self.all_reactions           = Reaction(fixed_bounds = self.medium, dbType=self.dbType)
        self.all_reactions.reactions = self.all_reactions.add_dict(self.exchange_reacs.reactions, self.db_reactions.reactions)
        self.all_reactions.reactions = self.all_reactions.add_dict(self.draft_reaction.reactions, self.all_reactions.reactions)

        if self.medium_file and self.medium:
            raise Exception("provide medium or path, not both") #If you give both it would be confusing
        if self.medium_file:
            print("Loading medium from: {}".format(self.medium_file))
            self.medium = self.load_medium()
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
            self.NN = NN(path = self.trainedNNPath)
            self.predicted_reactions = self.NN.predict( self.draft_reaction_ids )
            for p_reaction in self.predicted_reactions:
                self.weights[p_reaction]  = np.round(1-self.predicted_reactions[p_reaction], 10)


        # Add reactions from all_reactions to candidate_reactions, with cost = default_cost.
        for reaction in self.all_reactions.reactions:
            if (reaction not in self.draft_reaction_ids) and (reaction not in self.weights):
                self.weights[reaction] = default_cost

        if self.black_list is not None:
            if self.draft_reaction_ids.intersection(self.black_list):
                raise Exception("Blacklist cannot contain reactions from the draft model")
            self.remove_candidates(self.black_list)

        if self.grey_list is not None:
            self.set_weights({k:self.punish_cost for k in self.grey_list})

        # Delete reaction from candidate_reactions if it is present in the starting model.
        rem_draf = []
        for reaction in self.weights:
            if reaction in self.draft_reaction_ids:
                rem_draf.append(reaction)
        for reaction in rem_draf:
            del self.weights[reaction]

        #TEMP FIX for candidates not in db:
        ip = []
        for i in self.weights:
            if i not in self.all_reactions.reactions.keys():
                ip.append(i)
                print("WARNING: reaction {} in candidates but not in database, removing from candidates".format(i))
        for i in ip:
            del self.weights[i]

        if gapfill:
            self.gapfill()

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
        [metabolite_id] : {[reaction_id] : stoichiometry}

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
            if self.dbType == 'ModelSEED':
                if '_c0' in metab:
                    cobra_metabs[metab].compartment = 'c0'

                elif '_e0' in metab:
                    cobra_metabs[metab].compartment = 'e0'
            #Note the BiGG database has many compartments, I dont know if this works very well
            else:
                if metab[-2] == '_':
                    cobra_metabs[metab].compartment = metab[-1]

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

        # Make cobra metabolites and reaction objects.
        cobra_metabs = self.make_cobra_metabolites(metab_dict)
        cobra_reactions = [self.make_cobra_reaction(reaction_dict, cobra_metabs, e) for e in reactions_in_model]

        # Make cobra model.
        cobra_model = cobra.Model('tempmodel')
        cobra_model.add_reactions(cobra_reactions)
        cobra_model.objective = objective_name
        return cobra_model

    def gapfill(self,
                result_selection = "min_reactions",
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

        all_reacs = deepcopy(self.all_reactions.reactions)
        all_reacs_obj = Reaction()
        all_reacs_obj.reactions = all_reacs.copy()
        cand_reacs = self.weights.copy()
        self.result_selection = result_selection

        # Split bidirectional reactions into a forward and reverse reaction.
        all_reactions_split = Reaction()

        all_reactions_split.reactions = all_reacs_obj.split_all_bidirectional_reactions(all_reacs_obj.reactions)

        # Add reverse reactions to draft_reaction_ids_split and candidate_reactions.
        draft_reaction_ids_split = set()

        for reaction in all_reactions_split.reactions:
            forward_version = reaction.replace('_rv', '')

            if forward_version in self.draft_reaction_ids:
                draft_reaction_ids_split.add(reaction)

            # If forward version of a reverse reaction is in candidate_reactions.
            else:
                if '_rv' in reaction:
                    # Give reverse reaction same cost as forward version.
                    cand_reacs[reaction] = cand_reacs[forward_version]


        # Run gapfilling algorithm
        split_gapfill_result = self.binarySearch(all_reactions_split,
                                                 draft_reaction_ids_split,
                                                 cand_reacs,
                                                 self.objectiveName,
                                                 )

        # If the database and media did not result in a functional model
        if split_gapfill_result is None:
            # raise ValueError
            print("Media is too restrictive. No growing model can be found :( \n\n\n")
            return None


        gapfill_result = set([r.replace('_rv', '') for r in split_gapfill_result])

        self.added_reactions = set(gapfill_result) #All reactions that are added to the model during gapfilling.

        gapfill_result.update(self.draft_reaction_ids) # ?? Add the reactions from the draft model to the gapfill result?

        # Create cobra model
        metab_dict      = all_reacs_obj.get_gurobi_metabolite_dict()
        cobra_model     = self.make_cobra_model(all_reacs,
                                                metab_dict,
                                                gapfill_result,
                                                self.objectiveName)
        self.objective_value = cobra_model.optimize().objective_value
        print ('Objective value is %f.' %self.objective_value)

        # Refine model based on the gapfill findings
        if self.medium is not None:
            self.gapfilledModel = build_model.refine_model(cobra_model,
                                                           self.draftModel,
                                                           unscalled = list(self.medium.keys()), dbType=self.dbType)
        else:
            self.gapfilledModel = build_model.refine_model(cobra_model, self.draftModel, dbType=self.dbType)

        self.gapfilledModel.id = "{}_gapfilled".format(self.draftModel.id)
        print("NN gapfilling added {} new reactions".format(len(self.added_reactions)))
        if self.grey_list is not None:
            added_grey = len(self.added_reactions.union(self.grey_list))
            print('{} of which were in the grey list'.format(added_grey))

        print("The NN gapfilled model, comes with {} reactions and {} metabolites".format(len(self.gapfilledModel.reactions), len(self.gapfilledModel.metabolites)))
        return self.gapfilledModel

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
        beta  = 2 * len(M)

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

            delta = int((alpha + beta) / 2.0)

            gu_model = self.build_gurobi_model(reaction_dict, metabolite_dict, B, N, M, delta)

            if np.round(gu_model.getVarByName(B).X,6) > 0:
                # print ('Flux through biomass reaction is {:.8f}'.format(gu_model.getVarByName(B).X))

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
                # print ('Flux through objective reaction is {:.8f}'.format(gu_model.getVarByName(B).X))
                #alpha = delta
                pass
            print('\n\n', 'condition is currently: ', abs(alpha - beta), '\n\n')

        # List that has minimum nr of reactions, sum of cost or sum of flux, dependend on output.
        minimum_set = self.get_minimum(R, R_flux, R_cost, criteria = self.result_selection)

        return  minimum_set

    def remove_candidates(self, list_to_remove):
        """
        Function to manually remove candidates from all_reactions and weights.
        Input:
            list of reactions to remove from all reactions
        Output:
            list of reactions with reactions removed
        """
        for i in list_to_remove:
            if i in self.all_reactions.reactions.keys():
                del self.all_reactions.reactions[i]
            if i in self.weights.keys():
                del self.weights[i]

    #Function to make it a bit easier to load a medium
    #consider removing the dependence on pandas
    def load_medium(self, e_pf='_e0'):
        df = read_csv(self.medium_file, sep='\t')       #load df
        df['exchanges'] = 'EX_'+df['id']+e_pf           #create exchange_ids
        df2 = df.set_index('exchanges')                 #dictionary is easier
        med = {}
        for ex in df2.index:
            med[ex] = {'lower_bound':df2['minflux'][ex], 'upper_bound':df2['maxflux'][ex], 'metabolites':{ex[3:]:-1.0}}
        return med

    def set_weights(self, custom_weights):
        """
        Function to manually set weights.
        Input:
            Dictionary of reactions to weights
        Output:
            None
        """
        for k in custom_weights:
            if k in self.all_reactions.reactions:
                self.weights[k] = custom_weights[k]
            else:
                print("{} not found in db, ignored".format(k))

    #This can be used to reset the weights to NN predictions or full default
    def reset_weights(self, no_NN=False):
        self.weights={}
        if not no_NN:
            self.predicted_reactions = self.NN.predict( self.draft_reaction_ids )
            for p_reaction in self.predicted_reactions:
                self.weights[p_reaction]  = np.round(1-self.predicted_reactions[p_reaction], 10)
        for reaction in self.all_reactions.reactions:
            if (reaction not in self.draft_reaction_ids) and (reaction not in self.weights):
                self.weights[reaction] = self.default_cost
