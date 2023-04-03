# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 23:16:07 2021

@author: u0139894
"""

import cobra
import numpy as np
import os
import ast


class Reaction:
    
    def __init__(self, model_folder=None, model_list=None, model=None, biochem_input=None, fixed_bounds=None):
        '''
        General class to handle reaction sets from metabolic models.
        
        The class can handle the ModelSEED biochemistry databases, an SBML model, or a list of SBML models.

        Parameters
        ----------

        model_folder : str, optional
            Path to the folder with sbml format models from where to retrieve reactions.
            By default all files in the folder are considered models.
            If only a a number files are desired, provide a model list.
        
        model_list : list, optional
            List of models in sbml format from where to retrieve reactions.
            If no model_folder is provided, assumes that each element in the list
            contains a full path, otherwise it assumes that the elements are correct
            file names in the 'model_folder'.
            
        model : str, optional
            Path to sbml model to retrieve reactions.
            
        biochem_input : str, optional
            Path to the biochemistry reactions.tsv file 
            from the modelSeed db. I other db is used, this function needs to be 
            adjusted as it references specific indices in the reactions table.
        

        Returns
        -------
        None.

        '''
        
        self.fixed_bounds = fixed_bounds
        
        self.model_folder = model_folder
        
        self.model_list = model_list
        
        self.model = self.__get_model(model)
        
        self.biochem_input = self.__get_reactions_from_biochem_input(biochem_input)
        
        self.reactions = self.__get_reactions()


    def __get_model(self, model):
        
        '''
        load a cobra model obj
        
        parameters
        ---------
        model: str
        path to model
        
        returns
        -------
        model_obj
        '''
        
        if model is None:
            return None
        
        else:
            self.model=cobra.io.read_sbml_model(model)
            return self.model
    
    def __get_reactions_from_model(self, model):
        '''
        retrieve reactions from a cobra model obj.
        
        reactions[reaction_id]={upper_bound:, lower_bound:, metabolites:{metabolite:stoichiometry}}
        '''
        reactions={}
        for reaction in model.reactions:
            
            reactions[reaction.id] = {'lower_bound': (reaction.lower_bound and abs(reaction.lower_bound)/reaction.lower_bound or 0), 'upper_bound': (reaction.upper_bound and abs(reaction.upper_bound)/reaction.upper_bound or 0)}
            mets = reaction.metabolites
            reactions[reaction.id]['metabolites']={i.id:mets[i] for i in mets}
                
        return reactions
    
    def __get_reactions_from_biochem_input(self, biochem_input):
        '''
        reactions[reaction_id]={upper_bound:, lower_bound:, metabolites:{metabolite:stoichiometry}}
        '''
        if biochem_input is None:
            return None
        
        react_d = self.read_biochemistry_table(biochem_input)
        reactions = {}
        for reaction in react_d:
            
            #I am assuming bacterial models (c0 and e0 compartments) and modelSEED identifiers
            #as 'rxn00001'
            direction = react_d[reaction][1]
            
            
            reaction_id = reaction[0:8] + "_c0" 
            
            metabolites = react_d[reaction][0].split(';')
            
            
            #boundaries set to 1 and -1 make gapfilling cleaner,
           
            #the external media should also be in this range.
            
            if direction == '=':
                reactions[reaction_id] = {'lower_bound':-1.0, 'upper_bound':1.0}
                
            else:
                reactions[reaction_id] = {'lower_bound':0.0, 'upper_bound':1.0}
                
            
            mets={}
            
            for i in metabolites:
                
                    metabolite = i.split(':')
                    
                    if len(metabolite) > 1: #Check if metabolite is specified.
                        stoc = metabolite[0]
                        cpd = metabolite[1]
                        loc = metabolite[2]
                        if int(loc) == 0:
                            name = cpd + '_c0'
                            
                        elif int(loc) == 1:
                            name = cpd + '_e0'
                        
                        else: #if the loc of the metabolite is not specified, define it as cellular
                            name = cpd + '_c0'
                            
                        if direction == '<':
                            mets[name] = -float(stoc) #Flip the sign of the stoichiometry of reversed reactions.
                                
                        else:
                            mets[name] = float(stoc)
            reactions[reaction_id]['metabolites'] = {i:mets[i] for i in mets}
                    
               
        return reactions
    
    def add_dict(self, dict_x, dict_y):
        '''
        utility function to add the items of dict_y to dict_x. If the values of duplicated itemts diverge,
        the original dict_x is retained.
        '''
        _d = dict_y.copy()
        _d.update(dict_x)
        
        if self.fixed_bounds is not None: 
            for reaction in self.fixed_bounds:

                if reaction in _d: 
                    _d[reaction]['lower_bound'] = self.fixed_bounds[reaction]['lower_bound']
                    _d[reaction]['upper_bound'] = self.fixed_bounds[reaction]['upper_bound']
        
        return _d
    
    def __get_reactions(self):
        '''
        make the reaction dict based on the inputs to the class
        
        reactions[reaction_id]={upper_bound:, lower_bound:, metabolites:{metabolite:stoichiometry}}
        '''
        reaction_dict={}
        if self.model is not None:
            model_reactions = self.__get_reactions_from_model(self.model)
            reaction_dict = self.add_dict(reaction_dict, model_reactions)
            
        if self.model_folder is not None:
            path = self.model_folder
            files_in_folder = os.listdir(self.model_folder)
        else:
            path = None
        
        models = []
        
        if self.model_list is not None:
            
            if path is not None:
                models = [os.path.join(path, i) for i in self.model_list]
            else:
                models = self.model_list[:]
        
        else:
            if path is not None:
                models = [os.path.join(path, i) for i in os.listdir(path)]
        
        for i in models:
                mod=cobra.io.read_sbml_model(self.model_folder+i)
                mod_reac=self.__get_reactions_from_model(mod)
                reaction_dict= self.add_dict(reaction_dict, mod_reac)
        
            
        if self.biochem_input is not None:
            reaction_dict = self.add_dict(reaction_dict, self.biochem_input)
            
        
        if self.fixed_bounds is not None: 

            for reaction in self.fixed_bounds:

                if reaction in reaction_dict: 
                    reaction_dict[reaction]['lower_bound'] = self.fixed_bounds[reaction]['lower_bound']
                    reaction_dict[reaction]['upper_bound'] = self.fixed_bounds[reaction]['upper_bound']
            
        self.reactions = reaction_dict
        return self.reactions
    
    def get_gurobi_reaction_dict(self, list_of_reactions=None):
        '''
        Convert the reactions dict to a input for a gurobi model, for the reactions in a list_of_reactions, or all reactions in self.reactions.
        
        Reaction_dict = { reaction_id : (lower_bound, upper_bound) }
        '''
        if self.reactions is None:
            return None
        
        if list_of_reactions is None:
            return {i:(self.reactions[i]['lower_bound'],self.reactions[i]['upper_bound'] ) for i in self.reactions}
            
        else:
            return {i:(self.reactions[i]['lower_bound'],self.reactions[i]['upper_bound'] ) for i in list_of_reactions}
    
    def get_gurobi_metabolite_dict(self, list_of_reactions=None):
        '''
        Get a metabolite dict for use as input for a gurobi model, 
        for reactions in list_of_reactions, or all reactions in self.reactions.
        
        Metab_dict = { metabolite_id : { reaction_id : stoichiometry } }
        '''
        if self.reactions is None:
            return None
        
        metab_dict = {}
        if list_of_reactions is None:
            for reaction in self.reactions:
                for met in self.reactions[reaction]['metabolites']:
                    if met not in metab_dict:
                        metab_dict[met] = {reaction : self.reactions[reaction]['metabolites'][met]}
                        
                    else:
                        metab_dict[met][reaction] = self.reactions[reaction]['metabolites'][met]
        
        else:
            for reaction in list_of_reactions:
                for met in self.reactions[reaction]['metabolites']:
                    if met not in metab_dict:
                        metab_dict[met] = {reaction : self.reactions[reaction]['metabolites'][met]}
                        
                    else:
                        metab_dict[met][reaction] = self.reactions[reaction]['metabolites'][met]
        
        return metab_dict
    
    def split_bidirectional_reaction(self, reaction):
        
        if float(reaction['lower_bound']) == -1 and float(reaction['upper_bound']) == 1: #Check if bidirectional.
            reaction_f = {k:v for k,v in reaction.items()} #Deep copy of reaction.
            reaction_f['lower_bound'] = 0.0 #Make reaction unidirectional.
            reaction_r = {k:v for k,v in reaction_f.items()} #Deep copy of forward reaction.
            reaction_r['metabolites'] = {i:-reaction['metabolites'][i] for i in reaction['metabolites']}
            
            return reaction_f, reaction_r
        
        else:
            return reaction, False
        
    def split_all_bidirectional_reactions(self, reaction_dictionary):
        '''
        Split all bidirectional reactions in a dictionary (of similar structure to self.reactions), or all reactions in self.reactions.
        Splitted reactions become two unidirectional reactions, reaction_id and reaction_id_r.
        When using a dictionary, this does not change the input dictionary, only the output dictionary has splitted reactions.
        When no dictionary is provided to the function, the 
        reactions of the self.reactions dict are changed.
        '''
        
        
        _d = {}
        
        for reaction in reaction_dictionary:
            r1, r2 = self.split_bidirectional_reaction(reaction_dictionary[reaction])
            _d[reaction] = r1
            if r2:
                r2_name = reaction + '_r'
                _d[r2_name] = r2
                    
        return _d
        
    def read_biochemistry_table(self, filePath):
        """
        Reads a ModelSEED/PATRIC reaction database table, using the id, stoichiometry, direction and equation columns.
        Returns a dictionary where id maps to a list of the other columns.
        """
        biochem = {}
        
        with open(filePath, 'r') as f:
            f.readline() #Skip the header
            for line in f:
                splitted_line = line.strip().split("\t")
                biochem[splitted_line[0]] = [splitted_line[4], splitted_line[9], splitted_line[6]]
            
        return biochem
