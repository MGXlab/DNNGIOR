"""
Created on Thu March 9 2023

@author: meine
"""
import numpy as np
import pandas as pd
import cobra.core.model as cobra_model
from dnngior.reaction_class import Reaction
from dnngior.variables import *
import os
import sys
from math import exp
from pathlib import Path

class NN:
    def __init__(self, path=None, modeltype=None, custom=None):
        '''
        Light version of the model, saves space, uses only numpy and cobra and no tensorflow
        '''

        if custom:
            self.network   = custom[0]
            self.modeltype = custom[1]
            self.rxn_keys  = custom[2]
        else:
            if path:
                self.path=path

            elif modeltype:
                if modeltype == 'ModelSEED':
                    self.path = TRAINED_NN_MSEED
                elif modeltype == 'BiGG':
                    self.path = TRAINED_NN_BIGG
                else:
                    print("Modeltype: {} not recognized, defaulting to ModelSEED".format(modeltype))
                    self.path = TRAINED_NN_MSEED
            else:
                print("No path or modeltype provided, defaulting to ModelSEED")
                self.path = TRAINED_NN_MSEED
            self.__get_pseudo_network()




            #Function that loads the Neural network; path is path to .h5 file
    def __get_pseudo_network(self):

        f = np.load(self.path, allow_pickle=True)

        self.network = f['network']
        self.modeltype = str(f['modeltype'])
        self.rxn_keys = f['rxn_keys']

    def predict(self, input):
        """
        Function that makes a prediction based on input_data
        Input can be several things:
        multiple reaction sets:
        DataFrame -> DataFrame with predictions with index = self.rxn_keys
        array -> array with same positions
        single reaction set:
        Cobra model
        Dictionary of reactions mapping to 0 or 1
        List or set of reactions
        all -> dictionary {reaction: prediction}
        Exception:
        if input shape does not match with the input of the network
        """
    #Check if model -> get list of ids and convert
        if isinstance(input, cobra_model.Model):
            input2 = self.__convert_reaction_list(input.reactions.list_attr('id'))
        elif isinstance(input, Reaction):
            #check if reaction class
            input2 = self.__convert_reaction_list(set(input.reactions))
        elif isinstance(input, pd.DataFrame):
            input1b = input.reindex(self.rxn_keys).fillna(0.0)
            df_columns = input.columns
            #Transpose because rows need to be different models for the network
            input2 = np.asarray(input1b.T)
        elif isinstance(input, dict):
            #check if dictionary, get list of reactions and convert
            input2 = self.__convert_reaction_list([i for i in input if input[i]==1])
        elif isinstance(input, (list, set)):
            #check if list or set, convert
            input2 = self.__convert_reaction_list(input)
        elif np.isin(input, [0,1]).all():
            input2 = np.asarray(input.T)
        else:
            raise Exception("input type")

        #test for single input (trips up NN)
        if np.ndim(input2) == 1:
            single_input=True
            input2 = np.expand_dims(input2,axis=0)
        else:
            single_input=False

        if not input2.shape[1] == self.network[0][0].shape[0]:
            raise Exception("Input size ({}) does not match network ({})".format(input2.shape[1], len(self.rxn_keys)))

        a = input2
        for layer in self.network:
            a = a.clip(0)
            a = ((a @ layer[0]) + layer[1])
        prediction =  1 / (1 + np.exp(-a))#sigmoid(a)

        if single_input:
            prediction = dict(zip(self.rxn_keys, np.squeeze(prediction)))
        if isinstance(input, pd.DataFrame):
            prediction = pd.DataFrame(index=self.rxn_keys, columns=df_columns, data=prediction.T)
            if len(prediction.index) != len(input.index):
                print('Warning mismatch input vs prediction ({})'.format(len(prediction.index) - len(input.index)))
        return prediction

    #function that generates a binary input based on a list of reaction ids
    def __convert_reaction_list(self, reaction_set):
        """
        function that can be used to convert a set of reactions to a binary list with the order self.rxn_keys
        Input = the set of reactions that needs to be converted
        """
        try:
            b_input = []
            #I think that at this point this might be the only reason I would need modeltype, there are few things I am as annoyed with as the _c0
            if(self.modeltype=='ModelSEED'):
                reaction_list = [reaction[:8]+'_c0' if 'rxn' in reaction else reaction for reaction in reaction_set]
                self.rxn_keys  = [key[:8]+'_c0' if 'rxn' in key else key for key in self.rxn_keys ]
            else:
                reaction_list = list(reaction_set)
            for i in self.rxn_keys:
                if i in reaction_list:
                    b_input.append(1)
                else:
                    b_input.append(0)
            print("#reactions not found in NN-keys: ", len(set(reaction_set)) - sum(b_input), '/', len(reaction_set))
        except:
            raise Exception("Conversion failed")

        if(sum(b_input)==0):
            raise Exception("No reactions found")
        return np.array(b_input)
