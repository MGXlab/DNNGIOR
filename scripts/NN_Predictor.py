#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 13:19:13 2019

@author: meine
"""
import numpy as np
import pandas as pd
import os
import sys
import tensorflow as tf
from tensorflow.keras.models import Sequential
from pathlib import Path
import cobra


class NN:
    def __init__(self,modeltype='ModelSEED', path=None, rxn_keys_path=None):
        '''
        Kind of a wrapper for the Sequential keras model
        used to make predictions for models

        modeltype = kind of model, options are currently (ModelSEED or CarveMe)
        path  = path to NN (optional will default based on modeltype)
        rxn_keys_path =  path to the keys (optional based on modeltype)

        '''
        
        self.modeltype = modeltype
        
        #If you dont give a path to your own keys or NN I define a default path
        if path is None or rxn_keys_path is None:
            cwd = Path.cwd()
            sys.path.append(cwd)
            self.def_path = os.path.join(cwd.parent, 'files', 'NN')

        if path is None:
            if self.modeltype == 'ModelSEED':
                path = os.path.join(self.def_path, 'NN_MS.h5')
            elif self.modeltype == 'CarveMe':
                path = os.path.join(self.def_path, 'NN_CM.h5')
            else:
                raise Exception('No path or recognised modeltype provided')
            print('Loading default {} NN'.format(self.modeltype))
        else:
            print('Loading network at user provided path')

        if rxn_keys_path is None:
            if(self.modeltype == 'ModelSEED'):
                rxn_keys_path = os.path.join(self.def_path,'rxn_ids_ModelSEED.npy')
            elif(self.modeltype == 'CarveMe'):
                rxn_keys_path = os.path.join(self.def_path,'/rxn_ids_BiGG.npy')
            print('Using {} ids at {}'.format(self.modeltype, rxn_keys_path))
        else:
            print('Using ids at provided path')

        self.network = self.__get_network(path)
        self.rxn_keys = self.__get_ids(rxn_keys_path)


    #Function that loads the Neural network; path is path to .h5 file
    def __get_network(self, path):
        network = tf.keras.models.load_model(path, custom_objects={"custom_loss": 'binary_crossentropy'})
        if not isinstance(network, Sequential):
            raise Exception('Type: {} not supported'.format(type(network)))
        return network


    def __get_ids(self, rxn_keys_path):
        ids = np.load(rxn_keys_path, allow_pickle=True).astype('str')
        return ids

#Function that makes a prediction based on input_data using Neural Network (NN)
    def predict(self, input):
        #check if input = model
        if isinstance(input, cobra.core.model.Model):
            input = self.__convert_reaction_list(input.reactions.list_attr('id'))
        else:
            #if input is a pandas dataframe, reindex based on keys
            if (isinstance(input, pd.DataFrame)):
                input.reindex(self.rxn_keys)
                df_columns = input.columns
                input = input.T
            else:
                #if input is a dictionary, take keys and convert that reaction list
                if isinstance(input, dict):
                    input = self.__convert_reaction_list([i for i in input if input[i]==1])
                else:
                    #if set or list create convert to binary
                    if isinstance(input, (list, set)):
                        print('Converting to binary array:')
                        try:
                            input = self.__convert_reaction_list(input)
                        except:
                            raise Exception("Conversion failed")
                    #finally if it is already a binary array, nothing needs to be done
                    elif np.isin(input, [0,1]).all():
                        input = np.asarray(input.T)
                    else:
                        raise Exception("input type: {}".format(type(input))

        single_input=False
        #test for single input (trips up NN)
        if np.ndim(input) == 1:
            single_input=True
            input = np.expand_dims(input,axis=0)
        if input.shape[-1] == self.network.input_shape[-1]:
            prediction = np.zeros(input.shape)
            #   load network and make prediction
            t_result = self.network.predict(input)
            prediction = np.asarray(t_result)

            if single_input:
                prediction = dict(zip(self.rxn_keys, np.squeeze(prediction)))
        else:
            raise Exception("data has wrong shape: ", input.shape, 'instead of ', NN.input_shape)
        if isinstance(input, pd.DataFrame):
            prediction = pd.DataFrame(index=self.rxn_keys, columns=df_columns, data=prediction.T)
        return prediction

    #function that generates a binary input based on a list of reaction ids
    def __convert_reaction_list(self, reaction_set):
        b_input = []
        #need to take into account the annoying _c0 stuff
        if(self.modeltype=='ModelSEED'):
            reaction_list = [reaction[:8] for reaction in reaction_set]
            self.rxn_keys  = [key[:8] for key in self.rxn_keys]
        else:
            reaction_list = list(reaction_set)
        for i in self.rxn_keys:
            if i in reaction_list:
                b_input.append(1)
            else:
                b_input.append(0)
        print("#reactions not found in keys: ", len(set(reaction_set)) - sum(b_input), '/', len(reaction_set))
        if(sum(b_input)==0):
            raise Exception("No reactions found")
        return np.array(b_input)
