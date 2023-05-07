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
from tensorflow.python.keras.saving import hdf5_format
import h5py

class NN:
    def __init__(self,modeltype=None, path=None, custom=None):
        '''
        Wrapper for the Sequential keras model
        used to make predictions for models

        PARAMETERS
        -------
        modeltype: string, optional
            options are currently (ModelSEED or BiGG)
        path: TYPE, optional
            path to NN (optional will default based on modeltype)
        custom: tuple, list, optional
            temporary, I use this so NN_trainer can give a NN class as output
            might not be neccesary

        '''
        if custom:
            self.network = custom[0]
            self.rxn_keys = custom[1]
            self.modeltype = custom[2]
        else:
            if path is None:
                #this will only work if ran from scripts, should be a package way right?
                cwd = Path.cwd()
                sys.path.append(cwd)
                self.def_path = os.path.join(cwd.parent, 'files', 'NN')
                if modeltype:
                    self.modeltype = modeltype
                    if self.modeltype == 'ModelSEED':
                        path = os.path.join(self.def_path, 'NN_MS.h5')
                    elif self.modeltype == 'BiGG':
                        path = os.path.join(self.def_path, 'NN_BG.h5')
                else:
                    raise Exception('No path or recognised modeltype provided')
                print('Loading default {} NN'.format(self.modeltype))
            else:
                print('Loading network at user provided path')
            self.__get_network(path)



    #Function that loads the Neural network; path is path to .h5 file
    def __get_network(self, path):
        """
        Get network, modeltype and rxn_keys at path
        """
        with h5py.File(path, mode='r') as f:
            self.modeltype = f.attrs['modeltype']
            self.rxn_keys = f['rxn_keys'].asstr()[:]
        self.network = tf.keras.models.load_model(path, custom_objects={"custom_loss": 'binary_crossentropy'})

#
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
        if isinstance(input, cobra.core.model.Model):
            input = self.__convert_reaction_list(input.reactions.list_attr('id'))
        else:
            # check if DataFrame, reindex based on self.rxn_keys
            if (isinstance(input, pd.DataFrame)):
                input.reindex(self.rxn_keys)
                df_columns = input.columns
                #Transpose because rows need to be different models for the network
                input = input.T
            else:
                #check if dictionary, get list of reactions and convert
                if isinstance(input, dict):
                    input = self.__convert_reaction_list([i for i in input if input[i]==1])
                else:
                    #check if list or set, convert
                    if isinstance(input, (list, set)):
                        print('Converting to binary array:')
                        try:
                            input = self.__convert_reaction_list(input)
                        except:
                            raise Exception("Conversion failed")
                    #finally if input is already an array of 0s and 1s, Transpose to be the same as DataFrame with columns=models and rows=reactions
                    elif np.isin(input, [0,1]).all():
                        input = np.asarray(input.T)
                    else:
                        raise Exception("input type")

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
            raise Exception("data has wrong shape: ", input.shape, 'instead of ', self.network.input_shape)
        # if DataFrame is given as input, output DataFrame
        if isinstance(input, pd.DataFrame):
            prediction = pd.DataFrame(index=self.rxn_keys, columns=df_columns, data=prediction.T)
        return prediction

    #function that generates a binary input based on a list of reaction ids
    def __convert_reaction_list(self, reaction_set):
        """
        function that can be used to convert a set of reactions to a binary list with the order self.rxn_keys
        Input = the set of reactions that needs to be converted
        """
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
        print("#reactions not found in keys: ", len(set(reaction_set)) - sum(b_input), '/', len(reaction_set))
        if(sum(b_input)==0):
            raise Exception("No reactions found")
        return np.array(b_input)
