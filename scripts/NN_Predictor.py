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
path = Path.cwd()
sys.path.append(path)
trainedNN = os.path.join(path.parent, 'files', 'NN')


#Function that loads the Neural network, maybe unneccesary; path is path to .h5 file
def load_NN(path=None):
    '''
    

    Parameters
    ----------
    path : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    if path is None:
        print('Loading Default NN (ModelSEED)')
        path = os.path.join(trainedNN, 'NN_MS.h5')
    else:
        print('Loading network at user provided path')
    return tf.keras.models.load_model(path, custom_objects={"custom_loss": 'binary_crossentropy'})

def load_ids(nnpath=trainedNN):
    '''
    

    Parameters
    ----------
    path : TYPE, optional
        DESCRIPTION. The default is trainedNN + '/rxn_ids_ModelSEED.npy'.

    Returns
    -------
    ids : TYPE
        DESCRIPTION.

    '''
    nfile = os.path.join(nnpath, 'rxn_ids_ModelSEED.npy') 
    ids = np.load(nfile, allow_pickle=True).astype('str')
    return ids
#Function that makes a prediction based on input_data using Neural Network (NN)
def predict(input, trainedNN=None, rxn_ids=None):
    '''
    

    Parameters
    ----------
    input : TYPE
        DESCRIPTION.
    NN : TYPE, optional
        DESCRIPTION. The default is None.
    rxn_ids : TYPE, optional
        DESCRIPTION. The default is None.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    prediction : TYPE
        DESCRIPTION.

    '''
    if rxn_ids is None:
        rxn_ids = load_ids(trainedNN)
    if isinstance(input, cobra.core.model.Model):
        input = convert_reaction_list(input.reactions.list_attr('id'), rxn_ids, trainedNN)
    else:
        if (isinstance(input, pd.DataFrame)):
            input.reindex(rxn_ids)
            df_columns = input.columns
            input = input.T
        else:
            if isinstance(input, dict):
                input = convert_reaction_list([i for i in input if input[i]==1], rxn_ids, trainedNN)
            else:
                if not np.isin(input, [0,1]).all():
                    print('Converting to binary array:')
                    try:
                        input = convert_reaction_list(input, trainedNN = trainedNN)
                    except:
                        raise Exception("Conversion failed")
                else:
                    input = np.asarray(input.T)

    if trainedNN is None:
        trainedNN = load_NN()
    if isinstance(trainedNN, str):
            print('Loading network')
            trainedNN = load_NN(os.path.join(trainedNN, 'NN_MS.h5'))
    if not isinstance(trainedNN, Sequential):
        raise Exception('Type: {} not supported'.format(type(trainedNN)))


    single_input=False
    #test for single input (trips up trainedNN)
    if np.ndim(input) == 1:
        single_input=True
        input = np.expand_dims(input,axis=0)
    if input.shape[-1] == trainedNN.input_shape[-1]:
        prediction = np.zeros(input.shape)
        #   load network and make prediction
        t_result = trainedNN.predict(input)
        prediction = np.asarray(t_result)

        if single_input:
            prediction = dict(zip(rxn_ids, np.squeeze(prediction)))
    else:
        raise Exception("data has wrong shape: ", input.shape, 'instead of ', trainedNN.input_shape)
    if isinstance(input, pd.DataFrame):
        prediction = pd.DataFrame(index=rxn_ids, columns=df_columns, data=prediction.T)
    return prediction

#function that generates a binary input based on a list of reaction ids
def convert_reaction_list(reaction_set, NN_reaction_ids = None, trainedNN = None):
    '''
    

    Parameters
    ----------
    reaction_set : TYPE
        DESCRIPTION.
    NN_reaction_ids : TYPE, optional
        DESCRIPTION. The default is None.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    if NN_reaction_ids is None:
        if(list(reaction_set)[0][:3] == 'rxn'):
            model_type = 'ModelSEED'
            NN_reaction_ids = np.load(trainedNN + '/rxn_ids_ModelSEED.npy', allow_pickle=True).astype('str')
        else:
            model_type = 'BiGG'
            NN_reaction_ids = np.load(trainedNN + '/rxn_ids_bigg.npy', allow_pickle=True).astype('str')
        print('Using {} ids'.format(model_type))
    else:
        print('Using user-provided ids')
    b_input = []
    if(list(reaction_set)[0][:3] == 'rxn'):
        reaction_list = [reaction[0:8] + "_c0" for reaction in reaction_set]
    else:
        reaction_list = list(reaction_set)
    for i in NN_reaction_ids:
        if i in reaction_list:
            b_input.append(1)
        else:
            b_input.append(0)
    cheat = set([i for i in set(reaction_set) if not 'EX' in i])
    print("#reactions not in NN_rxn: ", len(cheat.difference(NN_reaction_ids)))
    if(sum(b_input)==0):
        raise Exception("No reactions found")
    return np.array(b_input)
