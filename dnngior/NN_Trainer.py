# -*- coding: utf-8 -*-
"""
Created on Mon 24 Jan 2022
@author: meine
"""
from tensorflow.python.keras.saving import hdf5_format
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Input
from tensorflow.keras.optimizers.legacy import Adam
from tensorflow.keras import backend as K
from tensorflow.keras.metrics import AUC, Precision, Recall

import h5py
import numpy as np
import pandas as pd
import os
import sys
from tensorflow import compat, config, dtypes

from dnngior.NN_Predictor import NN

# Tensorflow; please consider: https://www.tensorflow.org/api_docs/python/tf/compat/v1/disable_eager_execution
compat.v1.disable_eager_execution()

#FUNCTIONS
def noise_data(i, noise_0, noise_1, del_p, con_p):
    """
        Function to randomly change 0s and 1s.

        Parameters
        ----------
        i : numpy array, required
            an array of 0s and 1s which you want to noisify
        noise_0 : numpy array, requiredimport dnngior.NN_Predictor
            fraction of 0s to change to 1s
        noise_1 : numpy array, required
            fraction of 1s to change to 0s
        del_p : TYPE,
            deletion probability that can be used for weighted changes
        con_p : TYPE,
            contamination probability that can be used for weighted changes

        Returns
        -------
        numpy array
            an array with 0s and 1s changed

    """
    temp = i.copy()
    #get index of all positions with a 1 and 0
    a=np.arange(len(temp))[temp==0]
    b=np.arange(len(temp))[temp!=0]
    n0 = int(len(a)*noise_0)
    n1 = int(len(b)*noise_1)
    #can use weighted deletion or not if del_p == None
    if con_p is None:
        #shuffle and get the first n indeces
        np.random.shuffle(a)
        s0 = a[:n0]
    else:
        con_p_norm = con_p[a]/np.sum(con_p[a])
        s0 = np.random.choice(a,n0, replace=False,p=con_p_norm)
    if del_p is None:
        np.random.shuffle(b)
        s1 = b[:n1]
    else:
        del_p_norm = del_p[b]/np.sum(del_p[b])
        #select n indeces with biased choice
        s1 = np.random.choice(b,n1, replace=False,p=del_p_norm)
    #change selected 1s to 0s and vice versa
    temp[s0]= 1
    temp[s1]= 0
    o = temp
    return o

def generate_feature(data, nuplo, min_con, max_con, min_for, max_for, del_p, con_p):
    """
    Function to generate the dataset for training (feature).
        PARAMETERS:
        ----------
        data: array,
            This is the data based on which a training dataset (features and labels) will be created.
            Needs to be array of 0s and 1s corresponding to reaction sets of metabolic models
        nuplo: int
            create duplicates of input data
            default=30

        The omission and contamination rates will increase linearly from min to max,
        with stepsize determined by nuplo
        min_for, float
            minimum false omssion rate, default = 0.05
        max_for, float
            maximum false ommision rate, default = 0.55
        min_con, float
            minimum contanimation introduced, currently not in use, default = 0
        max_con, float
            maximum contamination introduced, currently not in use, default = 0
        del_p, list
            list of probabilities of deletion for reactions
        con_p, list
            list of probabilities of introduction for reactions

        Returns:
        -------------
        train_data: TYPE
            Data with reactions deleted or added for training
    """
    #copy training data nuplo times
    data_copy = np.repeat(np.copy(data), nuplo, axis=0)
    r, c = data_copy.shape
    train_data = np.zeros((r, c))
    con_train = min_con
    for_train = min_for
    #iterate over rows (species) and change presences of reactions
    for i in range(r):
        if not (i % data.shape[0]):
            for_train += (max_for-min_for)/nuplo
            con_train += (max_con-max_con)/nuplo
        train_data[i] = noise_data(data_copy[i],con_train,for_train, del_p, con_p)

    return train_data

def custom_weighted_loss(dI, bias, maskI):
    """
    Function to calculate weighted loss

    PARAMETERS:
    ----------
    dI : TYPE, required
        Input of the Neural Network
    bias: float, required
        Bias towards the 0 class
    maskI = boolean, required
        Determines wether the input positions are masked during loss calculation

    Returns:
    -------------
    custom_loss: TYPE
        A custom loss function for training
    """

    if not (K.int_shape(dI)[0] == 'None'):
        bf = K.ones(1, K.int_shape(dI)[0]) #array of ones with the size of dI
        if(maskI):
            w = bf - dI # 1-I
        else:
            w = bf
    def custom_loss(y_true, y_pred): #y_true is the label, y_pred is the prediction made by the network
        loss = w * K.binary_crossentropy(y_true, y_pred) # mask x loss
        return bias*(1-y_true)*loss+(1-bias)*y_true*loss # return the biased loss y_true are all cases where prediction shouold be 1, 1-y_true all cases where prediction should be one, can scale between these two classes
    return custom_loss

def train(data, modeltype,rxn_keys=None,labels = None,validation_split=0.0,nuplo=30, min_con=0, max_con=0, min_for=0.05, max_for=0.3, con_p=None, del_p = None, nlayers=1, nnodes=256,  nepochs=10, b_size=32, dropout=0.1, bias_0=0.3, maskI=True, save=True, output_path='dnngior_predictor.npz', return_history=False, return_full_network=False):
    """
        Most important function, creates actual NN, there are many optional parameters

        PARAMETERS:
        ----------
        data: DataFrame or array, required
            binary array of reactions presences, for DataFrame index is used as rxn_keys
            otherwise rxn_keys should be provided
        modeltype : string, (currently) required
            The modeltype of the training data,
        rxn_keys: list, optional
            Can be used if data is not a pandas dataframe but a numpy array. Default is None
        labels:
            User can specify labels, by default input data is used as labels

        TRAINING PARAMETERS:
        -------
        see generate_feature() ^

        NETWORK PARAMETERS
        -------------
        nlayers: int, optional
            number of hidden layers (layers that are not input or output)
            default=1
        nnodes: int, optional
            number of nodes per layer,
            default=256
        nepochs: int, optional
            how often the network needs to loop over all the data
            default=10
        b_size: int, optional
            batch_size (number of training examples that are simultaneously evaluated)
            default=32,
        dropout: float, optional
            parameter for training that can reduce overfitting
            default = 0.1,
        bias_0: float, optional
            default = 0.3,
        maskI: boolean, optional
            Determines wether the input positions are masked during loss calculation, default=True
            default=True
        validation_split: float, optional
            Splits the input data in training and validation
            default = 0 (no split)
            if used, the valdiation will not be used during training but instead to calculate scores after training

        SAVING PARAMETERS:

        save: boolean, optional
            Whether you want to save the network, default = True
        output_path: string,
            Where to save the network, file_extension that work are .h5 and .npz
            all other file_extensions defailt to npz (lite network)
            default='dnngior_predictor.npz'

        OPTIONAL RETURNS:

        return_history: boolean, optional
            If you want training history
            default = False
        return_full_network: boolean, optional
            if you want to return the lite_network or full tensorflow object
            default = False
       Returns:
        -------------
        trainedNN
            NN class containing network, rxn_keys and modeltype
        history: if history=True
            history of training
    """

    print("Num GPUs Available: ", len(config.list_physical_devices('GPU')))

    if os.path.exists(output_path):
        print("# WARNING: overwriting savefile")
    elif os.access(os.path.dirname(output_path), os.W_OK):
        print("Saving network at: {}".format(output_path))
    else:
        Exception("Can not save at: {}".format(output_path))

    if(isinstance(data, pd.DataFrame)):
        rxn_keys = data.index
        ndata = np.asarray(data, dtype=np.float32).T
    elif rxn_keys is None:
        raise(Exception('Provide DataFrame or rxn_keys'))

    #create feature from training data
    if(labels is None):
        feature = np.repeat(np.copy(ndata), nuplo, axis=0).astype(np.float32)
        print('using data as labels')
    else:
        if(isinstance(labels, pd.DataFrame)):
            rxn_keys = data.index
            nlabels = np.asarray(labels, dtype=np.float32).T
        else:
            nlabels = labels.astype(np.float32).T
        feature = np.repeat(np.copy(nlabels), nuplo, axis=0).astype(np.float32)
        print("using user provided labels")

    train_data = generate_feature(ndata, nuplo, min_con, max_con, min_for, max_for, del_p, con_p)

    print('dataset created')
    nmodels, nreactions = ndata.shape
    print('training on data with shape: {} with {} reactions'.format(train_data.shape, sum(train_data.flatten())))

    #Build sequential model, define architecture
    network = Sequential()
    network.add(Input((nreactions,)))
    for _ in range(nlayers):
        network.add(Dense(nnodes, activation='relu'))
        network.add(Dropout(dropout))
    network.add(Dense(nreactions, activation='sigmoid'))

    #compile the network, determine parameters and loss function
    network.compile(Adam(learning_rate=0.005, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.01, amsgrad=True),
                  loss=custom_weighted_loss(network.input, bias_0, maskI),
                  metrics=[AUC(),Precision(thresholds=0.5),Recall(thresholds=0.5)])
    #print summary of model
    network.summary()
    #train model, history can be used to observe training
    history = network.fit(train_data, feature, validation_split = validation_split, epochs = nepochs, shuffle=True, batch_size = b_size, verbose=1)
    pseudo_network = []
    for i in range(0, len(network.layers),2):
        pseudo_network.append(network.layers[i].get_weights())
    pseudo_network = np.asarray(pseudo_network, dtype=object)
    #save Network
    if(save):
        if(output_path.endswith('.h5')):
            with h5py.File(output_path, mode='w') as f:
                network.save(f)
                f.attrs['modeltype'] = modeltype
                f.create_dataset("rxn_keys", data =[n.encode("ascii", "ignore") for n in rxn_keys])
        else:
            if not output_path.endswith('.npz'):
                file_extension = output_path.split('.')[-1]
                print('{} not recognized, saving as .npz (lite) instead'.format(file_extension))
                output_path.replace(file_extension, '.npz')
            np.savez(output_path,network=pseudo_network, modeltype=modeltype,rxn_keys=rxn_keys)
    if return_full_network:
        trainedNN = NN(custom=[network,modeltype,rxn_keys])
    else:
        trainedNN = NN(custom=[pseudo_network,modeltype,rxn_keys])
    if return_history:
        return trainedNN, history
    else:
        return trainedNN
