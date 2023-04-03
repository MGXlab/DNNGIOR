# -*- coding: utf-8 -*-
"""
Created on Mon 24 Jan 2022
@author: meine
"""
from tensorflow.python.keras.saving import hdf5_format
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Input
from tensorflow.keras import optimizers
from tensorflow.keras import backend as K
from tensorflow.keras.metrics import AUC, Precision, Recall

import h5py
import numpy as np
import pandas as pd
import os
import sys
from tensorflow import compat, config, dtypes

import dnngior.NN_Predictor

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
        noise_0 : numpy array, required
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

def generate_training_set(data,nuplo, min_con, max_con, min_for, max_for, del_p, con_p):
    """
    Function to generate the dataset for training (feature).
        PARAMETERS:
        ----------
        data: array,
            DESCRIPTION
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
            minimum contanimation introduced, currently not used by me, default = 0
        max_con, float
            maximum contamination introduced, currently not used by me, default = 0
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

def train(data, modeltype,rxn_keys=None,labels = None,validation_split=0.0,nuplo=30, min_con=0, max_con=0, min_for=0.05, max_for=0.3, con_p=None, del_p = None, nlayers=1, nnodes=256,  nepochs=10, b_size=32, dropout=0.1, bias_0=0.3, maskI=True, save=False, name='noname', output_path='', return_history=False):
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
        see generate_training_set() ^

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

        SAVING PARAMETERS:

        save: boolean, optional
            Whether you want to save the network, default = False
        name: string, optional
            name of your network, default='noname'
        output_path: string,
            where output, default=''
        return_history: boolean, optional
            If you want training history

       Returns:
        -------------
        trainedNN
            NN class containing network, rxn_keys and modeltype
        history: type, if history=True
            history of training
    """

    print("Num GPUs Available: ", len(config.list_physical_devices('GPU')))

    if(isinstance(data, pd.DataFrame)):
        rxn_keys = data.index
        ndata = np.asarray(data, dtype=np.float32).T
    elif rxn_keys is None:
        raise(Exception('Provide DataFrame or rxn_keys'))

    #create feature and labels from training data
    if(labels is None):
        labels = np.repeat(np.copy(ndata), nuplo, axis=0).astype(np.float32)
        print('using data as labels')
    else:
        labels = np.repeat(np.copy(labels), nuplo, axis=0).astype(np.float32)
        print("using user provided labels")

    train_data = generate_training_set(ndata, nuplo, min_con, max_con, min_for, max_for, del_p, con_p)

    print('dataset created')
    nmodels, nreactions = ndata.shape
    print('training on data with shape: {} with {} reactions'.format(train_data.shape, sum(train_data.flatten())))

    #Build sequential model, define architecture
    model = Sequential()
    model.add(Input((nreactions,)))
    for _ in range(nlayers):
        model.add(Dense(nnodes, activation='relu'))
        model.add(Dropout(dropout))
    model.add(Dense(nreactions, activation='sigmoid'))

    #compile the network, determine parameters and loss function
    model.compile(optimizers.Adam(learning_rate=0.005, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.01, amsgrad=True),
                  loss=custom_weighted_loss(model.input, bias_0, maskI),
                  metrics=[AUC(),Precision(thresholds=0.5),Recall(thresholds=0.5)])
    #print summary of model
    model.summary()

    #train model, history can be used to observe training
    history = model.fit(train_data, labels, validation_split = validation_split, epochs = nepochs, shuffle=True, batch_size = b_size, verbose=1)
    #save Network
    if(save):
        model_path = os.path.join(output_path, "{}.h5".format(name))
        with h5py.File(model_path, mode='w') as f:
            model.save(f)
            f.attrs['modeltype'] = modeltype
            f.create_dataset("rxn_keys", data =[n.encode("ascii", "ignore") for n in rxn_keys])

    trainedNN = NN_Predictor.NN(custom=[model,rxn_keys,modeltype])
    if return_history:
        return trainedNN, history
    else:
        return trainedNN
