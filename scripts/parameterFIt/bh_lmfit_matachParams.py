# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 08:54:49 2022

@author: u0139894
"""

from pathlib import Path
import os
import sys

from scipy.interpolate import PchipInterpolator as CubicSpline
from lmfit import minimize, Parameters, fit_report


sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))


from mainClasses import *
from parseTable import *
from updateParameters import *
from readModelDB import *
from loadParameters import *
from general import *


def pseudoHuberLoss(y_true, y_pred, delta = 0.50):
    """
    Compute the Pseudo-Huber loss between y_true and y_pred with a given delta value.
    """
    
    choice = np.random.choice(np.arange(len(y_true)), size = int(0.7*len(y_true)), replace=False)
    
    y_t = y_true[choice]
    y_p = y_pred[choice]
    
    # Compute the squared error
    error = y_t - y_p
    squared_error = np.square(error)

    # Compute the loss for small errors
    small_error = delta**2 * (np.sqrt(1 + squared_error / delta**2) - 1)

    # For values of error greater than delta, use the linear part of the Pseudo-Huber loss
    large_error = delta * (np.sqrt(squared_error) - 0.5 * delta)

    # Return the mean loss
    return np.mean(np.where(np.abs(error) <= delta, small_error, large_error))



####################################################################
def writeOutput(lmfit_params, outputFile):
    with open(outputFile, 'w') as f:
        f.write('species\tparameter\tvalue\tmin\tmax\n')
        
        for i in lmfit_params:
            f.write(species + '\t' + i + '\t' + str(lmfit_params[i].value) + '\t' + str(lmfit_params[i].min) + '\t' + str(lmfit_params[i].max) + '\n')




def distance(lmfit_params, database, initialStates, splines, experimentLabel, species, lbd=0.1):
    
    conn = create_connection(database)
    
    assignBhParams(lmfit_params, conn)
    
    db = get_database(database)
    
    r = simulateExperiment(species, experimentLabel, lmfit_params, database, initialStates)
    
    
    distances = []
    
    for i in initialStates:
        if i=='live':
            distances.append(pseudoHuberLoss(splines['live'](r.time_simul), r.cellActive_dyn[0]))
        
        elif i=='dead':
            
            distances.append(pseudoHuberLoss(splines['dead'](r.time_simul), r.cellInactive_dyn[0]))
        
        elif i=='pH':
            distances.append(pseudoHuberLoss(splines['pH'](r.time_simul), r.pH_simul))
        
        else:
            distances.append(pseudoHuberLoss(splines[i](r.time_simul), r.met_simul[r.metabolome.metabolites.index(i)]))
    
    
    pars = np.array([inputParams[i].value for i in inputParams])
    
    pars_std = (pars-np.mean(pars))/np.std(pars)
    
    objV = sum(distances) + sum(pars_std**2)*lbd
   
    ssrSum = np.round(objV,3)
    
    print(ssrSum)
    if len(evals)>0:    
        if ssrSum<min(evals):
            writeOutput(lmfit_params, outputFile)
            plt.plot(evals)
            plt.title(str(ssrSum))
            plt.show()
    evals.append(ssrSum)
    return objV
    


##############SETUP###########################################
ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv') 

species = 'bh'

experimentLabel = 'bhbt'

strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh')

inputParams = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh0.tsv'))





outputFile = os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh0.tsv')


databaseName = 'modelDB_bhA.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')


database = os.path.join(databaseFolder, databaseName)

measuredStates = ['live', 
                  'dead',
                  'pH',
                  
                  'trehalose',
                  'pyruvate',
                  'glucose',
                  'lactate',
                  'acetate']




initialStates = ['live',
                 'trehalose',
                 'pyruvate',
                 'glucose',
                 'lactate',
                 'acetate'
                 ]

#####################################################################

####get initial values & splines
splines = {}

splines['bhbt'] = {i:get_spline(i, strainSummaryFolder, 'bhbt') for i in measuredStates}
splines['bhri'] = {i:get_spline(i, strainSummaryFolder, 'bhri') for i in measuredStates}
splines['bhbtri'] = {i:get_spline(i, strainSummaryFolder, 'bhbtri') for i in measuredStates}







evals = []


def minF(inputParams, database, initialStates, species, lbd):
    
    l = ['bhbt']#, 'bhri', 'bhbtri']
    
    expL = np.random.choice(l)
    
    return distance(inputParams, database, initialStates, splines[expL], expL, species, lbd = lbd)


# for i in range(10):
    
#     print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
    
    

#     selection = np.random.choice(list(inputParams.keys()), size = len(measuredStates), replace=False)
    
    
#     for param in inputParams:
#         if param not in selection:
#             inputParams[param].vary = False
#         else:
#             inputParams[param].vary = True
    
    
# #d = 
    
out = minimize(minF, params=inputParams, method='powell', kws = {'database' : database,
                                                                 'initialStates' : initialStates,
                                                                 'species' : species,
                                                                 'lbd':1.0})
    
