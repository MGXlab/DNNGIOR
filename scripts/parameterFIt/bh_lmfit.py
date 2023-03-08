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

#####################################################################



def getDistance(x,y, err=1, delta=2.0):
    size = int(len(x)*.5)
    
    k = np.random.choice(np.arange(len(x)), size=size)
    distance = (delta**2) * (((1+((y[k]-x[k])/delta)**2)**0.5)-1)
    

    return np.average(distance) #min(1.0, distance)




def pseudoHuberLoss(y_true, y_pred, delta = 1.0):
    """
    Compute the Pseudo-Huber loss between y_true and y_pred with a given delta value.
    """
    # Compute the squared error
    error = y_true - y_pred
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


####get initial values & splines
iValuesSplines = {i:(get_initialState(i, strainSummaryFolder, experimentLabel), get_spline(i, strainSummaryFolder, experimentLabel), get_StateStd(i, strainSummaryFolder, experimentLabel)) for i in measuredStates}


lmfit_params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh0.tsv'))


evals = []


def distance(lmfit_params, database = database):
    
    conn = create_connection(database)
    
    assignBhParams(lmfit_params, conn)
    
    db = get_database(database)
    
    wc = createMetabolome(db, 'wc')
    
    wc.metD['trehalose'].update(iValuesSplines['trehalose'][0])
    wc.metD['pyruvate'].update(iValuesSplines['pyruvate'][0])
    wc.metD['glucose'].update(iValuesSplines['glucose'][0])
    wc.metD['lactate'].update(iValuesSplines['lactate'][0])
    wc.metD['acetate'].update(iValuesSplines['acetate'][0])
    
    predictpH = getpH(wc.metabolites, ipH_path)
    pH =  predictpH(wc.get_concentration())
    
    wc_f = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    wc_f.metD['trehalose'].update(iValuesSplines['trehalose'][0])
    wc_f.metD['pyruvate'].update(iValuesSplines['pyruvate'][0])
    wc_f.metD['glucose'].update(iValuesSplines['glucose'][0])
    wc_f.metD['lactate'].update(iValuesSplines['lactate'][0])
    wc_f.metD['acetate'].update(iValuesSplines['acetate'][0])
    
    wc_r = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    wc_r.metD['trehalose'].update(iValuesSplines['trehalose'][0])
    wc_r.metD['pyruvate'].update(iValuesSplines['pyruvate'][0])
    wc_r.metD['glucose'].update(iValuesSplines['glucose'][0])
    wc_r.metD['lactate'].update(iValuesSplines['lactate'][0])
    wc_r.metD['acetate'].update(iValuesSplines['acetate'][0])
    
    bh_f = Microbiome({'bh':createBacteria(db, 'bh', 'wc')})
    bh_f.subpopD['bh.expA'].count = 0
    
    bh_r = Microbiome({'bh':createBacteria(db, 'bh', 'wc')})
    bh_r.subpopD['bh.expA'].count = iValuesSplines['live'][0]
    
    batchA = Pulse(wc_f, bh_f, 0, 120, 1000, 0, 0, 0, 0)
    
    r_bh = Reactor(bh_r, wc_r, [batchA], 60)
    
    r_bh.simulate()
    
    live_distance = pseudoHuberLoss(iValuesSplines['live'][1](r_bh.time_simul), r_bh.cellActive_dyn[0])
    
    dead_distance = pseudoHuberLoss(iValuesSplines['dead'][1](r_bh.time_simul), r_bh.cellInactive_dyn[0])
    
    pH_distance = pseudoHuberLoss(iValuesSplines['pH'][1](r_bh.time_simul), r_bh.pH_simul)
    
    trehalose_distance = pseudoHuberLoss(iValuesSplines['trehalose'][1](r_bh.time_simul), r_bh.met_simul[r_bh.metabolome.metabolites.index('trehalose')])
    
    pyruvate_distance = pseudoHuberLoss(iValuesSplines['pyruvate'][1](r_bh.time_simul), r_bh.met_simul[r_bh.metabolome.metabolites.index('pyruvate')])
    
    glucose_distance = pseudoHuberLoss(iValuesSplines['glucose'][1](r_bh.time_simul), r_bh.met_simul[r_bh.metabolome.metabolites.index('glucose')])
    
    lactate_distance = pseudoHuberLoss(iValuesSplines['lactate'][1](r_bh.time_simul), r_bh.met_simul[r_bh.metabolome.metabolites.index('lactate')])
    
    acetate_distance = pseudoHuberLoss(iValuesSplines['acetate'][1](r_bh.time_simul), r_bh.met_simul[r_bh.metabolome.metabolites.index('acetate')])
    
    distance = np.array([5*live_distance,
                         dead_distance,
                         pH_distance,
                         trehalose_distance,
                         pyruvate_distance,
                         5*glucose_distance,
                         lactate_distance,
                         acetate_distance])
    
    print(sum(distance))
    print(np.max(distance), '\t', np.mean(distance), live_distance)
  
    ssrMean = np.round(sum(distance),3)
        
    if len(evals)>0:    
        if ssrMean<min(evals):
            writeOutput(lmfit_params, outputFile)
            plt.plot(evals)
            plt.title(str(ssrMean))
            plt.show()
    evals.append(ssrMean)
    return sum(distance)
    
    
d = distance(lmfit_params, database)

out = minimize(distance, params=lmfit_params, method='bfgs')

