# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 11:16:19 2022

@author: u0139894
"""

import os
from scipy.interpolate import PchipInterpolator as CubicSpline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')

def parseTable(filePath):
    d = {}
    
    with open(filePath) as f:
        f.readline()
        for line in f:
            a = line.strip().split('\t')
            t = float(a[0])
            v = float(a[3])
            experiment = a[1]
            replicate = a[2]
            
            if experiment not in d:
                d[experiment] = {}
            
            if replicate not in d[experiment]:
                d[experiment][replicate] = {'time':[], 'measure': []}
            
            d[experiment][replicate]['time'].append(t)
            d[experiment][replicate]['measure'].append(v)
        
        return d





def getDFdict(parsedTable, measureName='measure', plot=True):
    d = {}
    
    
    if plot:
        for experiment in parsedTable:
            r = {'time':[], measureName:[]}
            
            replicates = list(parsedTable[experiment].keys())
            
            for i,timeP in enumerate(parsedTable[experiment][replicates[0]]['time']):
                
                for replicate in parsedTable[experiment]:
                    
                    r['time'].append(timeP)
                    r[measureName].append(parsedTable[experiment][replicate]['measure'][i])
            
            df =  pd.DataFrame.from_dict(r)
            d[experiment] = df.copy()
    
    else:
        for experiment in parsedTable:
            r = {}
            for replicate in parsedTable[experiment]:
                t = np.array(parsedTable[experiment][replicate]['time'])
                v = np.array(parsedTable[experiment][replicate]['measure'])
                
                r[replicate] = v
            
            r['time'] = t
            
            df = pd.DataFrame.from_dict(r)
            df = df.set_index('time')
            d[experiment] =df.copy()
            
    return d

# bt_od = parseTable('C:\\Users\\u0139894\\Documents\\GitHub\\SyntheticCommunity\\files\\strainSummaries\\bt\\od.txt')

# df_od = getDFdict(bt_od, 'OD', False)

def get_spline(stateName, strainSummaryFolder, experimentLabel):
    
    state = parseTable(os.path.join(strainSummaryFolder, stateName + '.tsv'))

    df_state = getDFdict(state, stateName, False)[experimentLabel]
    


    timeV = np.array(df_state.index)
    
    stateMean = np.array(df_state.mean(axis=1))
    

    return CubicSpline(timeV, stateMean, extrapolate=False)
    
def get_initialState(stateName, strainSummaryFolder, experimentLabel):
    
    state = parseTable(os.path.join(strainSummaryFolder, stateName + '.tsv'))

    df_state = getDFdict(state, stateName, False)[experimentLabel]
    
    
    stateMean = np.array(df_state.mean(axis=1))
    

    return stateMean[0]

def get_StateStd(stateName, strainSummaryFolder, experimentLabel):
    
    state = parseTable(os.path.join(strainSummaryFolder, stateName + '.tsv'))

    df_state = getDFdict(state, stateName, False)[experimentLabel]
    
    stateStd = np.array(df_state.var(axis=1))
    
    sum(stateStd**2)
    

    return (sum(stateStd**2))**0.5