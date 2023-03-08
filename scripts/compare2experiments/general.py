# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:09:43 2023

@author: danie
"""

import os
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')
import seaborn as sns

import plotly.io as pio
pio.renderers.default='browser'

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))


from mainClasses import *
from parseTable import *
from updateParameters import *
from readModelDB import *
from loadParameters import *


###setup####

def simulateExperiment(species, experimentLabel, params, dbPath, measuredStates):
    
    ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv')
    
    strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)
    
    
    initialStates = {i:get_initialState(i, strainSummaryFolder, experimentLabel) for i in measuredStates}

    conn = create_connection(dbPath)
    
    assignParameters2db(species, params, conn)     
        
    states = []
    
    db = get_database(dbPath)
    
    

    wc = createMetabolome(db, 'wc')
    
    
    for state in initialStates:
        if state!='live':
            wc.metD[state].update(initialStates[state])


    predictpH = getpH(wc.metabolites, ipH_path)
    
    pH =  predictpH(wc.get_concentration())


    wc_f = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    
    
    for state in initialStates:
        if state!='live':
            wc_f.metD[state].update(initialStates[state])
    
    

    wc_r = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    
    
    for state in initialStates:
        if state!='live':
            wc_r.metD[state].update(initialStates[state])
    



    species_f = Microbiome({species:createBacteria(db, species, 'wc')})
    if species == 'bh':
        species_f.subpopD['bh.expA'].count = 0
    elif species == 'ri':
        species_f.subpopD['ri.lag'].count = 0
    elif species == 'bt':
        species_f.subpopD['bt.lag'].count = 0
        

    species_r = Microbiome({species:createBacteria(db, species, 'wc')})
    if species == 'bh':
        species_r.subpopD['bh.expA'].count = initialStates['live']
    elif species == 'ri':
        species_r.subpopD['ri.lag'].count = initialStates['live']
    elif species == 'bt':
        species_r.subpopD['bt.lag'].count = initialStates['live']


    p1 = Pulse(wc_f, species_f, 0, 120, 10000, 0, 0, 0, 0)

    r_species = Reactor(species_r, wc_r, [p1], 60)

    r_species.simulate()
    
    return r_species

def makeExperimentPlot(species, 
             state, 
             stateType = 'metabolite',
             experiments = ['bhbtri'],
             lables = ['bh3'],
             colors = ['#ff0000'],
             simulObj = [None, None, None],
             alpha=1):
    
    fig, ax = plt.subplots()
    strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)
    
    stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
            
    for i,v in enumerate(experiments):
        stateDF = getDFdict(stateTable, state, True)[v]
        sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = colors[i], lw=.666, label=lables[i], alpha=alpha)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})
        
        if simulObj[i] is not None:
            
            if state == 'pH':
                ax.plot(simulObj[i].time_simul, simulObj[i].pH_simul, color=colors[i], label=lables[i] + ' simul', linestyle='--', lw=3)
            
            elif state == 'live':
                
                ax.plot(simulObj[i].time_simul, simulObj[i].cellActive_dyn[0], color=colors[i], label=lables[i] + ' simul', linestyle='--', lw=3)
            
            elif state == 'dead':
                
                ax.plot(simulObj[i].time_simul, simulObj[i].cellInactive_dyn[0], color=colors[i], label=lables[i] + ' simul', linestyle='--', lw=3)
                
            else:
                
                ax.plot(simulObj[i].time_simul, simulObj[i].met_simul[simulObj[i].metabolome.metabolites.index(state)], color=colors[i], label=lables[i] + ' simul', linestyle='--', lw=3)
                
           

    legend_properties = {'size':12}
    ax.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))

    if stateType =='metabolite':
        ax.set_ylabel(state + ' (mM)', fontsize=16)
        plt.title(state, fontsize=16)
    
    elif stateType =='cells':
        ax.set_ylabel('$10^5$ cells/uL', fontsize=16)
        plt.title(state + ' cells', fontsize=16)
    else:
        ax.set_ylabel('pH', fontsize=16)
        plt.title('pH', fontsize=16)
    ax.set_xlabel('Time (h)', fontsize=16)

    
    plt.tight_layout()

    return fig,ax



# species = 'bh'
# experiments = ['bhbt', 'bhri', 'bhbtri']
# labels = ['bh1', 'bh2', 'bh3']
# colors = ['#00ff26', '#003eff', '#ff0000']


# mod1 = simulateExperiment(species, experiments[0], 'bh_bhbt.txt', 'modelDB_bhbtri_bh1 - Copy.sqlite3', ['live', 'trehalose', 'glucose', 'pyruvate', 'lactate','acetate'])

# mod2 = simulateExperiment(species, experiments[1], 'bh_bhri.txt', 'modelDB_bhbtri_bh2 - Copy.sqlite3', ['live', 'trehalose', 'glucose', 'pyruvate', 'lactate','acetate'])

# mod3 = simulateExperiment(species, experiments[2], 'bh_bhbtri.txt', 'modelDB_bhbtri_bh3 - Copy.sqlite3', ['live', 'trehalose', 'glucose', 'pyruvate', 'lactate','acetate'])


# makeExperimentPlot(species, 'pH', 'pH', ['bhri'], labels, colors)
# makeExperimentPlot(species, 'live', 'cells', ['bhri'], labels, colors)
# makeExperimentPlot(species, 'trehalose', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
# makeExperimentPlot(species, 'glucose', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
# makeExperimentPlot(species, 'pyruvate', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
# makeExperimentPlot(species, 'lactate', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
# makeExperimentPlot(species, 'acetate', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
