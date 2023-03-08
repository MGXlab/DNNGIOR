# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:48:54 2023

@author: danie
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))

from general import *



species = 'bh'
experiments = ['bhbt', 'bhri', 'bhbtri']
labels = ['bh1', 'bh2', 'bh3']
colors = ['#00ff26', '#003eff', '#ff0000']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live',
          'dead',
          'pH',
          'trehalose',
          'pyruvate',
          'glucose',
          'acetate',
          'lactate']


stTypes = ['cells',
           'cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']



params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh0.tsv'))

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


r1 = simulateExperiment(species, 'bhbt', params, database, ['live',
                                                               'trehalose',
                                                               'pyruvate',
                                                               'glucose',
                                                               'lactate',
                                                               'acetate'])



params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh3.tsv'))

databaseName = 'modelDB_bhD.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')


database = os.path.join(databaseFolder, databaseName)


r2 = simulateExperiment(species, 'bhri', params, database, ['live',
                                                               'trehalose',
                                                               'pyruvate',
                                                               'glucose',
                                                               'lactate',
                                                               'acetate'])



r3 = simulateExperiment(species, 'bhbtri', params, database, ['live',
                                                               'trehalose',
                                                               'pyruvate',
                                                               'glucose',
                                                               'lactate',
                                                               'acetate'])



for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [r1, r2, r3], alpha=0.05)



