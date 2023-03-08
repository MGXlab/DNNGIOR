# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 13:34:26 2023

@author: danie
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))
from general import *


###########BH###############


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

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + '.png'), dpi = 150)
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + '.png'), dpi = 50)

#####################################