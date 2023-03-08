# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 17:32:56 2023

@author: danie
"""

from pathlib import Path
import os
import sys


sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))

from readModelDB import *
from mainClasses import *
from loadParameters import *

import plotly.io as pio
pio.renderers.default='browser'


ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv') 


databaseName = 'modelDB_bhA.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')


#update database params
conn = create_connection(os.path.join(databaseFolder, databaseName))

bh_params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh0.tsv'))

assignBhParams(bh_params, conn)





#Load database
db = get_database(os.path.join(databaseFolder, databaseName))


#getStarting pH
wc = createMetabolome(db, 'wc')
predictpH = getpH(wc.metabolites, ipH_path)
pH =  predictpH(wc.get_concentration())

#get the feed media and the reactor media
wc_f = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
wc_r = createMetabolome(db, 'wc', pH, pHFunc=predictpH)

#get the feed obj. Make it sterile
bhbtri_f = Microbiome({'bh':createBacteria(db, 'bh', 'wc'), 'bt':createBacteria(db, 'bt', 'wc'), 'ri':createBacteria(db, 'ri', 'wc')})
bhbtri_f.subpopD['bh.expA'].count = 0
bhbtri_f.subpopD['bt.lag'].count = 0
bhbtri_f.subpopD['ri.lag'].count = 0

#create the reactor obj, with starting populations
bhbtri_r = Microbiome({'bh':createBacteria(db, 'bh', 'wc'), 'bt':createBacteria(db, 'bt', 'wc'), 'ri':createBacteria(db, 'ri', 'wc')})
bhbtri_r.subpopD['bt.lag'].count = 0
bhbtri_r.subpopD['ri.lag'].count = 0

batchA = Pulse(wc_f, bhbtri_f, 0, 120, 10000, 0, 0, 0,0)

#simulate
r_bhbtri = Reactor(bhbtri_r, wc_r,[batchA], 60)
r_bhbtri.simulate()
r_bhbtri.makePlots()


