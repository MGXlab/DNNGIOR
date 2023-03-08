# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 13:29:17 2023

@author: danie
"""

import os
import sys
from pathlib import Path
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNetCV as EN
import scipy.stats as sts
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')
######add scripts to path####
experimentsPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries')
###################################

ipH_path = os.path.join(experimentsPath, 'bhbtri_ipH4.tsv')
metabolome = ['acetate', 'lactate', 'formate', 'butyrate', 'succinate']


y, x = [], []
with open(ipH_path) as f:
    header = f.readline().strip().split('\t')
    #get the pH index
    idx_ipH = header.index('pH')
    #indices of file header that match the metabolite list
    idx_metab = np.array([header.index(i.lower()) for i in metabolome if i.lower() in header])
    
    for line in f:
        a = line.strip().split('\t')
        x.append(np.array([float(a[z]) for z in idx_metab]))
        y.append(float(a[idx_ipH]))

x = np.array(x)
y = np.array(y)



X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.33)
    
#elastic net regression
model = EN(cv=5)
m = model.fit(X_train,y_train)

a = m.predict(X_test)

#plot
fig, ax = plt.subplots()
ax.set(xlabel='measured pH', ylabel='predicted pH', title='Quality of elasticNet pH predictions', xlim=(4.5,7.5),  ylim=(4.5,7.5))


sns.regplot(x=y_test, y = a, truncate=False)

props = dict(boxstyle='round', facecolor='#9ca3a6', alpha=0.5)

textstr = '\n'.join((
    "33% of the data (test-set):",
    "$r^2$ = {:.3f}".format(sts.pearsonr(y_test, a)[0]**2),
    "p = {:.2E}".format(sts.pearsonr(y_test, a)[1])))
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'pHregression.png'), dpi = 100)