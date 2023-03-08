# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 20:22:22 2023

@author: danie
"""

import os
import sys
from pathlib import Path
import numpy as np
from scipy.special import gammaln
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')


def pHS(pH, phOpt, pHalpha):
    pHbeta = (pHalpha-1)/phOpt
    
    maxima = gammaD(phOpt, pHalpha, pHbeta)


    return (1/maxima)*gammaD(pH, pHalpha, pHbeta)


def gammaD(pH, pHalpha, pHbeta):
    
    return np.exp(pHalpha*np.log(pHbeta) - gammaln(pHalpha) + (pHalpha-1)*np.log(pH) - pHbeta*pH)


pHvals = np.linspace(4.5,9.5,1000)

alphas = np.linspace(5,1000,1000)

idx = np.arange(0, 1000, 100)

fig, ax = plt.subplots()
ax.set(xlabel='pH', ylabel='subpopulation pH sensitivity', title='Impact of pH on growth')

col = np.linspace(0,1, len(alphas))


for i,v in enumerate(alphas):
    sensitivity = np.array([1-pHS(z, 7, v) for z in pHvals])
    if i in idx:
        
        ax.plot(pHvals, sensitivity, label = "{:.1f}".format(v), lw=
            3, color=plt.cm.winter(col[i]))
        
    else:
        ax.plot(pHvals, sensitivity, lw=.2, color=plt.cm.winter(col[i]))

ax.legend()


plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'pHsensitivity.png'), dpi = 100)