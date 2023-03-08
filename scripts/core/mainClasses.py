# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 12:55:23 2022

@author: u0139894
"""
import numpy as np
from scipy.special import gammaln
from scipy.integrate import solve_ivp as solver
from sklearn.linear_model import LinearRegression as LR
from sklearn.linear_model import ElasticNetCV as EN
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def getpH(metabolome, ipH_path):
    '''
    Parameters
    ----------
    metabolome : list
        list of metabolite names. For instance, the one retrieved from the `Metabolome.metabolites` class.
    ipH_path : str.
        path to file containing a column named 'pH' with pH values and columns with metabolite concentrations with names consistent with the names of the `metabolites` list. 
        Only columns for metabolites that influence the pH are needed.

    Returns
    -------
    func
        returns a function that can predict the pH from the `Metabolome` class.

    '''
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

    #simple linear regression model
    model = EN(cv=5)
    m = model.fit(x,y) 
    
    #indices of metabolites that match file header
    met_idx = np.array([i for i in range(len(metabolome)) if metabolome[i].lower() in  header])  

     
    
    def predictpH(metabolome_c):
        
        a = m.predict(metabolome_c[met_idx].reshape(1,-1))[0]
        a = max(3,a)
        a=min(10,a)
        
        
        return a
    return predictpH



#ipH_path = 'C:/Users/u0139894/Documents/GitHub/SyntheticCommunity/Tables/BTRI_ipH.txt'


class Metabolite:
    def __init__(self, name:str, concentration:float, formula:dict, color = ' #0093f5'):
        '''
        Parameters
        ----------
        name : str
            DESCRIPTION. metabolite name
        concentration : float
            DESCRIPTION. metabolite concentration
        formula : dict
            DESCRIPTION.metabolite elemental composition {'C':6, 'H':12, 'O':6}
        color : TYPE, optional
            DESCRIPTION. The default is  '#0093f5'.
            
        
        
        ```
        glucose = Metabolite('glucose', 8.0, {'C':6, 'H':12, 'O':6}, '#ff0000')
        glucose.add(-20)
        glucose.update(3)
        
        ```

        '''
        self.name = name
        
        self.formula = formula
        self.color = color
        
        self.carbons = 0
        self.hydrogens = 0
        
        if 'C' in self.formula:
            self.carbons = self.formula['C']
            
        if 'H' in self.formula:
            self.hydrogens = self.formula['H']
        
        self.concentration = None
        self.carbonMols = None
        self.hydrogenMols = None
        
        self.update(concentration) #in mM
    
    def update(self, concentration):
        self.concentration = concentration
        self.concentration = max(self.concentration,0)
        self.carbonMols = concentration*self.carbons
        self.hydrogenMols = concentration*self.hydrogens
    
    def add(self, concentration):
        self.concentration += concentration
        self.concentration = max(self.concentration,0)
        self.carbonMols = concentration*self.carbons
        self.hydrogenMols = concentration*self.hydrogens
    
    




class Metabolome:
    def __init__(self, metabolites : list, pH: float, pHFunc = None):
        '''
        Parameters
        ----------
        metabolites : list
            DESCRIPTION. list of `Metabolite` objects (created with the `Metabolite` class)
        pH : float
            DESCRIPTION.current value of the pH.
        pHFunc : TYPE, optional
            DESCRIPTION. A function describing the dependence of the pH on a vector of metabolite concentrations. This function can be generated with the `getpH` function. The default is None.

        Returns
        -------
        None.

        '''
        
        
        self.metD = {i.name: i for i in metabolites}
        self.metabolites = list(self.metD.keys())
        self.metabolites.sort()
        self.nmets = len(self.metabolites)
        self.pH = pH
        self.pHFunc = self.__getpHFunc(pHFunc)
    
        
        self.pH = self.pHFunc(self.get_concentration())
            
    
    def __getpHFunc(self, pHFunc):
        if pHFunc is None:
            def pHFunc(metC):
                return self.pH
            return pHFunc
        else:
            return pHFunc
            
        
    
    def get_concentration(self):
        self.concentration = np.array([self.metD[i].concentration for i in self.metabolites])
        return self.concentration
    
    def update(self, concentrationDict, add=True):
        
        if add:
            [self.metD[i].add(concentrationDict[i]) for i in concentrationDict]
            
        else:
            [self.metD[i].update(concentrationDict[i]) for i in concentrationDict]
        
        if self.pHFunc is not None:
            self.pH = self.pHFunc(self.get_concentration())
        
    

class FeedingTerm:
    def __init__(self, id : str, metDict : dict):
        '''
        Parameters
        ----------
        id : str
            DESCRIPTION. name of the functional term
        metDict : dict
            DESCRIPTION.[metabolite_id]:[yield, monodK]. All the metabolites in this dict will considered AND relationships (thus, multiplied)

        
        '''
        self.id = id
        self.metIDs = list(metDict.keys())
        self.metIDs.sort()
        self.yields = [metDict[metab][0] for metab in self.metIDs]
        self.monodKs = [metDict[metab][1] for metab in self.metIDs]
        self.intrinsicGrowth = self.__getIntrinsicGrowth()
        self.intrinsicMetabolism = self.__getIntrinsicMetabolism()
    
    def __getIntrinsicGrowth(self):
        
        def gr(metObj):
            
            metD = metObj.metD
            
            g = 1
            
            for i,v in enumerate(self.metIDs):
                
                if self.yields[i]>0:
                    g*= (metD[v].concentration/max((metD[v].concentration + self.monodKs[i]),0.0001))
            
            return g
        
        return gr
    
    def __getIntrinsicMetabolism(self):
        
        def metab(metObj):
            
            metD = metObj.metD
            
            omega = self.intrinsicGrowth(metObj)
            
            
            
            
            return -omega*np.array(self.yields)
            
        return metab
    
    
    

class Subpopulation:
    def __init__(self, name: str, count: float, species: str, mumax: float, feedingTerms : list, pHopt : float, pHalpha : float, state = 'active', color =  '#cf6f15'):
        '''
        Parameters
        ----------
        name : str
            DESCRIPTION.name of the subpopulation
        count : float
            DESCRIPTION.concentration of the subpopulation (in cells/10**5)
        species : str
            DESCRIPTION.species name
        mumax : float
            DESCRIPTION.maximum growth rate
        feedingTerms : list
            DESCRIPTION. list of feeding term objs
        pHopt : float
            DESCRIPTION. optimal pH
        pHalpha : float
            DESCRIPTION.pH sensitivity parameter
        state : str, optional
            DESCRIPTION. One of three state 'live', 'inactive' (pi positive), 'dead' (burst). The default is 'active'.
        color : TYPE, optional
            DESCRIPTION. The default is '#cf6f15'.

       
        '''
        
        
        
        self.name = name
        self.count = count
        self.species = species
        self.mumax = mumax
        self.state = state
        self.color = color
        
        self.feedingTerms = feedingTerms
        self.pHopt = pHopt
        self.pHalpha = pHalpha
        
        self.pHSensitivity = self.__getpHSensitivity()
        self.intrinsicGrowth = self.__getIntrGrowth()
        self.intrinsicMetabolism = self.__getIntrMetabolism()
    
    def __getpHSensitivity(self):
        self.pHbeta = (self.pHalpha-1)/self.pHopt
        
        maxima = self.gammaD(self.pHopt, self.pHalpha, self.pHbeta)
    
        def pHSensitivity(pH):
            return (1/maxima) * self.gammaD(pH, self.pHalpha, self.pHbeta)
        
        return pHSensitivity

    
    def __getIntrGrowth(self):
        def gr(metObj):
            growth = 0
            for fterm in self.feedingTerms:
                growth += fterm.intrinsicGrowth(metObj)
            
            return self.mumax*self.count*growth
            
        return gr
    
    def __getIntrMetabolism(self):
        def metabolism(metObj):
            
            metabV = np.zeros(metObj.nmets)
            
            for fterm in self.feedingTerms:
                metabV += fterm.intrinsicMetabolism(metObj)
            
            return self.mumax*self.count*metabV
            
        
        return metabolism
    

    @staticmethod
    def gammaD(x, alpha, beta):
        return np.exp(alpha*np.log(beta) - gammaln(alpha) + (alpha-1)*np.log(x) - beta*x)
        

    
        



class Bacteria:
    def __init__(self, species: str, subpopulations: dict, connections : dict, color = '#54f542'):
        '''
        Parameters
        ----------
        species : str
            DESCRIPTION. Species name
        subpopulations : dict
            DESCRIPTION. {nameOfSubpopulation: subpopulationObject}
        connections : dict
            DESCRIPTION. {subpopulationA: [[subpopulationB,
                               condition(metabolome)],
                               [subpopulationC,
                                condition(metabolome)]
                               
                               ]
                          }
        color : str, optional
            DESCRIPTION. The default is '#54f542'.

        

        '''
        
        self.species = species
        self.subpopulations = subpopulations
        self.connections = connections
        self.color = color
        
    
    def count(self):
        self.composition = {'active':0, 'inactive':0, 'dead':0}
        for subP in self.subpopulations:
            self.composition[self.subpopulations[subP].state]+=self.subpopulations[subP].count
        return self.composition
    
    def growth(self, metObj):
        growth = {subP:0 for subP in self.subpopulations}
        for subP in self.subpopulations:
            pop = self.subpopulations[subP]
            
            growth[subP] += pop.intrinsicGrowth(metObj) * pop.pHSensitivity(metObj.pH)
            
            for connection in self.connections[subP]:
                
                growth[connection[0]] += pop.count * connection[1](metObj) * connection[2]
                
                growth[subP] -= min(pop.count, pop.count * connection[1](metObj) * connection[2])
                
        return growth
    
    def metabolism(self, metObj):
        metV = np.zeros(metObj.nmets)
        
        for subP in self.subpopulations:
            pop = self.subpopulations[subP]
            metV+= pop.intrinsicMetabolism(metObj) * pop.pHSensitivity(metObj.pH)
        
        return metV
        
            
        




class Microbiome:
    def __init__(self, bacteria : dict):
        self.bacteria = bacteria
        self.species = list(self.bacteria.keys())
        self.subpopD = {}
        for bac in self.bacteria:
            self.subpopD = self.subpopD | self.bacteria[bac].subpopulations
        self.subpops = list(self.subpopD.keys())
        self.subpops.sort()
        self.nsubpops = len(self.subpops)
        self.nspecies = len(self.species)
        
    
    def growth(self, metObj):
        g = {}
        
        for bac in self.bacteria:
            b = self.bacteria[bac]
            g = g | b.growth(metObj)
        
        return g
    
    def metabolism(self, metObj):
        metabV = np.zeros(metObj.nmets)
        
        for bac in self.bacteria:
            metabV+= self.bacteria[bac].metabolism(metObj)
        
        return metabV
    
    def count(self):
        self.counts = {bac:self.bacteria[bac].count() for bac in self.bacteria}
        
        return self.counts
    
    def countSubpops(self):
        return np.array([max(0,self.subpopD[sp].count) for sp in self.subpops])
        
        

class Pulse:
    def __init__(self, metabolome, microbiome, t_start, t_end, n_steps, vin, vout, qin, qout):
        self.metabolome = metabolome
        self.microbiome = microbiome
        
        self.t_start = t_start
        self.t_end = t_end
        self.n_steps = n_steps
        self.range = np.linspace(self.t_start, self.t_end, self.n_steps)
        self.vin = vin
        self.vout = vout
        self.qin = qin
        self.qout = qout
        
        
        
    

class Reactor:
    def __init__(self, microbiome, metabolome, pulses : list, volume):
        
        self.microbiome = microbiome
        self.metabolome = metabolome
        self.pulses = pulses
        self.volume = volume
        self.pH = self.metabolome.pH
        self.nstates = 1 + self.metabolome.nmets + self.microbiome.nsubpops
    
    
    
    def get_states(self):
        vec = np.zeros(self.nstates)
        
        vec[0] = self.volume
        
        vec[1 : 1 + self.metabolome.nmets] = self.metabolome.get_concentration()
        vec[1 + self.metabolome.nmets::] = self.microbiome.countSubpops()
        vec = np.round(vec, 5)
        return vec
    
    def update_states(self, vec):
        self.volume = vec[0]
        
        
        self.volume = max(self.volume, 0)
        
        
        
        for idx,met in enumerate(self.metabolome.metabolites):
            self.metabolome.metD[met].update(vec[idx + 1])
        
        self.pH = self.metabolome.pHFunc(self.metabolome.get_concentration())
        self.metabolome.pH = self.pH
        for idx, subP in enumerate(self.microbiome.subpops):
            self.microbiome.subpopD[subP].count = vec[1 + self.metabolome.nmets + idx]
        
        
    
    def dvdt(self, pulseObj):
        
        return max(-self.volume, (pulseObj.qin-pulseObj.qout))
    
    
    
    def dsdt(self, metObj, pulseObj):
        #change due to influx
        #change due to metabolism
        
        current_concentration = self.metabolome.get_concentration()
        dsdt =  (pulseObj.qin/self.volume) * (pulseObj.metabolome.get_concentration()-current_concentration) + self.microbiome.metabolism(metObj)
        return np.maximum(-current_concentration, dsdt)
    
    def dxdt(self, metObj, pulseObj):
        sp_dxdt = self.microbiome.growth(metObj)
        #change due to influx
        #change due to growth
        current_population = self.microbiome.countSubpops()
        dxdt = (pulseObj.qin/self.volume) * (pulseObj.microbiome.countSubpops() - current_population) + np.array([sp_dxdt[sp] for sp in self.microbiome.subpops]) 
        
        return np.maximum(-current_population, dxdt)
    
    
    def ode(self, t, states, pulse):
        #1) update states
        
        self.update_states(states)
        
        
        derivatives = np.zeros(self.nstates)
        
        #volume
        derivatives[0] = self.dvdt(pulse)
        #pH
        
        #metabolites
        derivatives[1 : 1 + self.metabolome.nmets] = self.dsdt(self.metabolome, pulse)
        #subpopulations
        derivatives[1 + self.metabolome.nmets ::] = self.dxdt(self.metabolome, pulse)
        
        
        
        
        return derivatives
    
    def simulate(self):
        #1) perform the initial steps of a pulse (non-continuous)
        #2) integrate the continuous steps of the pulse
        # Solve ODE
        
        ts = np.empty((0))
        
        vol_dyn = np.empty((0))
        met_dyn = []
        subpop_dyn = []
        
        for pulse in self.pulses:
                
            #change in volume
            
            state = np.zeros(self.nstates)
            
            state[0] = max(0, self.volume + pulse.vin - pulse.vout)
            
            
            metReactor = self.metabolome.get_concentration()
            metPulse = pulse.metabolome.get_concentration()
            
            state[1 : 1 + self.metabolome.nmets] = ((self.volume*metReactor) + (metPulse*pulse.vin) - (metReactor*pulse.vout))/state[0]
            
            subpopReactor = self.microbiome.countSubpops()
            subpopPulse = pulse.microbiome.countSubpops()
            
            
            state[1 + self.metabolome.nmets::] = ((self.volume*subpopReactor) + (subpopPulse*pulse.vin) - (subpopReactor*pulse.vout))/state[0]
            
            self.update_states(state)
            
            
            
            # ode = solver.ode(self.ode)
            
            # ode.set_f_params(pulse)
            
            # ode.set_integrator('vode', nsteps=pulse.n_steps, method= 'bdf')
            
            # y_init = self.get_states()
            
            # # Time
           
            # t_step = (pulse.t_end - pulse.t_start)/pulse.n_steps
            
            # ode.set_initial_value(y_init,pulse.t_start)
            
            self.solution = solver(fun=self.ode, t_span =(pulse.t_start, pulse.t_end), y0 = self.get_states(), t_eval = pulse.range, args=[pulse])
            
            ts = np.concatenate([ts,self.solution.t])
            vol_dyn = np.concatenate([vol_dyn, self.solution.y[0]])
            
            met_dyn.append(self.solution.y[1 : 1 + self.metabolome.nmets])
            subpop_dyn.append(self.solution.y[1 + self.metabolome.nmets::])
            
            
            
            
            
        #     while ode.successful() and ode.t < pulse.t_end:
                
        #         ode.integrate(ode.t + t_step)
        #         ts.append(ode.t)
        #         pH_dyn.append(self.metabolome.pH)
        #         vol_dyn.append(self.volume)
        #         met_dyn.append(self.metabolome.get_concentration())
        #         subpop_dyn.append(self.microbiome.countSubpops())
        #         cell_dyn.append(self.microbiome.count())
           
        # #integration finished, store the results
        self.time_simul = np.array(ts)
        
        self.v_simul = np.array(vol_dyn)
        self.met_simul = np.hstack(met_dyn)
        self.pH_simul = np.array([self.metabolome.pHFunc(mc) for mc in self.met_simul.T])
        self.subpop_simul = np.hstack(subpop_dyn)
        
        cellActive_dyn = {i: np.zeros(len(self.subpop_simul[0])) for i in self.microbiome.species}
        cellInactive_dyn = {i: np.zeros(len(self.subpop_simul[0])) for i in self.microbiome.species}       
        cellDead_dyn = {i: np.zeros(len(self.subpop_simul[0])) for i in self.microbiome.species}
        
        for i,v in enumerate(self.microbiome.subpops):
            if self.microbiome.subpopD[v].state == 'active':
                cellActive_dyn[self.microbiome.subpopD[v].species] += self.subpop_simul[i]
            elif self.microbiome.subpopD[v].state == 'inactive':
                cellInactive_dyn[self.microbiome.subpopD[v].species] += self.subpop_simul[i]
            elif self.microbiome.subpopD[v].state == 'dead':
                cellDead_dyn[self.microbiome.subpopD[v].species] += self.subpop_simul[i]
        
        self.cellActive_dyn = np.array([cellActive_dyn[i] for i in self.microbiome.species])
        self.cellInactive_dyn = np.array([cellInactive_dyn[i] for i in self.microbiome.species])
        self.cellDead_dyn = np.array([cellDead_dyn[i] for i in self.microbiome.species])
                
                
        
     
        
    
    def makePlots(self, path=None):
        
        fig = make_subplots(rows = 2, cols = 2, subplot_titles = ['Metabolites', 'pH', 'Subpopulations', 'State'])
        
        for i,v in enumerate(self.metabolome.metabolites):
            fig.add_trace(go.Scatter(x=self.time_simul, y = self.met_simul[i], mode='lines', name=v, line=dict(color = self.metabolome.metD[v].color, simplify=True), opacity=0.8), row=1, col=1)
        
        fig.add_trace(go.Scatter(x=self.time_simul, y = self.pH_simul, mode='lines', name='pH', line=dict(color = 'rgb(57,255,20)', simplify=True), opacity=0.8), row=1, col=2)
        
        
        for i,v in enumerate(self.microbiome.subpops):
            if self.microbiome.subpopD[v].state=='active':
                opacity = 0.8
            else:
                opacity = 0.1
            fig.add_trace(go.Scatter(x=self.time_simul, y = self.subpop_simul[i], mode='lines', name=v, line=dict(color = self.microbiome.subpopD[v].color, simplify=True), opacity=opacity), row=2, col=1)
        
        
        for i,v in enumerate(self.microbiome.species):
            fig.add_trace(go.Scatter(x=self.time_simul, y = self.cellActive_dyn[i], mode='lines', name=v + '_active', line=dict(color = self.microbiome.bacteria[v].color, width = 2, simplify=True), opacity=1), row=2, col=2)
            fig.add_trace(go.Scatter(x=self.time_simul, y = self.cellInactive_dyn[i], mode='lines', name=v + '_inactive', line=dict(color = self.microbiome.bacteria[v].color, width = 2, simplify=True), opacity=.1), row=2, col=2)
            fig.add_trace(go.Scatter(x=self.time_simul, y = self.cellDead_dyn[i], mode='lines', name=v + '_dead', line=dict(color = self.microbiome.bacteria[v].color, width = 2, simplify=True), opacity=.1), row=2, col=2)
        if path is not None:
            fig.write_html(path)
        
        fig.show()
                    
        


        
        #metabolome_dynamics
        #active cells / inactive cells
        #pH dynamics
        #subpopulation dynamics
        
        
        
