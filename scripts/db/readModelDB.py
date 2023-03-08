# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:20:18 2022

@author: u0139894
"""


from pony.orm import *
from pathlib import Path
import sys
import os

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))
from mainClasses import *






######################################

def getTransitionFunction(statement):
    
    def tf(metObj):
        if statement == '""':
            return True
        return eval(statement)
    return tf



######################################

def get_database(filename='db/modelDB.sqlite3'):
    db = Database()
    db.bind(provider='sqlite', filename=filename)
    return db


@db_session()
def query_string(db, tableName, columnName, query):
    q = db.select("select * from " + tableName + " where " + columnName + " like " + "'%" +  query + "%'")
    
    return q

@db_session()
def query_value(db, tableName, columnName, query):
    q = db.select("select * from " + tableName + " where " + columnName + "='" + query + "'")
    return q



@db_session()
def createMetabolite(db, metName, concentration):
    formula = {i[2]:i[3] for i in query_value(db, 'metabolites2elements', 'metabolite', metName)}
      
    return Metabolite(name = metName, concentration = concentration, formula = formula, color = query_value(db, 'metabolites', 'id', metName)[0][1])



@db_session()
def createMetabolome(db, mediaName, pH0=6.4, pHFunc=None):
    metObjs = []
    metabolites = db.select("select * from " + mediaName) #name concentration
    
    for met in metabolites:
        metObjs.append(createMetabolite(db, met[0], met[1]))
    
    return Metabolome(metabolites = metObjs, pH = pH0, pHFunc = pHFunc)



@db_session()
def createFeedingTerm(db, feedingTermID, mediaName):
    metabolites = query_string(db, 'feedingTerms2metabolites', 'feedingTerm', feedingTermID)
    metabolome = createMetabolome(db, mediaName)
    term = {i[2]:(i[3], i[4]) for i in metabolites}
    for i in metabolome.metabolites:
        if i not in term:
            term[i] = (0,0)
    
    return FeedingTerm(id = feedingTermID, metDict = term)




@db_session()
def createSubpopulation(db, supopulationID, mediaName):
    spParams = query_value(db, 'subpopulations', 'id', supopulationID)[0]
    spD = [createFeedingTerm(db, i[2], mediaName) for i in query_value(db, 'subpopulations2feedingTerms', 'subpopulation', supopulationID)]
         
    return Subpopulation(name = spParams[0], count = spParams[5], species = spParams[1], mumax = spParams[2], feedingTerms = spD, pHopt = spParams[3], pHalpha = spParams[4], color=spParams[6], state = spParams[7])
    


@db_session()
def createBacteria(db, speciesID, mediaName):
    colors = {'bt':'#ff8300', 'ri':'#00b8ff', 'bh':'#FF10F0'}
    
    subpops = [createSubpopulation(db, i[0], mediaName) for i in query_value(db, 'subpopulations', 'species', speciesID)]
    spopsD = {i.name:i for i in subpops}
    
    transitions = {}
    for i in spopsD:
        transition = query_value(db, "subpopulations2subpopulations", "subpopulation_A", i)
        if len(transition)>0:
            relations = []    
            for z in transition:
                relations.append((z[2], getTransitionFunction(z[3]), z[4]))
            transitions[i] = relations
        else:
            transitions[i] = []
    return Bacteria(speciesID, spopsD, transitions, color = colors[speciesID])


# import os

# databaseName = 'modelDB.sqlite3'

# databaseFolder = 'C:\\Users\\danie\\Documents\\GitHub\\SyntheticCommunity\\scripts\\MODEL_hill\\db'

# db = get_database(os.path.join(databaseFolder, databaseName))

