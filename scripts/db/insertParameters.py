# -*- coding: utf-8 -*-
"""
Created on Sun May  8 11:09:00 2022

@author: u0139894
"""


import os
from pathlib import Path
import sqlite3
from sqlite3 import Error


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn


def getTypes(conn, tableName):
    q = "PRAGMA table_info(" + tableName + ");"
    record = q
    
    cur = conn.cursor()
    cur.execute(q)
    return cur.fetchall()
    



def insert2db(conn, tableName, parameterIDs, *params):
    q = ['?']*len(params)
    
    tps = getTypes(conn, tableName)
    
    values2add = []
    
    for i,v in enumerate(params):
        if tps[i][2] == 'REAL':
            values2add.append(float(v))
        elif tps[i][2] == 'INTEGER':
            values2add.append(int(v))
        else:
            values2add.append(v)
    
    query = """INSERT INTO """ + tableName + """ (""" + ','.join(parameterIDs) + """) VALUES (""" + ','.join(q) + """)"""
    
    record = values2add
    
    cur = conn.cursor()
    cur.execute(query, record)
    conn.commit()




def addFile(conn, filePath, tableName):
    
    with open(filePath) as f:
        params = f.readline().strip().split('\t')
        
        
        for line in f:
            a = line.strip().split('\t')
            insert2db(conn, tableName, params, *a)




databaseName = 'modelDB_bhD.sqlite3'

databaseFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

conn = sqlite3.connect(os.path.join(databaseFolder, databaseName))

dbParamsFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'dbsTemplateTables' )




addFile(conn, os.path.join(dbParamsFolder, 'elements.tsv'), 'elements')
addFile(conn, os.path.join(dbParamsFolder, 'metabolites.tsv'), 'metabolites')
addFile(conn, os.path.join(dbParamsFolder, 'metabolites2elements.tsv'), 'metabolites2elements')
addFile(conn, os.path.join(dbParamsFolder, 'wc.tsv'), 'wc')
addFile(conn, os.path.join(dbParamsFolder, 'species.tsv'), 'species')
addFile(conn, os.path.join(dbParamsFolder, 'feedingTerms.tsv'), 'feedingTerms')
addFile(conn, os.path.join(dbParamsFolder, 'feedingTerms2metabolites.tsv'), 'feedingTerms2metabolites')
addFile(conn, os.path.join(dbParamsFolder, 'subpopulations.tsv'), 'subpopulations')
addFile(conn, os.path.join(dbParamsFolder, 'subpopulations2subpopulations.tsv'), 'subpopulations2subpopulations')
addFile(conn, os.path.join(dbParamsFolder, 'subpopulations2feedingTerms.tsv'), 'subpopulations2feedingTerms')



    

