# -*- coding: utf-8 -*-
"""
Created on Sun May  8 11:09:00 2022

@author: u0139894
"""

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


def update_subpopulations(conn, newTup):
    """
    newTup: mumax, pHoptimal,phAlpha 
    """
    sql = ''' UPDATE subpopulations
              SET mumax = ? ,
                  pHoptimal = ? ,
                  pHalpha = ?
              WHERE id = ?'''
    cur = conn.cursor()
    cur.execute(sql, newTup)
    conn.commit()


def update_subpopulations2subpopulations(conn, newTup):
    """
    newTup: hill, rate, id 
    """
    sql = ''' UPDATE subpopulations2subpopulations
              SET hillFunc = ?,
                  rate = ?
              WHERE id = ?'''
    cur = conn.cursor()
    cur.execute(sql, newTup)
    conn.commit()


def update_feedingTerms2metabolites(conn, newTup):
    """
    newTup: yield, monodK, id 
    """
    sql = ''' UPDATE feedingTerms2metabolites
              SET yield = ? ,
                  monodK = ? 
              WHERE id = ?'''
    cur = conn.cursor()
    cur.execute(sql, newTup)
    conn.commit()

def update_wc(conn, newTup):
    """
    newTup: concentratio, metabolite 
    """
    sql = ''' UPDATE wc
              SET concentration = ?
              WHERE metabolite = ?'''
    cur = conn.cursor()
    cur.execute(sql, newTup)
    conn.commit()


# database = 'C:\\Users\\u0139894\\Documents\\GitHub\\SyntheticCommunity\\scripts\\MODEL\\db\\parameterFit\\modelDB.sqlite3'
# conn = create_connection(database)
# with conn:
#         update_subpopulations(conn, ('.1', 7.0, 250, 'bt.lag'))
#         update_subpopulations2subpopulations(conn, ("(metObj.metD['glucose'].concentration > 0.1) or (metObj.metD['pyruvate'].concentration > .001)", 0.1, 1))
#         update_feedingTerms2metabolites(conn, (2.0, 4.5,1))