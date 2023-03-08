# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 09:49:39 2022

@author: u0139894
"""

import os
from pathlib import Path
import sqlite3

#Replace this name for the desired name.
databaseName = 'modelDB_bhD.sqlite3'

#by default, the sqlite3 file is stored in /files/dbs
databaseFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

con = sqlite3.connect(os.path.join(databaseFolder, databaseName))

cur = con.cursor()


elements = """
CREATE TABLE IF NOT EXISTS  elements (
id TEXT PRIMARY KEY NOT NULL UNIQUE,
name TEXT NOT NULL UNIQUE,
MolecularWeight REAL NOT NULL)"""

cur.execute(elements)

metabolites = """
CREATE TABLE IF NOT EXISTS  metabolites (
id TEXT PRIMARY KEY NOT NULL UNIQUE,
color TEXT NOT NULL,
MolecularWeight REAL NOT NULL)"""
cur.execute(metabolites)

metabolites2elements = """
CREATE TABLE IF NOT EXISTS  metabolites2elements (
id INTEGER PRIMARY KEY AUTOINCREMENT,
metabolite TEXT NOT NULL,
element TEXT NOT NULL,
atoms INTEGER NOT NULL,
FOREIGN KEY (metabolite) REFERENCES metabolites (id),
FOREIGN KEY (element) REFERENCES elements (id))
"""
cur.execute(metabolites2elements)

wc = """
CREATE TABLE IF NOT EXISTS  wc (
metabolite TEXT PRIMARY KEY UNIQUE,
concentration REAL NOT NULL,
FOREIGN KEY (metabolite) REFERENCES metabolites (id))
"""
cur.execute(wc)


species = """
CREATE TABLE IF NOT EXISTS  species (
id TEXT PRIMARY KEY UNIQUE,
name TEXT NOT NULL,
genomeSize INTEGER,
geneNumber INTEGER,
patricID TEXT,
ncbiID TEXT)
"""
cur.execute(species)

feedingTerms = """
CREATE TABLE IF NOT EXISTS  feedingTerms (
id TEXT PRIMARY KEY UNIQUE,
name TEXT NOT NULL,
species TEXT NOT NULL,
FOREIGN KEY (species) REFERENCES species (id))
"""
cur.execute(feedingTerms)

feedingTerms2metabolites = """
CREATE TABLE IF NOT EXISTS  feedingTerms2metabolites (
id INTEGER PRIMARY KEY AUTOINCREMENT,
feedingTerm TEXT NOT NULL,
metabolite TEXT NOT NULL,
yield REAL NOT NULL,
monodK REAL NOT NULL,
FOREIGN KEY (feedingTerm) REFERENCES feedingTerms (id),
FOREIGN KEY (metabolite) REFERENCES metabolites (id))
"""
cur.execute(feedingTerms2metabolites)

subpopulations = """
CREATE TABLE IF NOT EXISTS  subpopulations (
id TEXT PRIMARY KEY UNIQUE,
species TEXT NOT NULL,
mumax REAL NOT NULL,
pHoptimal REAL NOT NULL,
pHalpha REAL NOT NULL,
count REAL NOT NULL,
color TEXT NOT NULL,
state TEXT NOT NULL,
FOREIGN KEY (species) REFERENCES species (id))
"""
cur.execute(subpopulations)

subpopulations2feedingTerms = """
CREATE TABLE IF NOT EXISTS  subpopulations2feedingTerms (
id INTEGER PRIMARY KEY AUTOINCREMENT,
subpopulation TEXT NOT NULL,
feedingTerm TEXT NOT NULL,
FOREIGN KEY (subpopulation) REFERENCES subpopulations (id),
FOREIGN KEY (feedingTerm) REFERENCES feedingTerms (id))
"""
cur.execute(subpopulations2feedingTerms)

subpopulations2subpopulations = """
CREATE TABLE IF NOT EXISTS  subpopulations2subpopulations (
id INTEGER PRIMARY KEY AUTOINCREMENT,
subpopulation_A TEXT NOT NULL,
subpopulation_B TEXT NOT NULL,
hillFunc TEXT NOT NULL,
rate REAL NOT NULL,
FOREIGN KEY (subpopulation_A) REFERENCES subpopulations (id),
FOREIGN KEY (subpopulation_B) REFERENCES subpopulations (id))
"""
cur.execute(subpopulations2subpopulations)


con.commit()