__version__ = "0.0.3"

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import tarfile
from pathlib import Path

path        = Path(os.path.dirname(__file__))
BiochemRoot = os.path.join(path, 'files', 'biochemistry')

rxns_file = "reactions.tsv"
ff = Path(os.path.join(BiochemRoot, rxns_file))
if not ff.is_file():
    file = tarfile.open(os.path.join(BiochemRoot, 'reactions.tsv.tar.gz'))
    file.extractall(BiochemRoot)

cpds_file = "compounds.tsv"
ff = Path(os.path.join(BiochemRoot, cpds_file))
if not ff.is_file():
    file = tarfile.open(os.path.join(BiochemRoot, 'compounds.tsv.tar.gz'))
    file.extractall(BiochemRoot)

import dnngior.files

from dnngior.variables import *

from dnngior.MSEED_compounds import Compounds
from dnngior.MSEED_reactions import Reactions
from dnngior.reaction_class  import Reaction

import dnngior.addExchanges 

from dnngior.gapfill_class import Gapfill

from dnngior.build_model import refine_model

from dnngior.NN_Predictor import NN
