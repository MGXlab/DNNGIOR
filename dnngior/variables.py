import os, sys


BASE = os.path.dirname(os.path.abspath(__file__))

FILES_PATH         = os.path.join(BASE, "files")
DEFAULT_TRAINED_NN = os.path.join(FILES_PATH, "NN")
MODELS_PATH        = os.path.join(FILES_PATH, "models")

BIOCHEMISTRY_PATH   = os.path.join(FILES_PATH, "biochemistry")
MODELSEED_REACTIONS = os.path.join(BIOCHEMISTRY_PATH, "reactions.tsv")

# sys.path.append(BASE)

