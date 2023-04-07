import os, sys


BASE = os.path.dirname(os.path.abspath(__file__))

FILES_PATH       = os.path.join(BASE, "files")
TRAINED_NN_PATH  = os.path.join(FILES_PATH, "NN")
TRAINED_NN_MSEED = os.path.join(TRAINED_NN_PATH, "NN_MS.h5")
TRAINED_NN_BIGG  = os.path.join(TRAINED_NN_PATH, "NN_BG.h5")

MODELS_PATH      = os.path.join(FILES_PATH, "models")

BIOCHEMISTRY_PATH   = os.path.join(FILES_PATH, "biochemistry")
MODELSEED_REACTIONS = os.path.join(BIOCHEMISTRY_PATH, "reactions.tsv")

