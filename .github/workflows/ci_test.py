import dnngior
import os

base_path = "/".join(os.path.abspath(__file__).split("/")[:-3])

# BiGG case
print("gapfill using a BiGG model")

# Draft model
draftModel = os.path.join(base_path, "docs/models/e_coli_core.xml")

# Gapfill it
gf_object = dnngior.Gapfill(draftModel, dbType="BiGG", objectiveName="BIOMASS_Ecoli_core_w_GAM")

# ModelSEED case
print("gapfill using a ModelSEED model")

# Draft model
draftModel = os.path.join(base_path, "docs/models/e_coli_core_Seed.sbml")

# Gapfill it
gf_object = dnngior.Gapfill(draftModel, objectiveName="bio1")
