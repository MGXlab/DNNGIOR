import dnngior
import os



base_path = "/".join(os.path.abspath(__file__).split("/")[:-3])

# Draft model
draftModel = os.path.join(base_path, "docs/models/e_coli_core.xml")


# Gapfill it
gf_object = dnngior.Gapfill(draftModel, objectiveName="BIOMASS_Ecoli_core_w_GAM")

