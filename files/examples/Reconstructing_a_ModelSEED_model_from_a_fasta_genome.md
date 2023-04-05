# Reconstructing a ModelSEED metabolic model from a genome fasta file



This example requires [ModelSEEDpy](https://github.com/ModelSEED/ModelSEEDpy), which can be installed by the following:



```python
!git clone https://github.com/ModelSEED/ModelSEEDpy

!pip install ModelSEEDpy/.
```

Next we get a protein fasta file from an ensembl genome

```python
!wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-56/fasta/bacteria_119_collection/blautia_hydrogenotrophica_dsm_10507_gca_000157975/pep/Blautia_hydrogenotrophica_dsm_10507_gca_000157975.ASM15797v1.pep.all.fa.gz
!gzip -d Blautia_hydrogenotrophica_dsm_10507_gca_000157975.ASM15797v1.pep.all.fa.gz
```



And build a model for this file



```python
# Set the path to your genome
print("Build MSGenome object")
ensembl_genome = MSGenome.from_fasta('Blautia_hydrogenotrophica_dsm_10507_gca_000157975.ASM15797v1.pep.all.fa', split = ' ')

rast = RastClient()
rast.annotate_genome(ensembl_genome)

print("Start building base model")
base_model = MSBuilder.build_metabolic_model(model_id = "Blautia hydrogenotrophica DSM 10507", 
                                             genome   = ensembl_genome, 
                                             index    = "0",
                                             classic_biomass = True, 
                                             gapfill_model   = False, 
                                             gapfill_media   = None, 
                                             annotate_with_rast = False,
                                             allow_all_non_grp_reactions = True
                                            )


model_name = "bh_ungapfilled_model.sbml"
cobra.io.write_sbml_model(cobra_model = base_model, filename = os.path.join('DNNGIOR', 'files', 'models', model_name))
print("An ungapfilled modelSEED reconstructed model in now available.")
```

**note**: this example is only instructive and is not guaranteed by the authors to work since it depends on third party software.