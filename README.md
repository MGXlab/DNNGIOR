# A novel way to gapfill metabolic models

## Installation instructions

To run the `dnngior` gapfiller, the [Gurobi](https://www.gurobi.com/) solver is mandatory.

```bash
pip install gurobipy
```

To use `gurobi`, you need a [license](https://www.gurobi.com/downloads/). If you are an acedemic, you may get a license for free.

Once you have successfully installed `gurobi`, you are ready to install the `dnngior` gapfiller.

```bash
pip install dnngior
```

Optionally, you may need to also get [Tensorflow](https://www.tensorflow.org/install) (or through [conda](https://anaconda.org/conda-forge/tensorflow)) 
in case you would like to use the `NN_Trainer`.

## How to use

Gapfilling models is done using the `Gapfill` class:
```python
import dnngior.gapfill_class.Gapfill  
Gapfill(path_to_model)
```

You may find examples of gap-filling a genome scale reconstruction (GEM) with `dnngior` with a complete or a defined medium in this [example notebook](tutorials/gapfilling_example.ipynb). `dnngior` can gapfill both ModelSEED and BiGG models, to gapfill BiGG models you need to specify modeltype. 

```python
Gapfill(path_to_BiGG_model, dbtype='BiGG')
```

## Custom Networks

By default `dnngior` uses an universally trained network capable of accurate predictions under most circumstances. If desired, it is possible to change the Neural Network you want to use during gapfilling:

```python
Gapfill(path_to_model, trainedNNPath=path_to_NN)
```

You can train your own Neural Network following this tutorial: [example training NN](tutorials/NN_training_example.ipynb).

Alternatively you can find additional custom Neural Networks for several taxonomic groups: [Custom Networks](docs/NN/custom_networks/). Upon request additional specially trained networks can be made available for specific biomes or taxonomic groups.

## Batch CLI

You can use this [CLI script](dnngior/fasta2model_CLI.py) to automate the workflow of creating metabolic models from protein sequence data assuming default settings. This allows easier generation of larger number of models.

You can run it in two modes:

Option 1: From FASTA files
`python fasta2model_CLI.py -f /path/to/fasta_folder -o /path/to/output [--medium /path/to/media_file]`

The script will create MSGenome objects from each .faa file, optionally annotate them via RAST, build base metabolic models, and perform gap-filling based on complete medium or your specified medium if provided. 

Option 2: From existing base models
`python fasta2model_CLI.py -m /path/to/base_models -o /path/to/output [--medium /path/to/media_file]`

Skip base model construction and only perform gap-filling on pre-existing .xml model files.

Base models are saved to output_folder/base_models/, gap-filled models to output_folder/gapfilled_models/, and gap-filling statistics are logged in a gf_data.tsv file showing reaction counts per model. The optional -e flag lets you specify growth conditions via a medium file.

See also `python fasta2model_CLI.py -h`

## License

Please see [License](LICENSE)


## Cite

The paper that accompanies the tool can be found here:\
https://www.cell.com/iscience/fulltext/S2589-0042(24)02574-4


