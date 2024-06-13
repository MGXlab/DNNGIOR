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

You may find examples of gap-filling a genome scale reconstruction (GEM) with `dnngior` with a complete or a defined medium in this [example notebook](tutorials/example.ipynb). `dnngior` can gapfill both ModelSEED and BiGG models, to gapfill BiGG models you need to specify modeltype. 

```python
Gapfill(path_to_BiGG_model, modeltype='BiGG')
```

## Custom Networks

By default `dnngior` uses an universally trained network capable of accurate predictions under most circumstances. If desired, it is possible to change the Neural Network you want to use during gapfilling:

```python
Gapfill(path_to_model, trainedNNPath=path_to_NN)
```

You can train your own Neural Network following this tutorial: [example training NN](tutorials/NN_training_example.ipynb).

Alternatively you can find additional custom Neural Networks for several taxonomic groups: [Custom Networks](docs/NN/custom_networks/). Upon request additional specially trained networks can be made available for specific biomes or taxonomic groups.


## License

Please see [License](LICENSE)


## Cite

The paper that will accompany the tool is currrently available as preprint:
https://www.biorxiv.org/content/10.1101/2023.07.10.548314v1


