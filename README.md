# A novel way to gapfill metabolic models

## Installation instructions

To run the `dnngior` gapfiller, the [Gurobi](https://www.gurobi.com/) solver is required (mandatory).

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

You may find examples of gapfilling a genome scale reconstruction (GEM) with `dnngior` with a complete or a defined medium in this [example notebook](tutorials/example.ipynb).


## License




## Cite



