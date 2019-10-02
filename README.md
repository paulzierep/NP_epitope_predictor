# NP_epitope_predictor

## Installation

### Clone the Repo

```
git clone https://github.com/paulzierep/NP_epitope_predictor.git
```

### Requirements (Conda)

Install conda:

[Conda installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

Activate a new environment:

```
conda activate **env_name**
```

Instal all packages:

```
conda install --file requirements.txt
```

### Requirements (Pip)

Unfortunatly rdkit is difficult to set up outside the coda env. Therefore this is not tested.
In general it should work to install rdkit with any method described here:

[rdkit installation](http://www.rdkit.org/docs/Install.html)

Then install the other packages with 

```
pip install -r requirements.txt
```

### Unzip the Data

In order to run a prediction it is not necessary to fit the predictor, it can run with pre-fitted data, which is provided in:

```
./tools/ML_data.zip
```

Unzip in same folder:

```
cd ./tools/
unzip ML_data.zip
```

## Run a prediction

An example how to run the prediction is given in the script:

```
cd ./tools/run_prediction.py
```

The predictor needs as input a correct SMILES string eg: "CCCCCCO".
If the smiles cannot be parsed the predictor returns None.

The predictor returns a dict which holds all needed information, eg.: Class, B_cell probability...

The utility functions: *Results_To_Json* creates a JSON object from the results dict. 
The json tree can be observed for example by copy-and-pasteing into http://jsonviewer.stack.hu/

This json can be tranvered into any desired html output.
