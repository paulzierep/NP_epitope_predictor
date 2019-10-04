
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

Install all packages:

```
conda install --file requirements.txt
```

### Requirements (Pip)

Unfortunately rdkit is difficult to set up outside the coda env. Therefore this is not tested.
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

### Input

An example how to run the prediction is given in the script:

```
cd ./tools/run_prediction.py
```

The predictor needs as input a correct SMILES string eg: "CCCCCCO".

There are two optional arguments which are currently included:

Argument | Parameter | Explanation
------------ | ----------|---
only_epitopes | True | The most similar epitopes are listed, 
only_epitopes | False | The most similar compounds from the entire ChEBI are listed
sort_order | "E" | Uses the eucleadien distance as similarity measure for the target compounds meaning the shown compounds have the most FPs in common with the query molecule 
sort_order | "T" | Uses the overall Tanimoto Coefficient as similarity measure, meaning the shown compound are most similar to the query molecule, this requires, that the similarity is computed between all compounds in the cluster and the query (takes a bit longer)

### Output

The predictor returns a dict which holds all needed information, eg.: Class, B_cell probability...

A tree of the dict is shown below. The ontology and fp_imp (FP importance) are pandas data frames. They can be converted into html, dict or json with (.to_html(),.to_dict() or to_json()).

```
error <class 'NoneType'>
input_svg <class 'str'>
cluster_info <class 'dict'>
	 Cluster <class 'int'>
	 Members <class 'int'>
	 Warning <class 'float'>
	 Name <class 'str'>
	 Ontology <class 'pandas.core.frame.DataFrame'>
classify_info_b <class 'dict'>
	 proba <class 'float'>
	 fp_imp <class 'pandas.core.frame.DataFrame'>
	 fitting_info <class 'NoneType'>
	 p_support <class 'int'>
	 n_support <class 'int'>
classify_info_t <class 'dict'>
	 proba <class 'float'>
	 fp_imp <class 'pandas.core.frame.DataFrame'>
	 fitting_info <class 'NoneType'>
	 p_support <class 'int'>
	 n_support <class 'int'>
```

The utility functions: *Results_To_Json* creates a JSON object from the results dict. 
The json tree can be observed for example by copy-and-pasting into http://jsonviewer.stack.hu/

The utility functions: *Results_To_Html* creates one possibility of how the HTML output could look like.
Although the incorporation of the json into a framework would be preferable. 

### Error handling

If the SMILES can not be parsed,
the parsing error is recorded and stored in the 
result dict in "error".
This can be shown to the user, as it is currently implemented in *Results_To_Json* and *Results_To_Html*.

## Training a new classifier

TODO

# Methods 

## Dataset

The entire ChEBI was assigned as negative samples (TODO explain why negative). ChEBI has curated and automatically 
assigned molecular samples. Only the curated samples were taken, marked with 3 stars in the ChEBI database.
Positive samples for B cell and T cell assays were downloaded from the IEDB via web-interface query (https://www.iedb.org/).

The positive samples where then marked as such in the ChEBI dataset.

All samples were parsed using RDKit and duplicated smiles were removed.

Final dataset:
ChEBI negative (42643), positive T Cell (579) and positive B Cell (2140).

## Molecular Fingerprints 

The molecules were encoded into vectors using Morgan fingerprints.

For the clustering folded bit fingerprints were computed. 

For the classification unfolded count fingerprints were used.
The unfolded count fingerprints were computed for each cluster separately.
As unfolded fingerprints lead to large vectors, only those fingerprints which occur
in at least 10 molecules were taken.

## Clustering

The molecules were clustered using k means (scikit-learn K-means) clustering, 
the optimal number of clusters was 
determined using the shoulder method. 
The clusters were described using BiNChE ontology enrichment analysis for each cluster.

## Classification

For each cluster a model was compiled which predicts the 
probability of a molecule to be a B cell or T cell epitope.

The models were created using random forests with 100 iterations (scikit-learn random forest classifier).

## Benchmark

Different sets of fingerprints were selected using a chi-squared feature selection approach
(scikit-learn SelectKBest).

The classifiers with each feature set were benchmarked using repeated (3) k-fold (5) cross validation.

The classifiers were compared against a classifier which uses a similarity based prediction approach.
This similarity classifiers were designed as follows: 
For each sample the similarity to all known epitopes 
in a cluster was computed using unfolded count fingerprints and tanimoto similarity.
The highest similarity score was then assigned as probability to this sample.

## Significant Fingerprints

Highly correlating features were removed. All samples that exceeded the person correlation coefficient (PCC)
with more then 0.8.

In order to highlight the fingerprints which are most significant to detect epitopes
the fingerprints distribution was investigated using the Mann–Whitney U test.

The null hypothesis was tested, that the fingerprint count from the epitope samples was
selected from the same distribution then all the ChEBI samples (without epitopes).

The p-values were Bonferroni corrected to account for the multiple comparison problem.
