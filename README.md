# Molecular feature and IE&EA prediction
This repository contains code for converting molecule stucture into four different representations(SYBYL, ECFP, ECFP&SYBYL, ECFP_num). Additionally, corresponding model had been trained by Material project database. These models can be used to predict ionization energy (IE) and electron affinity (EA) from molecular structure.

## Evironment Setting
### via Anaconda
Create a conda environment:

```
conda env create -f environment.yml
source activate molecule_feature_predict
```



## Requirements 
* python
* rdkit
* numpy
* openbabel
* argparse
* pandas
* os
* biopandas
* keras >=2.00 && <=2.0.7
* tensorflow==1.10.0
* scikit-learn

## Tutorial
Prepare a `.csv` file (for example, `test.csv`) which contains molecular SMILES, and using the command to create different molecular representations(`SYBYL.csv, ECFP.csv, ECFP_SYBYL.csv, ECFP_num.csv`):

```python featurize.py -i test.csv```

>if we want to change ECFP's the number of bits to 1024 (default = 2048) or/and the radius to 3 (default = 2), we can use the command:

>```python featurize.py -i test.csv -n 1024 -r 3```


After converting molecular structure into representations, we can predict IE and EA by using models which have been trained by Material Project database. All the models are stored in the folder `model`.

For example, if we want to use ECFP as representation to predict IE and EA, we can use the command:

```python predict.py -i ECFP.csv -m ECFP```


and it will generate  `ECFP_predict.csv` file,  which includes SMILES and IE and EA which are predicted by a model.

Here, we provide four different models :`ECFP, SYBYL, ECFP_SYBYL, ECFP_num`

so we can use these four different commands to predict IE and EA

>ECFP : ```python predict.py -i ECFP.csv -m ECFP```

>SYBYL : ```python predict.py -i SYBYL.csv -m SYBYL```

>ECFP_SYBYL : ```python predict.py -i ECFP_SYBYL.csv -m ECFP_SYBYL```

>ECFP_num : ```python predict.py -i ECFP_num.csv -m ECFP_num``` 


