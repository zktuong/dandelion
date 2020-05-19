# dandelion

![dandelion_logo](notebook/img/dandelion_logo.png)

## Intro
Hi there! I've put together a python package for analyzing single cell BCR/V(D)J data from 10x Genomics 5' solution! It streamlines the pre-processing of immcantation tools for single-cell BCR analysis and includes a couple of functions for visualization.

## Notebooks
Please see [notebooks](notebook/) for examples.

## Installation instructions

Try this first
```bash
# create a conda environment with specific modules
conda create --name dandelion python=3.7
conda activate dandelion

pip install git+https://github.com/zktuong/dandelion.git
```

If it fails, try installing the dependencies first.

I would reccomend instaling this in order
```bash
# create a conda environment with specific modules
conda create --name dandelion python=3.7
conda activate dandelion

# the following two are what's required by scanpy
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leidenalg 
# these are required by dandelion
conda install distance python-Levenshtein joblib jupyter

# Use pip to install the following with --no-cache-dir if necessary
pip install scanpy
pip install scrublet
pip install changeo

# and then lastly install this
pip install git+https://github.com/zktuong/dandelion.git
````
In case you need to link up kernel with jupyter notebook
```bash
python -m ipykernel install --user --name dandelion --display-name "Python (dandelion)"
```


dandelion also requires some R packages intalled.
```R
install.packages(c("alakazam", "tigger", "airr", "shazam", "ggplot2"))
```


## Required database
Last but not least, you will need download the database folder in this repository manually and place them somewhere accessible and export them as environmental variables in your `~/.bash_profile` so that dandelion and the blast programs can access them properly.

So for example, clone this repository into `~/Documents`
```bash
# set up environmental variables
export GERMLINE=~/Documents/dandelion/database/germlines/
export IGDATA=~/Documents/dandelion/database/igblast/
export BLASTDB=~/Documents/dandelion/database/blast/

# and now you should be good to go! You can remove the rest of the cloned folder (other then the database folder) in ~/Documents as the package is installed.
``` 


## External softwares
I have already included in this repository the binaries for the various blast executables and pip should pull them and install prorperl. However, you can also download them yourself and store the softwares somewhere more accessible if you like. Just make sure to set the path to them appropriately.
```bash
# download igblast and blast+ from
https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

# unpack where relevant and export the path to the softwares, e.g. ~/Documents/
echo 'export PATH=~/Documents/software/bin:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```


## Requirements
Python packages
```python
# conda
python==3.7.6 (conda-forge)
numpy==1.18.4 (conda-forge)
pandas==1.0.3 (conda-forge)
python-Levenshtein==0.12.0 (conda-forge)
distance==0.1.3 (conda-forge)
joblib==0.14.1 (conda-forge)
jupyter==1.0.0 (conda-forge)
scikit-learn==0.23.0 (conda-forge)
numba==0.48.0 (conda-forge)
pytables==3.6.1 (conda-forge)
seaborn==0.10.1 (conda-forge)
python-igraph==0.8.2 (conda-forge)
leidenalg==0.8.0 (conda-forge)

# pip
changeo==1.0.0
anndata==0.7.1
scanpy==1.4.6
scrublet==0.2.1
```
R packages
```R
alakazam_1.0.1
tigger_1.0.0
airr_1.2.0
shazam_1.0.0
ggplot2
```

## Acknowledgements
I would like to acknowledge the contributions from Dr. Ondrej Suschanek, Dr. Benjamin Stewart, Dr. Rachel Bashford-Rogers and Prof. Menna Clatworthy who helped with the initial conception of the project and for all discussions. 

I would also like to acknowledge Dr. Jongeun Park, Dr. Cecilia-Dominguez Conde, Dr. Hamish King, Dr. Krysztof Polanksi, Dr. Peng He who have had useful discussions. I would also like to thank my wife who helped name the package, because she thought the plots looked like a dandelion =D.

If there is any ideas, comments, suggestions, thing you would like to know more, please feel free to email me at kt16@sanger.ac.uk or post in the issue tracker and we can work something out.