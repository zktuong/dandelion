# dandelion

## Intro
Hi there! I've put together a python package for analyzing single cell BCR/V(D)J data from 10x Genomics 5' solution! It streamlines the pre-processing of immcantation tools for single-cell BCR analysis and includes a couple of functions for visualization.

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

## Installation instructions
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

# I use pip to install the following. use --no-cache-dir if necessary
pip install scanpy
pip install scrublet
pip install changeo

#  in case you need to link up kernel with jupyter notebook
python -m ipykernel install --user --name dandelion --display-name "Python (dandelion)"

# download igblast and blast+ from
https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

# unpack where relevant and export the path to the softwares, e.g. ~/Documents/
echo 'export PATH=~/Documents/software/bin:$PATH' >> ~/.bash_profile
source ~/.bash_profile

# set up environmental variables
export GERMLINE=/path/to/dandelion/database/germlines/
export IGDATA=/path/to/dandelion/database/igblast/
export BLASTDB=/path/to/dandelion/database/blast/
export PATH=/path/to/dandelion/bin:$PATH

# and now youd should be good to go!
``` 

## Acknowledgements
I would like to acknowledge the contributions from Dr. Ondrej Suschanek, Dr. Benjamin Stewart, Dr. Rachel Bashford-Rogers and Prof. Menna Clatworthy who helped with the initial conception of the project and for all discussions. 

I would also like to acknowledge Dr. Jongeun Park, Dr. Cecilia-Dominguez Conde, Dr. Hamish King, Dr. Krysztof Polanksi, Dr. Peng He who have had useful discussions. I would also like to thank my wife who helped name the package, because she thought the plots looked like a dandelion =D.

If there is any ideas, comments, suggestions, thing you would like to know more, please feel free to email me at kt16@sanger.ac.uk or post in the issue tracker and we can work something out.