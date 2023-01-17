[![](https://readthedocs.org/projects/sc-dandelion/badge/?version=latest)](https://sc-dandelion.readthedocs.io/en/latest/?badge=latest)
[![](https://img.shields.io/pypi/v/sc-dandelion?logo=PyPI)](https://pypi.org/project/sc-dandelion/)
[![](https://byob.yarr.is/zktuong/dandelion/master-version)](https://github.com/zktuong/dandelion/tree/master)
[![master](https://github.com/zktuong/dandelion/workflows/tests/badge.svg?branch=master)]((https://github.com/zktuong/dandelion/actions?query=workflow%3Atests))
[![codecov](https://codecov.io/gh/zktuong/dandelion/branch/master/graph/badge.svg?token=661BMU1FBO)](https://codecov.io/gh/zktuong/dandelion)
[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/zktuong/dandelion/blob/master/container/dandelion_singularity.ipynb)

![](docs/notebooks/img/dandelion_logo_illustration.png)

Hi there! I have put together a python package for analyzing single cell BCR/TCR data from 10x Genomics 5' solution! It streamlines the pre-processing, leveraging some tools from immcantation suite, and integrates with scanpy/anndata for single-cell BCR/TCR analysis. It also includes a couple of functions for visualization. 

## Citation
`dandelion` is now included in the the following manuscript published in [***Nature Medicine***](https://www.nature.com/articles/s41591-021-01329-2):

*Emily Stephenson, Gary Reynolds, Rachel A Botting, Fernando J Calero-Nieto, Michael Morgan, Zewen Kelvin Tuong, Karsten Bach, Waradon Sungnak, Kaylee B Worlock, Masahiro Yoshida, Natsuhiko Kumasaka, Katarzyna Kania, Justin Engelbert, Bayanne Olabi, Jarmila Stremenova Spegarova, Nicola K Wilson, Nicole Mende, Laura Jardine, Louis CS Gardner, Issac Goh, Dave Horsfall, Jim McGrath, Simone Webb, Michael W Mather, Rik GH Lindeboom, Emma Dann, Ni Huang, Krzysztof Polanski, Elena Prigmore, Florian Gothe, Jonathan Scott, Rebecca P Payne, Kenneth F Baker, Aidan T Hanrath, Ina CD Schim van der Loeff, Andrew S Barr, Amada Sanchez-Gonzalez, Laura Bergamaschi, Federica Mescia, Josephine L Barnes, Eliz Kilich, Angus de Wilton, Anita Saigal, Aarash Saleh, Sam M Janes, Claire M Smith, Nusayhah Gopee, Caroline Wilson, Paul Coupland, Jonathan M Coxhead, Vladimir Y Kiselev, Stijn van Dongen, Jaume Bacardit, Hamish W King, Anthony J Rostron, A John Simpson, Sophie Hambleton, Elisa Laurenti, Paul A Lyons, Kerstin B Meyer, Marko Z Nikolic, Christopher JA Duncan, Ken Smith, Sarah A Teichmann, Menna R Clatworthy, John C Marioni, Berthold Gottgens, Muzlifah Haniffa.* ***Single-cell multi-omics analysis of the immune response in COVID-19***. *Nature Medicine 2021.04.20; doi: https://dx.doi.org/10.1038/s41591-021-01329-2*

Please also cite the following pre-print if you use version 0.3.0 onwards:

*Chenqu Suo, Krzysztof Polanski, Emma Dann, Rik G.H. Lindeboom, Roser Vilarrasa Blasi, Roser Vento-Tormo, Muzlifah Haniffa, Kerstin B Meyer, Zewen Kelvin Tuong, Menna R. Clatworthy, Sarah Teichmann.* ***Single cell antigen receptor analysis reveals lymphocyte developmental origins***. *bioRxiv 2022.11.18.517068; doi: https://doi.org/10.1101/2022.11.18.517068*

## Overview

![](docs/notebooks/img/dandelion_overview.png)

Illustration of the `Dandelion` class slots

![](docs/notebooks/img/dandelion_class2.png)

Please refer to the [documentation](https://sc-dandelion.readthedocs.io/).

The raw files for the examples can be downloaded from 10X's Single Cell Immune Profiling datasets [website](https://support.10xgenomics.com/single-cell-vdj/datasets).

## Installation

### Singularity container

`dandelion` now comes ready in the form of a singularity container:
```bash
singularity pull library://kt16/default/sc-dandelion:latest
singularity shell --writable-tmpfs -B $PWD sc-dandelion_latest.sif
```
This will load up a container that has all the required dependencies installed.

This can be used for the preprocessing steps by navigating to the data folder and use:
```bash
singularity run -B $PWD sc-dandelion_latest.sif dandelion-preprocess
```

If you have multiple samples to process, it is reccomended to specify the `--meta` option with a `.csv` file detailing the samples:
```bash
singularity run -B $PWD sc-dandelion_latest.sif dandelion-preprocess --meta meta.csv
```

### Python package

Start off by creating a conda environment containing scanpy, following [official scanpy instructions](https://scanpy.readthedocs.io/en/stable/installation.html). Once done, run the following:

```bash
conda install -c conda-forge graph-tool
pip install sc-dandelion
```

Between this and the pipelines within the singularity container, you should be covered for most of your needs.

## Basic requirements
Python packages
```python
# conda
python>=3.7 (conda-forge)
numpy>=1.18.4 (conda-forge)
pandas>=1.0.3 (conda-forge)
distance>=0.1.3 (conda-forge)
jupyter (conda-forge) # if running via a notebook
scikit-learn>=0.23.0 (conda-forge)
numba>=0.48.0 (conda-forge)
pytables>=3.6.1 (conda-forge)
seaborn>=0.10.1 (conda-forge)
leidenalg>=0.8.0 (conda-forge)
plotnine>=0.6.0 (conda-forge)
graph-tool>=2.3.5 (conda-forge) # optional

# Other executables (through conda)
blast>=2.10.1 (bioconda)
igblast>=1.15.0 (bioconda)

# pip
anndata>=0.7.1
scanpy>=1.4.6
scrublet>=0.2.1
changeo>=1.0.0
presto>=0.6.0
polyleven>=0.5
networkx>=2.4
rpy2>=3.4.2
```

## Acknowledgements
I would like to acknowledge the contributions from Dr. Chenqu Suo, Dr. Krysztof Polanksi, Dr. Sarah Teichmann and Prof. Menna Clatworthy, who helped with the initial conception of the project and for all discussions.

I would also like to acknowledge Dr. Ondrej Suschanek, Dr. Benjamin Stewart, Dr. Rachel Bashford-Rogers, Dr. Jongeun Park,  Dr. Cecilia-Dominguez Conde, Dr. Kirsten Stewart, Dr. Hamish King and Dr. Peng He with whom I have had very useful discussions. I would also like to thank my wife who helped name the package, because she thought the plots looked like a dandelion =D.

## Support
 
Support is provided on a voluntary basis, as time permits.

If there are any ideas, comments, suggestions, thing you would like to know more etc., please feel free to email me at z.tuong@uq.edu.au or post in the issue tracker and I will get back to you.
