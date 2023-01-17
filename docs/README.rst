|Docs| |PyPI| |Master| |MasterTest| |CodeCov| |Colab|

|logo|

Hi there! I have put together a python package for analyzing single cell
BCR/TCR data from 10x Genomics 5' solution! It streamlines the
pre-processing, leveraging some tools from immcantation suite, and
integrates with scanpy/anndata for single-cell BCR/TCR analysis. It also
includes a couple of functions for visualization. Try it out on |Colab| !

Overview
--------

|overview|

Illustration of the ``Dandelion`` class slots

|class|

Please refer to the
`documentation <https://sc-dandelion.readthedocs.io/>`__ or the
notebooks
`here <https://nbviewer.jupyter.org/github/zktuong/dandelion/tree/latest/docs/notebooks/>`__:

The raw files for the examples can be downloaded from 10X's Single Cell
Immune Profiling datasets
`website <https://support.10xgenomics.com/single-cell-vdj/datasets>`__.


.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/pxz31b1iIFY" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>


Installation
------------

Singularity container
~~~~~~~~~~~~~~~~~~~~~

``dandelion`` now comes ready in the form of a singularity container 
which has all the required dependencies installed:

.. code:: bash

    singularity pull library://kt16/default/sc-dandelion:latest
    singularity shell --writable-tmpfs -B $PWD sc-dandelion_latest.sif


This will load up a conda-environment that has all the required 
dependencies installed.

This can be used for the preprocessing steps by navigating to the data 
folder and use:

.. code:: bash

    singularity run -B $PWD sc-dandelion_latest.sif dandelion-preprocess

Python package
~~~~~~~~~~~~~~

Start off by creating a conda environment containing scanpy, following
`official scanpy instructions <https://scanpy.readthedocs.io/en/stable/installation.html>`__.
Once done, run the following:

.. code:: bash

    conda install -c conda-forge graph-tool
    pip install sc-dandelion


Between this and the pipelines within the singularity container, you 
should be covered for most of your needs.

Basic requirements
------------------

Python packages

.. code:: python

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


Acknowledgements
----------------

I would like to acknowledge the contributions from Dr. Chenqu Suo, Dr. 
Krysztof Polanksi, Dr. Sarah Teichmann and Prof. Menna Clatworthy, who 
helped with the initial conception of the project and for all discussions.

I would also like to acknowledge Dr. Ondrej Suschanek,
Dr. Benjamin Stewart, Dr. Rachel Bashford-Rogers, Dr. Jongeun Park, 
Dr. Cecilia-Dominguez Conde, Dr. Kirsten Stewart, Dr. Hamish King and 
Dr. Peng He with whom I have had very useful discussions. I would also 
like to thank my wife who helped name the package, because she thought 
the plots looked like a dandelion =D.

Support
-------

Support is provided on a voluntary basis, as time permits.

If there are any ideas, comments, suggestions, thing you would like to
know more etc., please feel free to email me at z.tuong@uq.edu.au or
post in the issue tracker and I will get back to you.

Citation
--------

``dandelion`` is now included in the the following manuscript:

.. [Stephenson2021] Stephenson *et al.* (2021),
   *Single-cell multi-omics analysis of the immune response in COVID-19*,
   `Nature Medicine <https://www.nature.com/articles/s41591-021-01329-2>`__.

*Emily Stephenson, Gary Reynolds, Rachel A Botting, Fernando J
Calero-Nieto, Michael Morgan, Zewen Kelvin Tuong, Karsten Bach, Waradon
Sungnak, Kaylee B Worlock, Masahiro Yoshida, Natsuhiko Kumasaka,
Katarzyna Kania, Justin Engelbert, Bayanne Olabi, Jarmila Stremenova
Spegarova, Nicola K Wilson, Nicole Mende, Laura Jardine, Louis CS
Gardner, Issac Goh, Dave Horsfall, Jim McGrath, Simone Webb, Michael W
Mather, Rik GH Lindeboom, Emma Dann, Ni Huang, Krzysztof Polanski, Elena
Prigmore, Florian Gothe, Jonathan Scott, Rebecca P Payne, Kenneth F
Baker, Aidan T Hanrath, Ina CD Schim van der Loeff, Andrew S Barr, Amada
Sanchez-Gonzalez, Laura Bergamaschi, Federica Mescia, Josephine L
Barnes, Eliz Kilich, Angus de Wilton, Anita Saigal, Aarash Saleh, Sam M
Janes, Claire M Smith, Nusayhah Gopee, Caroline Wilson, Paul Coupland,
Jonathan M Coxhead, Vladimir Y Kiselev, Stijn van Dongen, Jaume
Bacardit, Hamish W King, Anthony J Rostron, A John Simpson, Sophie
Hambleton, Elisa Laurenti, Paul A Lyons, Kerstin B Meyer, Marko Z
Nikolic, Christopher JA Duncan, Ken Smith, Sarah A Teichmann, Menna R
Clatworthy, John C Marioni, Berthold Gottgens, Muzlifah Haniffa.*
**Single-cell multi-omics analysis of the immune response in
COVID-19**. *Nature Medicine 2021.04.20; doi:
https://dx.doi.org/10.1038/s41591-021-01329-2*

Please also cite the following pre-print if you use version 0.3.0 onwards:

.. [Suo2022] Suo *et al.* (2022),
   *Single cell antigen receptor analysis reveals lymphocyte developmental origins*,
   `bioRxiv <https://doi.org/10.1101/2022.11.18.517068>`__.


*Chenqu Suo, Krzysztof Polanski, Emma Dann, Rik G.H. Lindeboom, Roser Vilarrasa Blasi, Roser Vento-Tormo, Muzlifah Haniffa, Kerstin B Meyer, Zewen Kelvin Tuong, Menna R. Clatworthy, Sarah Teichmann.* **Single cell antigen receptor analysis reveals lymphocyte developmental origins**. *bioRxiv 2022.11.18.517068; doi: https://doi.org/10.1101/2022.11.18.517068*


If you use the pre-processing tools/functions, please cite the relevant manuscripts from the immcantation suite, including:

.. [changeo]
*Gupta NT, Vander Heiden JA, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH.* **Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data.** *Bioinformatics 31, 3356-8 (2015). doi: https://doi.org/10.1093/bioinformatics/btv359*

.. [tigger]
*Gadala-Maria D, Yaari G, Uduman M, Kleinstein SH.* **Automated analysis of high-throughput B cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles.** *Proceedings of the National Academy of Sciency of the United States of America, E862-70.*

References
----------

.. [Bashford-Rogers2013] Bashford-Rogers *et al.* (2013),
   *Network properties derived from deep sequencing of human B-cell receptor repertoires delineate B-cell populations*,
   `Genome Research <https://genome.cshlp.org/content/23/11/1874>`__.

.. [Bashford-Rogers2019] Bashford-Rogers *et al.* (2019),
   *Analysis of the B cell receptor repertoire in six immune-mediated diseases*,
   `Nature <https://www.nature.com/articles/s41586-019-1595-3>`__.

.. [Dann2022] Dann *et al.* (2022),
   *Differential abundance testing on single-cell data using k-nearest neighbor graphs*,
   `Nature Biotechnology <https://doi.org/10.1038/s41587-021-01033-z>`__.
   `GitHub <https://github.com/emdann/milopy>`__.

.. [Gadala-Maria2015] Gadala-Maria *et al.* (2015),
   *Automated analysis of high-throughput B cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles*,
   `Proceedings of the National Academy of Sciency of the United States of America <https://www.pnas.org/content/112/8/E862>`__.

.. [Gupta2015] Gupta *et al.* (2015),
   *Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data*,
   `Bioinformatics <https://academic.oup.com/bioinformatics/article/31/20/3356/195677>`__.

.. [Setty2019] Setty *et al.* (2019)
   *Characterization of cell fate probabilities in single-cell data with Palantir*,
   `Nature Biotechnology <https://doi.org/10.1038/s41587-019-0068-4>`__.
   `GitHub <https://github.com/dpeerlab/Palantir>`__.

.. [Sleckman1998] Sleckman *et al.* (1998)
   *Assembly of productive T cell receptor delta variable region genes exhibits allelic inclusion*,
   `Journal of Experimental Medicine <https://rupress.org/jem/article-lookup/doi/10.1084/jem.188.8.1465>`__.

.. [Stephenson2021] Stephenson *et al.* (2021),
   *Single-cell multi-omics analysis of the immune response in COVID-19*,
   `Nature Medicine <https://www.nature.com/articles/s41591-021-01329-2>`__.

.. [Sturm2020] Sturm *et al.* (2020),
   *Scirpy: a Scanpy extension for analyzing single-cell T-cell receptor-sequencing data*,
   `Bioinformatics <https://academic.oup.com/bioinformatics/article/36/18/4817/5866543>`__.
   `GitHub <https://github.com/icbi-lab/scirpy>`__.

.. [Suo2022] Suo *et al.* (2022),
   *Single cell antigen receptor analysis reveals lymphocyte developmental origins*,
   `bioRxiv <https://doi.org/10.1101/2022.11.18.517068>`__.

.. [Wolf2018] Wolf *et al.* (2018),
   *Scanpy: large-scale single-cell gene expression data analysis*,
   `Genome Biology <https://doi.org/10.1186/s13059-017-1382-0>`__.
   `GitHub <https://github.com/theislab/scanpy>`__.


.. |Docs| image:: https://readthedocs.org/projects/sc-dandelion/badge/?version=latest
   :target: https://sc-dandelion.readthedocs.io/en/latest/?badge=latest
.. |PyPI| image:: https://img.shields.io/pypi/v/sc-dandelion?logo=PyPI
   :target: https://pypi.org/project/sc-dandelion/
.. |Master| image:: https://byob.yarr.is/zktuong/dandelion/master-version
   :target: https://github.com/zktuong/dandelion/tree/master
.. |MasterTest| image:: https://github.com/zktuong/dandelion/workflows/tests/badge.svg?branch=master
   :target: https://github.com/zktuong/dandelion/actions/workflows/tests.yml
.. |CodeCov| image:: https://codecov.io/gh/zktuong/dandelion/branch/master/graph/badge.svg?token=661BMU1FBO
   :target: https://codecov.io/gh/zktuong/dandelion
.. |Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/zktuong/dandelion/blob/master/container/dandelion_singularity.ipynb
.. |logo| image:: notebooks/img/dandelion_logo_illustration.png
.. |overview| image:: notebooks/img/dandelion_overview.png
.. |class| image:: notebooks/img/dandelion_class2.png
