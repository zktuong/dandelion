Dandelion for TCR gamma/delta reannotation
==========================================

.. figure:: img/dandelion_logo_illustration.png
   :alt: dandelion_logo

   dandelion_logo

In Cell Ranger 3.1.0, the VDJ algorithm was changed to favour TCR
alpha/beta annotation. Since then, calling gamma/delta chains has become
challenging, and 10X support recommends using Cell Ranger 3.0.2 when
working with gamma/delta-rich libraries.

However, the contigs themselves are still accurately reconstructed, just
not annotated correctly. It may be desirable to use a newer Cell Ranger
version for access to some previously unavailable run options, like
specifying custom enrichment primers. In those cases, the contigs can be
reannotated via ``dandelion`` to yield functional output.

Just follow `standard
protocol <https://sc-dandelion.readthedocs.io/en/latest/notebooks/singularity_preprocessing.html>`__
for preparing and running the preprocessing. The parameterisation
recommendation is applicable here as well:

.. code:: bash

   singularity run -B $PWD /path/to/sc-dandelion_latest.sif dandelion-preprocess \
        --chain TR \
        --file_prefix all \
        --filter_to_high_confidence

