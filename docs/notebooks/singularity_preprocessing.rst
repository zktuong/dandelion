Dandelion preprocessing with Singularity
========================================

.. figure:: img/dandelion_logo_illustration.png
   :alt: dandelion_logo

   dandelion_logo

Arguably the greatest strength of the dandelion package is a streamlined
preprocessing setup making use of a variety of specialised single cell
VDJ algorithms:

-  V(D)J gene reannotation with ``igblastn`` and parsed to AIRR format
   with
   `changeo’s <https://changeo.readthedocs.io/en/stable/examples/10x.html>`__
   ``MakeDB.py``
   `[Gupta2015] <https://academic.oup.com/bioinformatics/article/31/20/3356/195677>`__,
   with the pipeline strengthened by running ``blastn`` in parallel and
   using the best alignments
-  Reassigning heavy chain IG V gene alleles with
   `TIgGER <https://tigger.readthedocs.io/en/stable/>`__
   `[Gadala-Maria15] <https://www.pnas.org/content/112/8/E862>`__
-  Reassigning IG constant region calls by blasting against a curated
   set of highly specific C gene sequences

However, running this workflow requires a high number of dependencies
and databases, which can be troublesome to set up. As such, we’ve put
together a Singularity container that comes pre-configured with all of
the required software and resources, allowing you to run the
pre-processing pipeline with a single call and easy installation.

Setup and running
-----------------

Once you have `Singularity
installed <https://sylabs.io/guides/3.0/user-guide/installation.html>`__,
you can download the Dandelion container. This command will create
``sc-dandelion_latest.sif``, note its location.

::

   singularity pull library://kt16/default/sc-dandelion:latest

In order to prepare your BCR data for ingestion, create a folder for
each sample you’d like to analyse, name it with your sample ID, and
store the Cell Ranger ``filtered_contig_annotations.csv`` and
``filtered_contig.fasta`` output files inside.

::

   5841STDY7998693
   ├── filtered_contig_annotations.csv
   └── filtered_contig.fasta

You can then navigate to the directory holding all your sample folders
and run Dandelion pre-processing like so:

.. code:: bash

   singularity run -B $PWD /path/to/sc-dandelion_latest.sif dandelion-preprocess [optional arguments follow here]

If you’re running TR data rather than IG data, specify ``--chain TR``.

If you wish to process files that have a different prefix than
``filtered``, e.g. ``all_contig_annotations.csv`` and
``all_contig.fasta``, provide the desired file prefix with
``--file_prefix``. In that case, be sure that your input folder contains
those files rather than the filtered ones. You can also provide the
``--filter_to_high_confidence`` flag to only keep the contigs that Cell
Ranger has called as high confidence.

Recommended parameterisation
----------------------------

If in possession of gene expression data that the BCR data will be
integrated with, the following parameterisation is likely to yield the
best results:

.. code:: bash

   singularity run -B $PWD /path/to/sc-dandelion_latest.sif dandelion-preprocess \
        --file_prefix all \
        --filter_to_high_confidence

Part of the Cell Ranger VDJ filtering criteria is whether the algorithm
thinks the contig is in a cell or not, for which you will have superior
information based on the gene expression data. The other half of the
Cell Ranger VDJ filtering process, requiring the contig to be high
confidence, is retained by providing the ``--filter_to_high_confidence``
flag.

Optional arguments
------------------

By default, this workflow will analyse all provided IG samples jointly
with TIgGER to maximise inference power, and in the event of multiple
input folders will prepend the sample IDs to the cell barcodes to avoid
erroneously merging barcodes overlapping between samples at this stage.
TIgGER should be ran on a per-individual level. If running the workflow
on multiple individuals’ worth of data at once, or wanting to flag the
cell barcodes in a non-default manner, information can be provided to
the script in the form of a CSV file passed through the ``--meta``
argument:

-  The first row of the CSV needs to be a header identifying the
   information in the columns, and the first column needs to contain
   sample IDs.
-  Barcode flagging can be controlled by an optional
   ``prefix``/``suffix`` column. The pipeline will then add the
   specified prefixes/suffixes to the barcodes of the samples. This may
   be desirable, as corresponding gene expression samples are likely to
   have different IDs, and providing the matched ID will pre-format the
   BCR output to match the GEX nomenclature.
-  Individual information for TIgGER can be specified in an optional
   ``individual`` column. If specified, TIgGER will be ran for each
   unique value present in the column, pooling the corresponding
   samples.

It’s possible to just pass a prefix/suffix or individual information. An
excerpt of a sample CSV file that could be used on input:

::

   sample,suffix,individual
   5841STDY7998693,5841STDY7991475,A37
   5841STDY7998694,5841STDY7991476,A37
   5841STDY7998695,5841STDY7991477,A37
   WSSS_A_LNG9030827,WSSS_A_LNG8986832,A51
   WSSS8090101,WSSS8015042,A40
   WSSS8090102,WSSS8015043,A40
   [...]

The delimiter between the barcode and the prefix/suffix can be
controlled with the ``--sep`` argument. By default, the workflow will
strip out the trailing ``"-1"`` from the Cellranger ouput barcode names;
pass ``--keep_trailing_hyphen_number`` if you don’t want to do that.
Pass ``--clean_output`` if you want to remove intermediate files and
just keep the primary output. The intermediate files may be useful for
more detailed inspection.

Output
------

The main file of interest will be
``dandelion/filtered_contig_dandelion.tsv``, stored in a new subfolder
each sample folder. This is an AIRR formatted export of the corrected
contigs, which can be used for downstream analysis by both dandelion
itself, and other packages like
`scirpy <https://icbi-lab.github.io/scirpy/generated/scirpy.io.read_airr.html>`__
`[Sturm2020] <https://academic.oup.com/bioinformatics/article/36/18/4817/5866543>`__
and changeo
`[Gupta2015] <https://academic.oup.com/bioinformatics/article/31/20/3356/195677>`__.

The file above features a contig space filtered with immcantation. If
this is not of interest to you and you wish to see the full contig space
as provided on input, refer to
``dandelion/tmp/filtered_contig_iblastn_db-all.tsv``.

The plots showing the impact of TIgGER are in
``<tigger>/<tigger>_reassign_alleles.pdf``, for each TIgGER folder (one
per unique individual if using ``--meta``, ``tigger`` otherwise). The
impact of C gene reannotation is shown in
``dandelion/data/assign_isotype.pdf`` for each sample.

If you’re interested in more detail about the pre-processing this
offers, or wish to use the workflow in a more advanced manner (e.g. by
using your own databases), proceed to the pre-processing section of the
tutorial.

