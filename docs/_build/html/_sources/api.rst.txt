Preprocessing: `pp`
===================
.. module:: dandelion.preprocessing

.. autosummary::
   :toctree: modules

   format_fasta
   format_fastas
   assign_isotype
   assign_isotypes
   reannotate_genes
   reassign_alleles
   create_germlines
   filter_bcr
   quantify_mutations
   calculate_threshold

Preprocessing (external): `pp.external`
=======================================
.. module:: dandelion.preprocessing.external

.. autosummary::
   :toctree: modules

   assigngenes_igblast
   makedb_igblast
   parsedb_heavy
   parsedb_light
   parsedb_light
   creategermlines
   tigger_genotype
   recipe_scanpy_qc

Tools: `tl`
===========
.. module:: dandelion.tools

.. autosummary::
   :toctree: modules

   find_clones
   define_clones
   clone_size
   clone_overlap
   transfer
   generate_network
   clone_degree
   clone_centrality
   extract_edge_weights
   clone_diversity
   clone_rarefaction

Plotting: `pl`
==============
.. module:: dandelion.plotting

.. autosummary::
   :toctree: modules

   clone_rarefaction
   clone_network
   barplot
   stackedbarplot
   spectratype
   clone_overlap

Utilities: `utl`
================
.. module:: dandelion.utilities

.. autosummary::
   :toctree: modules

   makeblastdb
   load_data
   Dandelion
   update_metadata
   read_h5
   read_pkl

