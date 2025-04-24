Preprocessing: `pp`
===================
.. module:: dandelion.preprocessing

.. autosummary::
   :toctree: .

   assign_isotype
   assign_isotypes
   check_contigs
   create_germlines
   filter_contigs
   format_fasta
   format_fastas
   reannotate_genes
   reassign_alleles

Tools: `tl`
===========
.. module:: dandelion.tools

.. autosummary::
   :toctree: .


   transfer
   find_clones
   define_clones
   generate_network
   clone_centrality
   clone_degree
   clone_diversity
   clone_overlap
   clone_rarefaction
   clone_size
   extract_edge_weights
   productive_ratio
   vj_usage_pca
   setup_vdj_pseudobulk
   vdj_pseudobulk
   pseudobulk_gex
   pseudotime_transfer
   project_pseudotime_to_cell
   bin_expression
   chatterjee_corr


Plotting: `pl`
==============
.. module:: dandelion.plotting

.. autosummary::
   :toctree: .

   barplot
   clone_network
   clone_overlap
   clone_rarefaction
   productive_ratio
   spectratype
   stackedbarplot

Utilities: `utl`
================
.. module:: dandelion.utilities

.. autosummary::
   :toctree: .

   concat
   load_data
   read_h5ddl
   read_pkl
   read_airr
   read_10x_airr
   read_10x_vdj
   read_parse_airr
   read_bd_airr
   to_scirpy
   from_scirpy
   makeblastdb


Dandelion
=========
.. module:: dandelion

.. autosummary::
   :toctree: .
   :recursive:

   Dandelion

dandelion.Dandelion
-------------------

.. currentmodule:: dandelion.Dandelion

.. autoclass:: Dandelion
   :members:
   :undoc-members:
   :show-inheritance:

.. autosummary::
   :toctree: .
   :recursive:

   add_cell_prefix
   add_cell_suffix
   add_sequence_prefix
   add_sequence_suffix
   copy
   reset_ids
   simplify
   store_germline_reference
   update_metadata
   update_plus
   write
   write_10x
   write_airr
   write_h5ddl
   write_h5ddl_legacy
   write_pkl

Logging
=======
.. module:: dandelion.logging

.. autosummary::
   :toctree: .

   print_header
   print_versions

External
========

scanpy
------
.. module:: dandelion.external.scanpy
.. autosummary::
   :toctree: .

   recipe_scanpy_qc

Immmcantation
-------------

Wrappers for tools in Immcantation pipeline.

changeo
~~~~~~~
.. module:: dandelion.external.immcantation.changeo

.. autosummary::
   :toctree: .

   assigngenes_igblast
   creategermlines
   makedb_igblast
   parsedb_heavy
   parsedb_light

tigger
~~~~~~
.. module:: dandelion.external.immcantation.tigger

.. autosummary::
   :toctree: .

   tigger_genotype


shazam
~~~~~~
.. module:: dandelion.external.immcantation.shazam

.. autosummary::
   :toctree: .

   calculate_threshold
   quantify_mutations

scoper
~~~~~~
.. module:: dandelion.external.immcantation.scoper

.. autosummary::
   :toctree: .

   identical_clones
   hierarchical_clones
   spectral_clones
