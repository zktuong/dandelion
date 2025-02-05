Preprocessing: `pp`
===================
.. module:: dandelion.preprocessing

.. autosummary::
   :toctree: modules

   assign_isotype
   assign_isotypes
   calculate_threshold
   check_contigs
   create_germlines
   filter_contigs
   format_fasta
   format_fastas
   quantify_mutations
   reannotate_genes
   reassign_alleles

Preprocessing (external): `pp.external`
=======================================
.. module:: dandelion.preprocessing.external

.. autosummary::
   :toctree: modules

   assigngenes_igblast
   creategermlines
   makedb_igblast
   parsedb_heavy
   parsedb_light
   recipe_scanpy_qc
   tigger_genotype

Tools: `tl`
===========
.. module:: dandelion.tools

.. autosummary::
   :toctree: modules


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
   :toctree: modules

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
   :toctree: modules

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
   :toctree: modules

   Dandelion

.. module:: dandelion.Dandelion

.. autosummary::
   :toctree: modules

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
=========
.. module:: dandelion.logging

.. autosummary::
   :toctree: modules

   print_header
   print_versions
