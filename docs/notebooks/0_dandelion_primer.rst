*Dandelion* class
=================

.. figure:: img/dandelion_logo_illustration.png
   :alt: dandelion_logo

   dandelion_logo

Much of the functions and utility of the ``dandelion`` package revolves
around the ``Dandelion`` class object. The class will act as an
intermediary object for storage and flexible interaction with other
tools. This section will run through a quick primer to the ``Dandelion``
class.

**Import modules**

.. code:: ipython3

    import os
    os.chdir(os.path.expanduser('/Users/kt16/Downloads/dandelion_tutorial/'))
    import dandelion as ddl
    ddl.logging.print_versions()


.. parsed-literal::

    dandelion==0.2.4.dev58 pandas==1.4.2 numpy==1.21.6 matplotlib==3.5.2 networkx==2.8.4 scipy==1.8.1


.. code:: ipython3

    vdj = ddl.read_h5ddl('dandelion_results.h5ddl')
    vdj


.. parsed-literal::

    /var/folders/nb/wrd6px6171j52lqpmkljt6vw000l2l/T/ipykernel_44796/3199686786.py:1: DeprecationWarning: read_h5 is a deprecated in 0.2.2 and will be removed in 0.4.0. read_h5ddl will be the recommended way to read.




.. parsed-literal::

    Dandelion class object with n_obs = 2773 and n_contigs = 5609
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id', 'fwr1_VJ', 'fwr1_VDJ', 'mu_count_VDJ', 'mu_count_VJ', 'mu_count', 'junction_length_VDJ', 'junction_length_VJ', 'junction_aa_length_VDJ', 'junction_aa_length_VJ', 'np1_length_VDJ', 'np1_length_VJ', 'np2_length_VDJ'
        layout: layout for 2773 vertices, layout for 1067 vertices
        graph: networkx graph of 2773 vertices, networkx graph of 1067 vertices 



Basically, the object can be summarized in the following illustration:
|dandelion_class <|

.. |dandelion_class <| image:: img/dandelion_class2.png

Essentially, the ``.data`` slot holds the AIRR contig table while the
``.metadata`` holds a collapsed version that is compatible with
combining with ``AnnData``\ ’s ``.obs`` slot. You can retrieve these
slots like a typical class object; for example, if I want the metadata:

.. code:: ipython3

    vdj.metadata




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>clone_id</th>
          <th>clone_id_by_size</th>
          <th>sample_id</th>
          <th>locus_VDJ</th>
          <th>locus_VJ</th>
          <th>productive_VDJ</th>
          <th>productive_VJ</th>
          <th>v_call_genotyped_VDJ</th>
          <th>d_call_VDJ</th>
          <th>j_call_VDJ</th>
          <th>...</th>
          <th>productive_B_VJ</th>
          <th>duplicate_count_B_VDJ</th>
          <th>duplicate_count_B_VJ</th>
          <th>isotype</th>
          <th>isotype_status</th>
          <th>locus_status</th>
          <th>chain_status</th>
          <th>rearrangement_status_VDJ</th>
          <th>rearrangement_status_VJ</th>
          <th>changeo_clone_id</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACCTGTCATATCGG</th>
          <td>None</td>
          <td>None</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>None</td>
          <td>IGK</td>
          <td>None</td>
          <td>T</td>
          <td>None</td>
          <td>None</td>
          <td>None</td>
          <td>...</td>
          <td>T</td>
          <td>NaN</td>
          <td>68.0</td>
          <td>None</td>
          <td>None</td>
          <td>Orphan IGK</td>
          <td>Orphan VJ</td>
          <td>None</td>
          <td>standard</td>
          <td></td>
        </tr>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACCTGTCCGTTGTC</th>
          <td>B_36_3_2_153_2_2</td>
          <td>2191</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV1-69</td>
          <td>IGHD3-22</td>
          <td>IGHJ3</td>
          <td>...</td>
          <td>T</td>
          <td>51.0</td>
          <td>43.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>9_0</td>
        </tr>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACCTGTCGAGAACG</th>
          <td>B_40_1_1_181_1_1</td>
          <td>1172</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>IGH</td>
          <td>IGL</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV1-2</td>
          <td>None</td>
          <td>IGHJ3</td>
          <td>...</td>
          <td>T</td>
          <td>47.0</td>
          <td>90.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGL</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>152_1</td>
        </tr>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACCTGTCTTGAGAC</th>
          <td>B_174_4_3_202_1_1</td>
          <td>1086</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV5-51</td>
          <td>None</td>
          <td>IGHJ3</td>
          <td>...</td>
          <td>T</td>
          <td>80.0</td>
          <td>22.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>325_2</td>
        </tr>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACGGGAGCGACGTA</th>
          <td>B_53_2_1_22_2_7</td>
          <td>1398</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>IGH</td>
          <td>IGL</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV4-4</td>
          <td>IGHD6-13</td>
          <td>IGHJ3</td>
          <td>...</td>
          <td>T</td>
          <td>18.0</td>
          <td>14.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGL</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>293_3</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTCCTCAGCAATATG</th>
          <td>B_82_2_1_41_2_8</td>
          <td>384</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV2-5</td>
          <td>IGHD5/OR15-5b,IGHD5/OR15-5a</td>
          <td>IGHJ4,IGHJ5</td>
          <td>...</td>
          <td>T</td>
          <td>41.0</td>
          <td>71.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>1780_1193</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTCCTCAGCGCTTAT</th>
          <td>B_148_6_5_99_1_3</td>
          <td>400</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV3-30</td>
          <td>IGHD4-17</td>
          <td>IGHJ6</td>
          <td>...</td>
          <td>T</td>
          <td>11.0</td>
          <td>28.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>840_1194</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTCCTCAGGGAAACA</th>
          <td>B_70_1_1_68_4_13</td>
          <td>381</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV4-59</td>
          <td>IGHD6-13</td>
          <td>IGHJ2</td>
          <td>...</td>
          <td>T</td>
          <td>14.0</td>
          <td>159.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>1224_1195</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTGCGCCATACCATG</th>
          <td>B_68_7_1_114_2_6</td>
          <td>380</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGL</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV1-69</td>
          <td>IGHD2-15</td>
          <td>IGHJ6</td>
          <td>...</td>
          <td>T</td>
          <td>32.0</td>
          <td>28.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGL</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>1821_1196</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTGGTTGTAGGCATG</th>
          <td>B_186_5_3_178_3_2</td>
          <td>379</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGL</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV3-23</td>
          <td>None</td>
          <td>IGHJ4</td>
          <td>...</td>
          <td>T</td>
          <td>22.0</td>
          <td>36.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGL</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>1958_1197</td>
        </tr>
      </tbody>
    </table>
    <p>2773 rows × 34 columns</p>
    </div>



slicing
~~~~~~~

You can slice the ``Dandelion`` object via the ``.data`` or
``.metadata`` via their indices, with the behavior similar to how it is
in pandas ``DataFrame`` and ``AnnData``.

Slicing ``.data``
^^^^^^^^^^^^^^^^^

.. code:: ipython3

    vdj[vdj.data['clone_id'] == 'B_36_3_2_153_2_2']




.. parsed-literal::

    Dandelion class object with n_obs = 1 and n_contigs = 2
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id'
        layout: layout for 1 vertices, layout for 0 vertices
        graph: networkx graph of 1 vertices, networkx graph of 0 vertices 



.. code:: ipython3

    vdj[vdj.data_names.isin(['sc5p_v2_hs_PBMC_10k_AAACCTGTCATATCGG_contig_1','sc5p_v2_hs_PBMC_10k_AAACCTGTCCGTTGTC_contig_2','sc5p_v2_hs_PBMC_10k_AAACCTGTCCGTTGTC_contig_1','sc5p_v2_hs_PBMC_10k_AAACCTGTCGAGAACG_contig_1','sc5p_v2_hs_PBMC_10k_AAACCTGTCGAGAACG_contig_2',])]




.. parsed-literal::

    Dandelion class object with n_obs = 3 and n_contigs = 5
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id'
        layout: layout for 3 vertices, layout for 0 vertices
        graph: networkx graph of 3 vertices, networkx graph of 0 vertices 



slicing ``.metadata``
^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    vdj[vdj.metadata['productive_VDJ'].isin(['T','T|T'])]




.. parsed-literal::

    Dandelion class object with n_obs = 2585 and n_contigs = 5372
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id'
        layout: layout for 2585 vertices, layout for 1067 vertices
        graph: networkx graph of 2585 vertices, networkx graph of 1067 vertices 



.. code:: ipython3

    vdj[vdj.metadata_names == 'vdj_v1_hs_pbmc3_TTTCCTCAGCGCTTAT']




.. parsed-literal::

    Dandelion class object with n_obs = 1 and n_contigs = 2
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id'
        layout: layout for 1 vertices, layout for 1 vertices
        graph: networkx graph of 1 vertices, networkx graph of 1 vertices 



copy
~~~~

You can deep copy the ``Dandelion`` object to another variable which
will inherit all slots:

.. code:: ipython3

    vdj2 = vdj.copy()
    vdj2.metadata




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>clone_id</th>
          <th>clone_id_by_size</th>
          <th>sample_id</th>
          <th>locus_VDJ</th>
          <th>locus_VJ</th>
          <th>productive_VDJ</th>
          <th>productive_VJ</th>
          <th>v_call_genotyped_VDJ</th>
          <th>d_call_VDJ</th>
          <th>j_call_VDJ</th>
          <th>...</th>
          <th>productive_B_VJ</th>
          <th>duplicate_count_B_VDJ</th>
          <th>duplicate_count_B_VJ</th>
          <th>isotype</th>
          <th>isotype_status</th>
          <th>locus_status</th>
          <th>chain_status</th>
          <th>rearrangement_status_VDJ</th>
          <th>rearrangement_status_VJ</th>
          <th>changeo_clone_id</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACCTGTCATATCGG</th>
          <td>None</td>
          <td>None</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>None</td>
          <td>IGK</td>
          <td>None</td>
          <td>T</td>
          <td>None</td>
          <td>None</td>
          <td>None</td>
          <td>...</td>
          <td>T</td>
          <td>NaN</td>
          <td>68.0</td>
          <td>None</td>
          <td>None</td>
          <td>Orphan IGK</td>
          <td>Orphan VJ</td>
          <td>None</td>
          <td>standard</td>
          <td></td>
        </tr>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACCTGTCCGTTGTC</th>
          <td>B_36_3_2_153_2_2</td>
          <td>2191</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV1-69</td>
          <td>IGHD3-22</td>
          <td>IGHJ3</td>
          <td>...</td>
          <td>T</td>
          <td>51.0</td>
          <td>43.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>9_0</td>
        </tr>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACCTGTCGAGAACG</th>
          <td>B_40_1_1_181_1_1</td>
          <td>1172</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>IGH</td>
          <td>IGL</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV1-2</td>
          <td>None</td>
          <td>IGHJ3</td>
          <td>...</td>
          <td>T</td>
          <td>47.0</td>
          <td>90.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGL</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>152_1</td>
        </tr>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACCTGTCTTGAGAC</th>
          <td>B_174_4_3_202_1_1</td>
          <td>1086</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV5-51</td>
          <td>None</td>
          <td>IGHJ3</td>
          <td>...</td>
          <td>T</td>
          <td>80.0</td>
          <td>22.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>325_2</td>
        </tr>
        <tr>
          <th>sc5p_v2_hs_PBMC_10k_AAACGGGAGCGACGTA</th>
          <td>B_53_2_1_22_2_7</td>
          <td>1398</td>
          <td>sc5p_v2_hs_PBMC_10k</td>
          <td>IGH</td>
          <td>IGL</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV4-4</td>
          <td>IGHD6-13</td>
          <td>IGHJ3</td>
          <td>...</td>
          <td>T</td>
          <td>18.0</td>
          <td>14.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGL</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>293_3</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTCCTCAGCAATATG</th>
          <td>B_82_2_1_41_2_8</td>
          <td>384</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV2-5</td>
          <td>IGHD5/OR15-5b,IGHD5/OR15-5a</td>
          <td>IGHJ4,IGHJ5</td>
          <td>...</td>
          <td>T</td>
          <td>41.0</td>
          <td>71.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>1780_1193</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTCCTCAGCGCTTAT</th>
          <td>B_148_6_5_99_1_3</td>
          <td>400</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV3-30</td>
          <td>IGHD4-17</td>
          <td>IGHJ6</td>
          <td>...</td>
          <td>T</td>
          <td>11.0</td>
          <td>28.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>840_1194</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTCCTCAGGGAAACA</th>
          <td>B_70_1_1_68_4_13</td>
          <td>381</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGK</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV4-59</td>
          <td>IGHD6-13</td>
          <td>IGHJ2</td>
          <td>...</td>
          <td>T</td>
          <td>14.0</td>
          <td>159.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGK</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>1224_1195</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTGCGCCATACCATG</th>
          <td>B_68_7_1_114_2_6</td>
          <td>380</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGL</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV1-69</td>
          <td>IGHD2-15</td>
          <td>IGHJ6</td>
          <td>...</td>
          <td>T</td>
          <td>32.0</td>
          <td>28.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGL</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>1821_1196</td>
        </tr>
        <tr>
          <th>vdj_v1_hs_pbmc3_TTTGGTTGTAGGCATG</th>
          <td>B_186_5_3_178_3_2</td>
          <td>379</td>
          <td>vdj_v1_hs_pbmc3</td>
          <td>IGH</td>
          <td>IGL</td>
          <td>T</td>
          <td>T</td>
          <td>IGHV3-23</td>
          <td>None</td>
          <td>IGHJ4</td>
          <td>...</td>
          <td>T</td>
          <td>22.0</td>
          <td>36.0</td>
          <td>IgM</td>
          <td>IgM</td>
          <td>IGH + IGL</td>
          <td>Single pair</td>
          <td>standard</td>
          <td>standard</td>
          <td>1958_1197</td>
        </tr>
      </tbody>
    </table>
    <p>2773 rows × 34 columns</p>
    </div>



Retrieving entries with ``update_metadata``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``.metadata`` slot in Dandelion class automatically initializes
whenever the ``.data`` slot is filled. However, it only returns a
standard number of columns that are pre-specified. To retrieve other
columns from the ``.data`` slot, we can update the metadata with
``ddl.update_metadata`` and specify the options ``retrieve`` and
``retrieve_mode``.

The following modes determine how the retrieval is completed:

``split and unique only`` - splits the retrieval into VDJ and VDJ
chains. A ``|`` will separate **unique** element.

``merge and unique only`` - smiliar to above but merged into a single
column.

``split`` - split retrieval into **individual** columns for each contig.

``merge`` - merge retrieval into a **single** column where a ``|`` will
separate **every** element.

For numerical columns, there’s additional options:

``split and sum`` - splits the retrieval into VDJ and VDJ chains and sum
separately.

``split and average`` - smiliar to above but average instead of sum.

``sum`` - sum the retrievals into a single column.

``average`` - averages the retrievals into a single column.

If ``retrieve_mode`` is not specified, it will default to
``split and unique only``

**Example: retrieving fwr1 sequences**

.. code:: ipython3

    ddl.update_metadata(vdj, retrieve = 'fwr1')
    vdj




.. parsed-literal::

    Dandelion class object with n_obs = 2773 and n_contigs = 5609
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id', 'fwr1_VJ', 'fwr1_VDJ'
        layout: layout for 2773 vertices, layout for 1067 vertices
        graph: networkx graph of 2773 vertices, networkx graph of 1067 vertices 



Note the additional ``fwr1`` VDJ and VJ columns in the metadata slot.

By default, ``dandelion`` will not try to merge numerical columns as it
can create mixed dtype columns.

There is a new class function now that will try and retrieve frequently
used columns such as ``np1_length``, ``np2_length``:

.. code:: ipython3

    vdj.update_plus()
    vdj


.. parsed-literal::

    /Users/kt16/miniconda3/envs/dandelion/lib/python3.9/site-packages/numpy/core/fromnumeric.py:3440: RuntimeWarning: Mean of empty slice.
    /Users/kt16/miniconda3/envs/dandelion/lib/python3.9/site-packages/numpy/core/_methods.py:189: RuntimeWarning: invalid value encountered in double_scalars




.. parsed-literal::

    Dandelion class object with n_obs = 2773 and n_contigs = 5609
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id', 'fwr1_VJ', 'fwr1_VDJ', 'mu_count_VDJ', 'mu_count_VJ', 'mu_count', 'junction_length_VDJ', 'junction_length_VJ', 'junction_aa_length_VDJ', 'junction_aa_length_VJ', 'np1_length_VDJ', 'np1_length_VJ', 'np2_length_VDJ'
        layout: layout for 2773 vertices, layout for 1067 vertices
        graph: networkx graph of 2773 vertices, networkx graph of 1067 vertices 



concatenating multiple objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a simple function to concatenate (append) two or more
``Dandelion`` class, or ``pandas`` dataframes. Note that this operates
on the ``.data`` slot and not the ``.metadata`` slot.

.. code:: ipython3

    # for example, the original dandelion class has 2773 unique cell barcodes and 5609 contigs
    vdj




.. parsed-literal::

    Dandelion class object with n_obs = 2773 and n_contigs = 5609
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id', 'fwr1_VJ', 'fwr1_VDJ', 'mu_count_VDJ', 'mu_count_VJ', 'mu_count', 'junction_length_VDJ', 'junction_length_VJ', 'junction_aa_length_VDJ', 'junction_aa_length_VJ', 'np1_length_VDJ', 'np1_length_VJ', 'np2_length_VDJ'
        layout: layout for 2773 vertices, layout for 1067 vertices
        graph: networkx graph of 2773 vertices, networkx graph of 1067 vertices 



.. code:: ipython3

    # now it has 16827 (5609*3) contigs instead, and the metadata should also be properly populated
    vdj_concat = ddl.concat([vdj, vdj, vdj])
    vdj_concat




.. parsed-literal::

    Dandelion class object with n_obs = 2773 and n_contigs = 16827
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ'



``ddl.concat`` also lets you add in your custom prefixes/suffixes to
append to the sequence ids. If not provided, it will add ``-0``, ``-1``
etc. as a suffix if it detects that the sequence ids are not unique.

read/write
~~~~~~~~~~

``Dandelion`` class can be saved using ``.write_h5ddl`` and
``.write_pkl`` functions with accompanying compression methods.
``write_h5ddl`` primarily uses pandas ``to_hdf`` library and
``write_pkl`` just uses pickle. ``read_h5ddl`` and ``read_pkl``
functions will read the respective file formats accordingly.

.. code:: ipython3

    %time vdj.write_h5ddl('dandelion_results.h5ddl', complib = 'bzip2')


.. parsed-literal::

    <timed eval>:1: DeprecationWarning: write_h5 is a deprecated in 0.2.2 and will be removed in 0.4.0. write_h5ddl will be the recommended way to save.


.. parsed-literal::

    CPU times: user 4.5 s, sys: 154 ms, total: 4.65 s
    Wall time: 4.68 s


If you see any warnings above, it’s due to mix dtypes somewhere in the
object. So do some checking if you think it will interfere with
downstream usage.

.. code:: ipython3

    %time vdj_1 = ddl.read_h5ddl('dandelion_results.h5ddl')
    vdj_1


.. parsed-literal::

    <timed exec>:1: DeprecationWarning: read_h5 is a deprecated in 0.2.2 and will be removed in 0.4.0. read_h5ddl will be the recommended way to read.


.. parsed-literal::

    CPU times: user 1.34 s, sys: 148 ms, total: 1.48 s
    Wall time: 1.45 s




.. parsed-literal::

    Dandelion class object with n_obs = 2773 and n_contigs = 5609
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id', 'fwr1_VJ', 'fwr1_VDJ', 'mu_count_VDJ', 'mu_count_VJ', 'mu_count', 'junction_length_VDJ', 'junction_length_VJ', 'junction_aa_length_VDJ', 'junction_aa_length_VJ', 'np1_length_VDJ', 'np1_length_VJ', 'np2_length_VDJ'
        layout: layout for 2773 vertices, layout for 1067 vertices
        graph: networkx graph of 2773 vertices, networkx graph of 1067 vertices 



The read/write times using ``pickle`` can be situationally faster/slower
and file sizes can also be situationally smaller/larger (depending on
which compression is used).

.. code:: ipython3

    %time vdj.write_pkl('dandelion_results.pkl.gz')


.. parsed-literal::

    CPU times: user 8.5 s, sys: 54.3 ms, total: 8.55 s
    Wall time: 8.63 s


.. code:: ipython3

    %time vdj_2 = ddl.read_pkl('dandelion_results.pkl.gz')
    vdj_2


.. parsed-literal::

    CPU times: user 204 ms, sys: 21.3 ms, total: 226 ms
    Wall time: 236 ms




.. parsed-literal::

    Dandelion class object with n_obs = 2773 and n_contigs = 5609
        data: 'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'stop_codon', 'vj_in_frame', 'locus', 'junction_length', 'np1_length', 'np2_length', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'd_sequence_start', 'd_sequence_end', 'd_germline_start', 'd_germline_end', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', 'j_germline_end', 'v_score', 'v_identity', 'v_support', 'd_score', 'd_identity', 'd_support', 'j_score', 'j_identity', 'j_support', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'cdr1', 'cdr2', 'cdr3', 'cell_id', 'c_call', 'consensus_count', 'duplicate_count', 'v_call_10x', 'd_call_10x', 'j_call_10x', 'junction_10x', 'junction_10x_aa', 'v_call_genotyped', 'germline_alignment_d_mask', 'sample_id', 'j_support_igblastn', 'j_score_igblastn', 'j_call_igblastn', 'j_call_blastn', 'j_identity_blastn', 'j_alignment_length_blastn', 'j_number_of_mismatches_blastn', 'j_number_of_gap_openings_blastn', 'j_sequence_start_blastn', 'j_sequence_end_blastn', 'j_germline_start_blastn', 'j_germline_end_blastn', 'j_support_blastn', 'j_score_blastn', 'j_sequence_alignment_blastn', 'j_germline_alignment_blastn', 'cell_id_blastn', 'j_source', 'd_support_igblastn', 'd_score_igblastn', 'd_call_igblastn', 'd_call_blastn', 'd_identity_blastn', 'd_alignment_length_blastn', 'd_number_of_mismatches_blastn', 'd_number_of_gap_openings_blastn', 'd_sequence_start_blastn', 'd_sequence_end_blastn', 'd_germline_start_blastn', 'd_germline_end_blastn', 'd_support_blastn', 'd_score_blastn', 'd_sequence_alignment_blastn', 'd_germline_alignment_blastn', 'd_source', 'c_sequence_alignment', 'c_germline_alignment', 'c_sequence_start', 'c_sequence_end', 'c_score', 'c_identity', 'c_call_10x', 'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa', 'v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'mu_count', 'ambiguous', 'rearrangement_status', 'clone_id', 'changeo_clone_id'
        metadata: 'clone_id', 'clone_id_by_size', 'sample_id', 'locus_VDJ', 'locus_VJ', 'productive_VDJ', 'productive_VJ', 'v_call_genotyped_VDJ', 'd_call_VDJ', 'j_call_VDJ', 'v_call_genotyped_VJ', 'j_call_VJ', 'c_call_VDJ', 'c_call_VJ', 'junction_VDJ', 'junction_VJ', 'junction_aa_VDJ', 'junction_aa_VJ', 'v_call_genotyped_B_VDJ', 'd_call_B_VDJ', 'j_call_B_VDJ', 'v_call_genotyped_B_VJ', 'j_call_B_VJ', 'productive_B_VDJ', 'productive_B_VJ', 'duplicate_count_B_VDJ', 'duplicate_count_B_VJ', 'isotype', 'isotype_status', 'locus_status', 'chain_status', 'rearrangement_status_VDJ', 'rearrangement_status_VJ', 'changeo_clone_id', 'fwr1_VJ', 'fwr1_VDJ', 'mu_count_VDJ', 'mu_count_VJ', 'mu_count', 'junction_length_VDJ', 'junction_length_VJ', 'junction_aa_length_VDJ', 'junction_aa_length_VJ', 'np1_length_VDJ', 'np1_length_VJ', 'np2_length_VDJ'
        layout: layout for 2773 vertices, layout for 1067 vertices
        graph: networkx graph of 2773 vertices, networkx graph of 1067 vertices 



