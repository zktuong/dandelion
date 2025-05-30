#!/usr/bin/env python
import math
import os
import re
import sys

import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc

from anndata import AnnData
from changeo.Gene import getGene
from collections import defaultdict, Counter
from distance import hamming
from itertools import product
from scanpy import logging as logg
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
from subprocess import run
from tqdm import tqdm
from typing import Literal

from dandelion.tools._network import *
from dandelion.utilities._core import *
from dandelion.utilities._io import *
from dandelion.utilities._utilities import *


def find_clones(
    vdj_data: Dandelion | pd.DataFrame,
    identity: dict[str, float] | float = 0.85,
    key: str | None = None,
    by_alleles: bool = False,
    key_added: str | None = None,
    recalculate_length: bool = True,
    verbose: bool = True,
    **kwargs,
) -> Dandelion:
    """
    Find clones based on VDJ chain and VJ chain CDR3 junction hamming distance.

    Parameters
    ----------
    vdj_data : Dandelion | pd.DataFrame
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file
        after clones have been determined.
    identity : dict[str, float] | float, optional
        junction similarity parameter. Default 0.85. If provided as a dictionary, please use the following
        keys:'ig', 'tr-ab', 'tr-gd'.
    key : str | None, optional
        column name for performing clone clustering. `None` defaults to 'junction_aa'.
    by_alleles : bool, optional
        whether or not to collapse alleles to genes. `None` defaults to False.
    key_added : str | None, optional
        If specified, this will be the column name for clones. `None` defaults to 'clone_id'
    recalculate_length : bool, optional
        whether or not to re-calculate junction length, rather than rely on parsed assignment (which occasionally is
        wrong). Default is True
    verbose : bool, optional
        whether or not to print progress.
    **kwargs
        Additional arguments to pass to `Dandelion.update_metadata`.

    Returns
    -------
    Dandelion
        `Dandelion` object with clone_id annotated in `.data` slot and `.metadata` initialized.

    Raises
    ------
    ValueError
        if `key` not found in Dandelion.data.
    """
    start = logg.info("Finding clonotypes")
    pd.set_option("mode.chained_assignment", None)
    if isinstance(vdj_data, Dandelion):
        dat_ = load_data(vdj_data.data)
    else:
        dat_ = load_data(vdj_data)

    clone_key = key_added if key_added is not None else "clone_id"
    dat_[clone_key] = ""

    dat = dat_.copy()
    if "ambiguous" in dat_:
        dat = dat_[dat_["ambiguous"] == "F"].copy()

    locus_log = {"ig": "B", "tr-ab": "abT", "tr-gd": "gdT"}
    locus_dict1 = {"ig": ["IGH"], "tr-ab": ["TRB"], "tr-gd": ["TRD"]}
    locus_dict2 = {"ig": ["IGK", "IGL"], "tr-ab": ["TRA"], "tr-gd": ["TRG"]}
    DEFAULTIDENTITY = {"ig": 0.85, "tr-ab": 1, "tr-gd": 1}

    key_ = key if key is not None else "junction_aa"  # default

    if key_ not in dat.columns:
        raise ValueError("key {} not found in input table.".format(key_))

    locuses = ["ig", "tr-ab", "tr-gd"]

    # quick check
    for locus in locuses:
        locus_1 = locus_dict1[locus]
        locus_2 = locus_dict2[locus]

        dat_vj = dat[dat["locus"].isin(locus_2)].copy()
        dat_vdj = dat[dat["locus"].isin(locus_1)].copy()

        chain_check = check_chains(dat_vdj=dat_vdj, dat_vj=dat_vj)
        if all(~chain_check["All VDJ"]) and all(~chain_check["All VJ"]):
            locuses.remove(locus)

    if len(locuses) > 0:
        for locusx in locuses:
            locus_1 = locus_dict1[locusx]
            locus_2 = locus_dict2[locusx]
            if isinstance(identity, dict):
                if locusx not in identity:
                    identity.update({locusx: DEFAULTIDENTITY[locusx]})
                    warnings.warn(
                        UserWarning(
                            "Identity value for {} chains ".format(locusx)
                            + "not specified in provided dictionary. "
                            + "Defaulting to {} for {} chains.".format(
                                DEFAULTIDENTITY[locusx], locusx
                            )
                        )
                    )
                identity_ = identity[locusx]
            else:
                identity_ = identity

            dat_vj = dat[dat["locus"].isin(locus_2)].copy()
            dat_vdj = dat[dat["locus"].isin(locus_1)].copy()
            chain_check = check_chains(dat_vdj=dat_vdj, dat_vj=dat_vj)
            if dat_vdj.shape[0] > 0:
                vj_len_grp_vdj, seq_grp_vdj = group_sequences(
                    dat_vdj,
                    junction_key=key_,
                    recalculate_length=recalculate_length,
                    by_alleles=by_alleles,
                    locus=locusx,
                )
                cid_vdj = group_pairwise_hamming_distance(
                    clonotype_vj_len_group=vj_len_grp_vdj,
                    clonotype_sequence_group=seq_grp_vdj,
                    identity=identity_,
                    locus=locus_log[locusx],
                    chain="VDJ",
                    verbose=verbose,
                )
                clone_dict_vdj = rename_clonotype_ids(
                    clonotype_groups=cid_vdj,
                    cell_type=locus_log[locusx] + "_VDJ_",
                )
                # add it to the original dataframes
                dat_vdj[clone_key] = pd.Series(clone_dict_vdj)
                # dat[clone_key].update(pd.Series(dat_vdj[clone_key]))
                for i, row in dat_vdj.iterrows():
                    if i in dat.index:
                        dat.at[i, clone_key] = row[clone_key]
            if dat_vj.shape[0] > 0:
                vj_len_grp_vj, seq_grp_vj = group_sequences(
                    dat_vj,
                    junction_key=key_,
                    recalculate_length=recalculate_length,
                    by_alleles=by_alleles,
                    locus=locusx,
                )
                cid_vj = group_pairwise_hamming_distance(
                    clonotype_vj_len_group=vj_len_grp_vj,
                    clonotype_sequence_group=seq_grp_vj,
                    identity=identity_,
                    locus=locus_log[locusx],
                    chain="VJ",
                    verbose=verbose,
                )
                clone_dict_vj = rename_clonotype_ids(
                    clonotype_groups=cid_vj,
                    cell_type=locus_log[locusx] + "_VJ_",
                )
                refine_clone_assignment(
                    dat=dat,
                    clone_key=clone_key,
                    clone_dict_vj=clone_dict_vj,
                    verbose=verbose,
                )
                # dat_[clone_key].update(pd.Series(dat[clone_key]))
            for i, row in dat.iterrows():
                if i in dat_.index:
                    dat_.at[i, clone_key] = row[clone_key]

    # dat_[clone_key].replace('', 'unassigned')
    if os.path.isfile(str(vdj_data)):
        data_path = Path(vdj_data)
        write_airr(dat_, data_path.parent / (data_path.stem + "_clone.tsv"))

    logg.info(
        " finished",
        time=start,
        deep=(
            "Updated Dandelion object: \n"
            "   'data', contig-indexed AIRR table\n"
            "   'metadata', cell-indexed observations table\n"
        ),
    )
    if isinstance(vdj_data, Dandelion):
        if vdj_data.germline is not None:
            germline_ = vdj_data.germline
        else:
            germline_ = None
        if vdj_data.layout is not None:
            layout_ = vdj_data.layout
        else:
            layout_ = None
        if vdj_data.graph is not None:
            graph_ = vdj_data.graph
        else:
            graph_ = None
        if vdj_data.threshold is not None:
            threshold_ = vdj_data.threshold
        else:
            threshold_ = None
        if ("clone_id" in vdj_data.data.columns) and (key_added is None):
            # TODO: need to check the following bits if it works properly if only heavy chain tables are provided
            vdj_data.__init__(
                data=dat_,
                germline=germline_,
                layout=layout_,
                graph=graph_,
            )
            vdj_data.update_metadata(reinitialize=True, **kwargs)
        elif ("clone_id" in vdj_data.data.columns) and (key_added is not None):
            vdj_data.__init__(
                data=dat_,
                germline=germline_,
                layout=layout_,
                graph=graph_,
            )
            vdj_data.update_metadata(
                reinitialize=True,
                clone_key="clone_id",
                retrieve=clone_key,
                retrieve_mode="merge and unique only",
                **kwargs,
            )
        else:
            vdj_data.__init__(
                data=dat_,
                germline=germline_,
                layout=layout_,
                graph=graph_,
                clone_key=clone_key,
            )
            vdj_data.update_metadata(
                reinitialize=True, clone_key=clone_key, **kwargs
            )
        vdj_data.threshold = threshold_

    else:
        out = Dandelion(
            data=dat_,
            clone_key=clone_key,
            retrieve=clone_key,
            retrieve_mode="merge and unique only",
            **kwargs,
        )
        return out


def transfer(
    adata: AnnData,
    dandelion: Dandelion,
    expanded_only: bool = False,
    neighbors_key: str | None = None,
    rna_key: str | None = None,
    vdj_key: str | None = None,
    clone_key: str | None = None,
    collapse_nodes: bool = False,
    overwrite: bool | list[str] | str | None = None,
) -> None:
    """
    Transfer data in `Dandelion` slots to `AnnData` object, updating the `.obs`, `.uns`, `.obsm` and `.obsp`slots.

    Parameters
    ----------
    adata : AnnData
        `AnnData` object.
    dandelion : Dandelion
        `Dandelion` object.
    expanded_only : bool, optional
        Whether or not to transfer the embedding with all cells with BCR (False) or only for expanded clones (True).
    neighbors_key : str | None, optional
        key for 'neighbors' slot in `.uns`.
    rna_key : str | None, optional
        prefix for stashed RNA connectivities and distances.
    vdj_key : str | None, optional
        prefix for stashed VDJ connectivities and distances.
    clone_key : str | None, optional
        column name of clone/clonotype ids. Only used for integration with scirpy.
    collapse_nodes : bool, optional
        Whether or not to transfer a cell x cell or clone x clone connectivity matrix into `.uns`. Only used for
        integration with scirpy.
    overwrite : bool | list[str] | str | None, optional
        Whether or not to overwrite existing anndata columns. Specifying a string indicating column name or
        list of column names will overwrite that specific column(s).
    """
    start = logg.info("Transferring network")
    # always overwrite with whatever columns are in dandelion's metadata:
    for x in dandelion.metadata.columns:
        if x not in adata.obs.columns:
            adata.obs[x] = pd.Series(dandelion.metadata[x])
        elif overwrite is True:
            adata.obs[x] = pd.Series(dandelion.metadata[x])
        if type_check(dandelion.metadata, x):
            adata.obs[x] = adata.obs[x].replace(np.nan, "No_contig")
        if adata.obs[x].dtype == "bool":
            adata.obs[x] = [str(x) for x in adata.obs[x]]

    if (overwrite is not None) and (overwrite is not True):
        if not type(overwrite) is list:
            overwrite = [overwrite]
        for ow in overwrite:
            adata.obs[ow] = pd.Series(dandelion.metadata[ow])
            if type_check(dandelion.metadata, ow):
                adata.obs[ow] = adata.obs[ow].replace(np.nan, "No_contig")

    if dandelion.graph is not None:
        if expanded_only:
            G = dandelion.graph[1]
        else:
            G = dandelion.graph[0]
        logg.info("converting matrices")
        distances = nx.to_pandas_adjacency(
            G, dtype=np.float32, weight="weight", nonedge=np.nan
        )
        df_distances = distances.reindex(
            index=adata.obs_names, columns=adata.obs_names
        )
        # convert to connectivity
        distances = distances.apply(lambda x: np.maximum(1e-45, 1 / np.exp(x)))
        df_connectivities = distances.reindex(
            index=adata.obs_names, columns=adata.obs_names
        )

        df_connectivities = df_connectivities.values
        df_connectivities[np.isnan(df_connectivities)] = 0
        df_distances = df_distances.values
        df_distances[np.isnan(df_distances)] = 0

        df_connectivities_ = csr_matrix(df_connectivities, dtype=np.float32)
        df_distances = csr_matrix(df_distances, dtype=np.float32)

        logg.info("Updating anndata slots")
        if neighbors_key is None:
            neighbors_key = "neighbors"
        if neighbors_key not in adata.uns:
            skip_stash = True
            adata.uns[neighbors_key] = {}

        rna_neighbors_key = "rna_" + neighbors_key
        vdj_neighbors_key = "vdj_" + neighbors_key
        if rna_neighbors_key not in adata.uns:
            adata.uns[rna_neighbors_key] = adata.uns[neighbors_key].copy()

        if rna_key is None:
            r_connectivities_key = "rna_connectivities"
            r_distances_key = "rna_distances"
        else:
            r_connectivities_key = rna_key + "_connectivitites"
            r_distances_key = rna_key + "_distances"

        if vdj_key is None:
            b_connectivities_key = "vdj_connectivities"
            b_distances_key = "vdj_distances"
        else:
            b_connectivities_key = vdj_key + "_connectivitites"
            b_distances_key = vdj_key + "_distances"

        # stash_rna_connectivities:
        if r_connectivities_key not in adata.obsp:
            if "skip_stash" not in locals():
                try:
                    adata.obsp[r_connectivities_key] = adata.obsp[
                        "connectivities"
                    ].copy()
                    adata.obsp[r_distances_key] = adata.obsp["distances"].copy()
                except:
                    adata.obsp[r_connectivities_key] = adata.uns[neighbors_key][
                        "connectivities"
                    ]
                    adata.obsp[r_distances_key] = adata.uns[neighbors_key][
                        "distances"
                    ]

        # always overwrite the bcr slots
        adata.obsp["connectivities"] = df_connectivities_.copy()
        adata.obsp["distances"] = df_distances.copy()
        adata.obsp[b_connectivities_key] = adata.obsp["connectivities"].copy()
        adata.obsp[b_distances_key] = adata.obsp["distances"].copy()

        # create the dictionary that will enable the use of scirpy's plotting.
        clonekey = clone_key if clone_key is not None else "clone_id"

        if not collapse_nodes:
            df_connectivities_[df_connectivities_.nonzero()] = 1
            cell_indices = {
                str(i): np.array([k])
                for i, k in zip(range(0, len(adata.obs_names)), adata.obs_names)
            }
        else:
            cell_indices = Tree()
            for x, y in adata.obs[clonekey].items():
                if y not in [
                    "",
                    "unassigned",
                    np.nan,
                    "NaN",
                    "NA",
                    "nan",
                    "None",
                    None,
                    "none",
                ]:
                    cell_indices[y][x].value = 1
            cell_indices = {
                str(x): np.array(list(r))
                for x, r in zip(
                    range(0, len(cell_indices)), cell_indices.values()
                )
            }
            df_connectivities_ = np.zeros(
                [len(cell_indices), len(cell_indices)]
            )
            np.fill_diagonal(df_connectivities_, 1)
            df_connectivities_ = csr_matrix(df_connectivities_)

        adata.uns[clonekey] = {
            # this is a symmetrical, pairwise, sparse distance matrix of clonotypes
            # the matrix is offset by 1, i.e. 0 = no connection, 1 = distance 0
            "distances": df_connectivities_,
            # '0' refers to the row/col index in the `distances` matrix
            # (numeric index, but needs to be strbecause of h5py)
            # np.array(["cell1", "cell2"]) points to the rows in `adata.obs`
            "cell_indices": cell_indices,
        }

    tmp = adata.obs.copy()
    if dandelion.graph is not None:
        if dandelion.layout is not None:
            if expanded_only:
                coord = pd.DataFrame.from_dict(
                    dandelion.layout[1], orient="index"
                )
            else:
                coord = pd.DataFrame.from_dict(
                    dandelion.layout[0], orient="index"
                )
            for x in coord.columns:
                tmp[x] = coord[x]

            X_vdj = np.array(tmp[[0, 1]], dtype=np.float32)
            adata.obsm["X_vdj"] = X_vdj

        logg.info(
            " finished",
            time=start,
            deep=(
                "updated `.obs` with `.metadata`\n"
                "added to `.uns['"
                + neighbors_key
                + "']` and `.uns['"
                + clonekey
                + "']`\n"
                "and `.obsp`\n"
                "   'distances', clonotype-weighted adjacency matrix\n"
                "   'connectivities', clonotype-weighted adjacency matrix"
            ),
        )
    else:
        logg.info(
            " finished", time=start, deep=("updated `.obs` with `.metadata`\n")
        )


def define_clones(
    vdj_data: Dandelion | pd.DataFrame | str,
    dist: float | None = None,
    action: Literal["first", "set"] = "set",
    model: Literal[
        "ham",
        "aa",
        "hh_s1f",
        "hh_s5f",
        "mk_rs1nf",
        "mk_rs5nf",
        "hs1f_compat",
        "m1n_compat",
    ] = "ham",
    norm: Literal["len", "mut", "none"] = "len",
    doublets: Literal["drop", "count"] = "drop",
    fileformat: Literal["changeo", "airr"] = "airr",
    ncpu: int | None = None,
    outFilePrefix: int | None = None,
    key_added: int | None = None,
    out_dir: Path | str | None = None,
    additional_args: list[str] = [],
) -> Dandelion:
    """
    Find clones using changeo's `DefineClones.py <https://changeo.readthedocs.io/en/stable/tools/DefineClones.html>`__.

    Only callable for BCR data at the moment.

    Parameters
    ----------
    vdj_data : Dandelion | pd.DataFrame | str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after
        clones have been determined.
    dist : float | None, optional
        The distance threshold for clonal grouping. If None, the value will be retrieved from the Dandelion class
        `.threshold` slot.
    action : Literal["first", "set"], optional
        Specifies how to handle multiple V(D)J assignments for initial grouping. Default is 'set'.
        The “first” action will use only the first gene listed. The “set” action will use all gene assignments and
        construct a larger gene grouping composed of any sequences sharing an assignment or linked to another sequence
        by a common assignment (similar to single-linkage).
    model : Literal["ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "mk_rs5nf", "hs1f_compat", "m1n_compat", ], optional
        Specifies which substitution model to use for calculating distance between sequences. Default is 'ham'.
        The “ham” model is nucleotide Hamming distance and “aa” is amino acid Hamming distance. The “hh_s1f” and
        “hh_s5f” models are human specific single nucleotide and 5-mer content models, respectively, from Yaari et al,
        2013. The “mk_rs1nf” and “mk_rs5nf” models are mouse specific single nucleotide and 5-mer content models,
        respectively, from Cui et al, 2016. The “m1n_compat” and “hs1f_compat” models are deprecated models provided
        backwards compatibility with the “m1n” and “hs1f” models in Change-O v0.3.3 and SHazaM v0.1.4. Both 5-mer
        models should be considered experimental.
    norm : Literal["len", "mut", "none"], optional
        Specifies how to normalize distances. Default is 'len'. 'none' (do not normalize), 'len' (normalize by length),
        or 'mut' (normalize by number of mutations between sequences).
    doublets : Literal["drop", "count"], optional
        Option to control behaviour when dealing with heavy chain 'doublets'. Default is 'drop'. 'drop' will filter out
        the doublets while 'count' will retain only the highest umi count contig.
    fileformat : Literal["changeo", "airr"], optional
        Format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
    ncpu : int | None, optional
        Number of cpus for parallelization. Default is 1, no parallelization.
    outFilePrefix : int | None, optional
        If specified, the out file name will have this prefix. `None` defaults to 'dandelion_define_clones'
    key_added : int | None, optional
        Column name to add for define_clones.
    out_dir : Path | str | None, optional
        If specified, the files will be written to this directory.
    additional_args : list[str], optional
        Additional arguments to pass to `DefineClones.py`.

    Returns
    -------
    Dandelion
        `Dandelion` object with clone_id annotated in `.data` slot and `.metadata` initialized.

    Raises
    ------
    ValueError
        if .threshold not found in `Dandelion`.
    """
    start = logg.info("Finding clones")
    if ncpu is None:
        nproc = 1
    else:
        nproc = ncpu

    if key_added is None:
        clone_key = "clone_id"
    else:
        clone_key = key_added

    if isinstance(vdj_data, Dandelion):
        dat_ = load_data(vdj_data.data)
    else:
        dat_ = load_data(vdj_data)
    if "ambiguous" in dat_:
        dat = dat_[dat_["ambiguous"] == "F"].copy()
    else:
        dat = dat_.copy()
    dat_h = dat[dat["locus"] == "IGH"]
    dat_l = dat[dat["locus"].isin(["IGK", "IGL"])]

    if os.path.isfile(str(vdj_data)):
        vdj_path = Path(vdj_data)
        tmpFolder = vdj_path.parent / "tmp"
        outFolder = vdj_path.parent
    elif out_dir is not None:
        vdj_path = Path(out_dir)
        tmpFolder = vdj_path / "tmp"
        outFolder = vdj_path
    else:
        import tempfile

        outFolder = Path(tempfile.TemporaryDirectory().name)
        tmpFolder = outFolder / "tmp"

    for _ in [outFolder, tmpFolder]:
        _.mkdir(parents=True, exist_ok=True)

    if "vdj_path" in locals():
        h_file1 = tmpFolder / (vdj_path.stem + "_heavy-clone.tsv")
        h_file2 = outFolder / (vdj_path.stem + "_heavy-clone.tsv")
        l_file = tmpFolder / (vdj_path.stem + "_light.tsv")
        outfile = outFolder / (vdj_path.stem + "_clone.tsv")
    else:
        out_FilePrefix = (
            "dandelion_define_clones"
            if outFilePrefix is None
            else outFilePrefix
        )
        h_file1 = tmpFolder / (out_FilePrefix + "_heavy-clone.tsv")
        h_file2 = outFolder / (out_FilePrefix + "_heavy-clone.tsv")
        l_file = tmpFolder / (out_FilePrefix + "_light.tsv")
        outfile = outFolder / (out_FilePrefix + "_clone.tsv")
    write_airr(dat_h, h_file1)
    write_airr(dat_l, l_file)
    v_field = (
        "v_call_genotyped" if "v_call_genotyped" in dat.columns else "v_call"
    )
    if dist is None:
        if isinstance(vdj_data, Dandelion):
            if vdj_data.threshold is not None:
                dist_ = vdj_data.threshold
            else:
                raise ValueError(
                    "Threshold value in Dandelion object is None. Please run calculate_threshold first"
                )
        else:
            raise ValueError(
                "Distance value is None. Please provide a distance value (float)"
            )
    else:
        dist_ = dist

    cmd = [
        "DefineClones.py",
        "-d",
        str(h_file1),
        "-o",
        str(h_file2),
        "--act",
        action,
        "--model",
        model,
        "--norm",
        norm,
        "--dist",
        str(dist_),
        "--nproc",
        str(nproc),
        "--vf",
        v_field,
    ]
    cmd = cmd + additional_args

    def clusterLinkage(cell_series, group_series):
        """
        Return a dictionary of {cell_id : cluster_id}.

        that identifies clusters of cells by analyzing their shared
        features (group_series) using single linkage.

        Arguments:
        cell_series (iter): iter of cell ids.
        group_series (iter): iter of group ids.

        Returns:
        dict:  dictionary of {cell_id : cluster_id}.

        """
        # assign initial clusters
        # initial_dict = {cluster1: [cell1], cluster2: [cell1]}
        initial_dict = {}
        for cell, group in zip(cell_series, group_series):
            try:
                initial_dict[group].append(cell)
            except KeyError:
                initial_dict[group] = [cell]

        # naive single linkage clustering (ON^2 best case, ON^3 worst case) ...ie for cells with multiple light chains
        # cluster_dict = {cluster1: [cell1, cell2]}, 2 cells belong in same group if they share 1 light chain
        while True:
            cluster_dict = {}
            for i, group in enumerate(initial_dict.keys()):
                cluster_dict[i] = initial_dict[group]
                for cluster in cluster_dict:
                    # if initial_dict[group] and cluster_dict[cluster] share common cells, add initial_dict[group] to
                    # cluster
                    if cluster != i and any(
                        cell in initial_dict[group]
                        for cell in cluster_dict[cluster]
                    ):
                        cluster_dict[cluster] = (
                            cluster_dict[cluster] + initial_dict[group]
                        )
                        del cluster_dict[i]
                        break
            # break if clusters stop changing, otherwise restart
            if len(cluster_dict.keys()) == len(initial_dict.keys()):
                break
            else:
                initial_dict = cluster_dict.copy()

        # invert cluster_dict for return
        assign_dict = {
            cell: k for k, v in cluster_dict.items() for cell in set(v)
        }

        return assign_dict

    # TODO: might need to remove this function to drop requirement to maintain this as a dependency internally
    def _lightCluster(heavy_file, light_file, out_file, doublets, fileformat):
        """
        Split heavy chain clones based on light chains.

        Arguments:
        heavy_file (str): heavy chain input file.
        light_file (str): light chain input file.
        out_file (str): heavy chain output file.
        doublets (str): method for handling multiple heavy chains per cell. one of 'drop' or 'count'.
        format (str): file format. one of 'changeo' or 'airr'.
        """
        # Set column names
        if fileformat == "changeo":
            cell_id = "cell_id"
            clone_id = "clone_id"
            v_call = "v_call"
            j_call = "j_call"
            junction_length = "junction_length"
            umi_count = "umicount"
        elif fileformat == "airr":
            cell_id = "cell_id"
            clone_id = "clone_id"
            v_call = "v_call"
            j_call = "j_call"
            junction_length = "junction_length"
            umi_count = "umi_count"
        else:
            sys.exit("Invalid format %s" % fileformat)

        # read in heavy and light DFs
        heavy_df = pd.read_csv(
            heavy_file, dtype="object", na_values=["", "None", "NA"], sep="\t"
        )
        light_df = pd.read_csv(
            light_file, dtype="object", na_values=["", "None", "NA"], sep="\t"
        )

        # column checking
        expected_heavy_columns = [
            cell_id,
            clone_id,
            v_call,
            j_call,
            junction_length,
            umi_count,
        ]
        if set(expected_heavy_columns).issubset(heavy_df.columns) is False:
            raise ValueError(
                "Missing one or more columns in heavy chain file: "
                + ", ".join(expected_heavy_columns)
            )
        expected_light_columns = [
            cell_id,
            v_call,
            j_call,
            junction_length,
            umi_count,
        ]
        if set(expected_light_columns).issubset(light_df.columns) is False:
            raise ValueError(
                "Missing one or more columns in light chain file: "
                + ", ".join(expected_light_columns)
            )

        # Fix types
        try:
            heavy_df[junction_length] = heavy_df[junction_length].astype("int")
            light_df[junction_length] = light_df[junction_length].astype("int")
        except:
            heavy_df[junction_length] = heavy_df[junction_length].replace(
                np.nan, pd.NA
            )
            light_df[junction_length] = light_df[junction_length].replace(
                np.nan, pd.NA
            )
            heavy_df[junction_length] = heavy_df[junction_length].astype(
                "Int64"
            )
            light_df[junction_length] = light_df[junction_length].astype(
                "Int64"
            )

        # filter multiple heavy chains
        if doublets == "drop":
            heavy_df = heavy_df.drop_duplicates(cell_id, keep=False)
            if heavy_df.empty is True:
                raise ValueError(
                    "Empty heavy chain data, after doublets drop. Are you combining experiments "
                    "in a single file? If so, split your data into multiple files."
                )
        elif doublets == "count":
            heavy_df[umi_count] = heavy_df[umi_count].astype("int")
            heavy_df = heavy_df.groupby(cell_id, sort=False).apply(
                lambda x: x.nlargest(1, umi_count)
            )

        # transfer clone IDs from heavy chain df to light chain df
        clone_dict = {
            v[cell_id]: v[clone_id]
            for k, v in heavy_df[[clone_id, cell_id]].T.to_dict().items()
        }
        light_df = light_df.loc[
            light_df[cell_id].apply(lambda x: x in clone_dict.keys()),
        ]
        light_df[clone_id] = light_df.apply(
            lambda row: clone_dict[row[cell_id]], axis=1
        )

        # generate a "cluster_dict" of CELL:CLONE dictionary from light df  (TODO: use receptor object V/J gene names)
        cluster_dict = clusterLinkage(
            light_df[cell_id],
            light_df.apply(
                lambda row: getGene(row[v_call])
                + ","
                + getGene(row[j_call])
                + ","
                + str(row[junction_length])
                + ","
                + row[clone_id],
                axis=1,
            ),
        )

        # add assignments to heavy_df
        heavy_df = heavy_df.loc[
            heavy_df[cell_id].apply(lambda x: x in cluster_dict.keys()), :
        ]
        heavy_df[clone_id] = (
            heavy_df[clone_id]
            + "_"
            + heavy_df.apply(
                lambda row: str(cluster_dict[row[cell_id]]), axis=1
            )
        )

        # write heavy chains
        write_airr(heavy_df, out_file)
        return (heavy_df, light_df)

    logg.info("Running command: %s\n" % (" ".join(cmd)))
    run(cmd)

    h_df, l_df = _lightCluster(
        h_file2, l_file, outfile, doublets=doublets, fileformat=fileformat
    )

    h_df = load_data(h_df)
    # create a dictionary for cell_id : clone_id from h_df
    linked_clones = dict(zip(h_df["cell_id"], h_df["clone_id"]))

    # create a clone_reference
    clone_ref = list(set(h_df["clone_id"]))
    clone_ref = [c.split("_")[1] if c is not np.nan else c for c in clone_ref]
    l_df = load_data(l_df)

    for x in l_df.index:
        if l_df.loc[x, "clone_id"] in clone_ref:
            l_df.at[x, "clone_id"] = linked_clones[l_df.loc[x, "cell_id"]]
        else:
            try:
                l_df.at[x, "clone_id"] = l_df.loc[x, "cell_id"] + "_notlinked"
            except:
                pass

    cloned_ = pd.concat([h_df, l_df])
    # transfer the new clone_id to the heavy + light file
    dat_[str(clone_key)] = pd.Series(cloned_["clone_id"])
    dat_[str(clone_key)] = dat_[str(clone_key)].fillna("")
    if isinstance(vdj_data, Dandelion):
        germline_ = vdj_data.germline if vdj_data.germline is not None else None
        layout_ = vdj_data.layout if vdj_data.layout is not None else None
        graph_ = vdj_data.graph if vdj_data.graph is not None else None
        threshold_ = (
            vdj_data.threshold if vdj_data.threshold is not None else None
        )
        if ("clone_id" in vdj_data.data) and (clone_key is not None):
            vdj_data.__init__(
                data=dat_,
                germline=germline_,
                layout=layout_,
                graph=graph_,
                initialize=True,
                retrieve=clone_key,
                retrieve_mode="merge and unique only",
            )
        elif ("clone_id" not in vdj_data.data) and (clone_key is not None):
            vdj_data.__init__(
                data=dat_,
                germline=germline_,
                layout=layout_,
                graph=graph_,
                initialize=True,
                clone_key=clone_key,
                retrieve=clone_key,
                retrieve_mode="merge and unique only",
            )
        else:
            vdj_data.__init__(
                data=dat_,
                germline=germline_,
                layout=layout_,
                graph=graph_,
                initialize=True,
                clone_key=clone_key,
            )
        vdj_data.threshold = threshold_
    else:
        if ("clone_id" in dat_.columns) and (clone_key is not None):
            out = Dandelion(
                data=dat_,
                retrieve=clone_key,
                retrieve_mode="merge and unique only",
            )
        elif ("clone_id" not in dat_.columns) and (clone_key is not None):
            out = Dandelion(data=dat_, clone_key=clone_key)
        else:
            out = Dandelion(data=dat_)
        return out
    logg.info(
        " finished",
        time=start,
        deep=(
            "Updated Dandelion object: \n"
            "   'data', contig-indexed AIRR table\n"
            "   'metadata', cell-indexed observations table\n"
        ),
    )


def tabuluate_clone_sizes(
    metadata_: pd.DataFrame, clonesize_dict: dict, clonekey: str
) -> pd.Series:
    """Tabulate clone sizes."""
    return pd.Series(
        dict(
            zip(
                metadata_.index,
                [
                    str(y) if pd.notnull(y) else str(0)
                    for y in [
                        (
                            sorted(
                                list(
                                    {clonesize_dict[c_] for c_ in c.split("|")}
                                ),
                                key=lambda x: (
                                    int(x.split(">= ")[1])
                                    if type(x) is str
                                    else int(x)
                                ),
                                reverse=True,
                            )[0]
                            if "|" in c
                            else clonesize_dict[c]
                        )
                        for c in metadata_[str(clonekey)]
                    ]
                ],
            )
        )
    )


def clone_size(
    vdj_data: Dandelion,
    max_size: int | None = None,
    clone_key: str | None = None,
    key_added: str | None = None,
) -> None:
    """
    Quantify size of clones.

    Parameters
    ----------
    vdj_data : Dandelion
        `Dandelion` object
    max_size : int | None, optional
        The maximum size before value gets clipped. If None, the value will be returned as a numerical value.
    clone_key : str | None, optional
        Column name specifying the clone_id column in metadata.
    key_added : str | None, optional
        column name where clone size is tabulated into.
    """
    start = logg.info("Quantifying clone sizes")

    metadata_ = vdj_data.metadata.copy()

    clonekey = "clone_id" if clone_key is None else clone_key

    tmp = metadata_[str(clonekey)].str.split("|", expand=True).stack()
    tmp = tmp.reset_index(drop=False)
    tmp.columns = ["cell_id", "tmp", str(clonekey)]

    clonesize = tmp[str(clonekey)].value_counts()
    if "None" in clonesize.index:
        clonesize.drop("None", inplace=True)

    if max_size is not None:
        clonesize_ = clonesize.astype("object")
        for i in clonesize.index:
            if clonesize.loc[i] >= max_size:
                clonesize_.at[i] = ">= " + str(max_size)
        clonesize_ = clonesize_.astype("category")
    else:
        clonesize_ = clonesize.copy()

    clonesize_dict = dict(clonesize_)
    clonesize_dict.update({"None": np.nan})

    col_key, col_key_suffix = "", ""
    col_key = str(clonekey) if key_added is None else key_added
    col_key_suffix = (
        "_size" if max_size is None else "_size_max_" + str(max_size)
    )

    vdj_data.metadata[col_key + col_key_suffix] = tabuluate_clone_sizes(
        metadata_, clonesize_dict, clonekey
    )
    if max_size is not None:
        vdj_data.metadata[col_key + col_key_suffix] = vdj_data.metadata[
            col_key + col_key_suffix
        ].astype("category")
    else:
        try:
            vdj_data.metadata[col_key + col_key_suffix] = [
                float(x) for x in vdj_data.metadata[col_key + col_key_suffix]
            ]
        except ValueError:
            # this happens if there are multiple clonotypes associated to the cell
            pass

    logg.info(
        " finished",
        time=start,
        deep=(
            "Updated Dandelion object: \n"
            "   'metadata', cell-indexed clone table"
        ),
    )


def clone_overlap(
    vdj_data: Dandelion | AnnData,
    groupby: str,
    min_clone_size: int | None = None,
    weighted_overlap: bool = False,
    clone_key: str | None = None,
) -> pd.DataFrame:
    """
    A function to tabulate clonal overlap for input as a circos-style plot.

    Parameters
    ----------
    vdj_data : Dandelion | AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        column name in obs/metadata for collapsing to columns in the clone_id x groupby data frame.
    min_clone_size : int | None, optional
        minimum size of clone for plotting connections. Defaults to 2 if left as None.
    weighted_overlap : bool, optional
        if True, instead of collapsing to overlap to binary, overlap will be returned as the number of cells.
        In the future, there will be the option to use something like a jaccard index.
    clone_key : str | None, optional
        column name for clones. `None` defaults to 'clone_id'.

    Returns
    -------
    pd.DataFrame
        clone_id x groupby overlap :class:`pandas.core.frame.DataFrame'.

    Raises
    ------
    ValueError
        if min_clone_size is 0.
    """
    start = logg.info("Finding clones")
    if isinstance(vdj_data, Dandelion):
        data = vdj_data.metadata.copy()
    elif isinstance(vdj_data, AnnData):
        data = vdj_data.obs.copy()

    if min_clone_size is None:
        min_size = 2
    else:
        min_size = int(min_clone_size)

    if clone_key is None:
        clone_ = "clone_id"
    else:
        clone_ = clone_key

    # get rid of problematic rows that appear because of category conversion?
    allgroups = list(data[groupby].unique())
    data = data[
        ~(
            data[clone_].isin(
                [np.nan, "nan", "NaN", "No_contig", "unassigned", "None", None]
            )
        )
    ]

    # prepare a summary table
    datc_ = data[clone_].str.split("|", expand=True).stack()
    datc_ = pd.DataFrame(datc_)
    datc_.reset_index(drop=False, inplace=True)
    datc_.columns = ["cell_id", "tmp", clone_]
    datc_.drop("tmp", inplace=True, axis=1)
    datc_ = datc_[
        ~(
            datc_[clone_].isin(
                [
                    "",
                    np.nan,
                    "nan",
                    "NaN",
                    "No_contig",
                    "unassigned",
                    "None",
                    None,
                ]
            )
        )
    ]
    dictg_ = dict(data[groupby])
    datc_[groupby] = [dictg_[l] for l in datc_["cell_id"]]

    overlap = pd.crosstab(datc_[clone_], datc_[groupby])
    for x in allgroups:
        if x not in overlap:
            overlap[x] = 0

    if min_size == 0:
        raise ValueError("min_size must be greater than 0.")
    if not weighted_overlap:
        if min_size > 2:
            overlap[overlap < min_size] = 0
            overlap[overlap >= min_size] = 1
        elif min_size == 2:
            overlap[overlap >= min_size] = 1

    overlap.index.name = None
    overlap.columns.name = None

    if isinstance(vdj_data, AnnData):
        vdj_data.uns["clone_overlap"] = overlap.copy()
        logg.info(
            " finished",
            time=start,
            deep=("Updated AnnData: \n" "   'uns', clone overlap table"),
        )
    else:
        return overlap


def clustering(
    distance_dict: dict, threshold: float, sequences_dict: dict[str, str]
) -> dict:
    """Clustering the sequences."""
    out_dict = {}
    # find out the unique indices in this subset
    i_unique = list(set(flatten(distance_dict)))
    # for every pair of i1,i2 is their dictance smaller than the thresholdeshold?
    i_pair_d = {
        (i1, i2): (
            distance_dict[(i1, i2)] <= threshold
            if (i1, i2) in distance_dict
            else False
        )
        for i1, i2 in product(i_unique, repeat=2)
    }
    i_pair_d.update(
        {
            (i2, i1): (
                distance_dict[(i2, i1)] <= threshold
                if (i2, i1) in distance_dict
                else False
            )
            for i1, i2 in product(i_unique, repeat=2)
        }
    )
    # so which indices should not be part of a clone?
    canbetogether = defaultdict(list)
    for ii1, ii2 in product(i_unique, repeat=2):
        if i_pair_d[(ii1, ii2)] or i_pair_d[(ii2, ii1)]:
            if (ii1, ii2) in distance_dict:
                canbetogether[ii1].append((ii1, ii2))
                canbetogether[ii2].append((ii1, ii2))
            elif (ii2, ii1) in distance_dict:
                canbetogether[ii2].append((ii2, ii1))
                canbetogether[ii1].append((ii2, ii1))
        else:
            if (ii1, ii2) or (ii2, ii1) in distance_dict:
                canbetogether[ii1].append(())
                canbetogether[ii2].append(())
    for x in canbetogether:
        canbetogether[x] = list({y for y in canbetogether[x] if len(y) > 0})
    # convert the indices to sequences
    for x in canbetogether:
        if len(canbetogether[x]) > 0:
            out_dict[sequences_dict[x]] = tuple(
                sorted(
                    set(
                        list(
                            [
                                sequences_dict[y]
                                for y in flatten(canbetogether[x])
                            ]
                        )
                        + [sequences_dict[x]]
                    )
                )
            )
        else:
            out_dict[sequences_dict[x]] = tuple([sequences_dict[x]])
    return out_dict


def productive_ratio(
    adata: AnnData,
    vdj: Dandelion,
    groupby: str,
    groups: list[str] | None = None,
    locus: Literal["TRB", "TRA", "TRD", "TRG", "IGH", "IGK", "IGL"] = "TRB",
):
    """
    Compute the cell-level productive/non-productive contig ratio.

    Only the contig with the highest umi count in a cell will be used for this
    tabulation.

    Returns inplace `AnnData` with `.uns['productive_ratio']`.

    Parameters
    ----------
    adata : AnnData
        AnnData object holding the cell level metadata (`.obs`).
    vdj : Dandelion
        Dandelion object holding the repertoire data (`.data`).
    groupby : str
        Name of column in `AnnData.obs` to return the row tabulations.
    groups : list[str] | None, optional
        Optional list of categories to return.
    locus : Literal["TRB", "TRA", "TRD", "TRG", "IGH", "IGK", "IGL"], optional
        One of the accepted locuses to perform the tabulation
    """
    start = logg.info("Tabulating productive ratio")
    vdjx = vdj[(vdj.data.cell_id.isin(adata.obs_names))].copy()
    if "ambiguous" in vdjx.data:
        tmp = vdjx[
            (vdjx.data.locus == locus) & (vdjx.data.ambiguous == "F")
        ].copy()
    else:
        tmp = vdjx[(vdjx.data.locus == locus)].copy()

    if groups is None:
        if is_categorical(adata.obs[groupby]):
            groups = list(adata.obs[groupby].cat.categories)
        else:
            groups = list(set(adata.obs[groupby]))
    df = tmp.data.drop_duplicates(subset="cell_id")
    dict_df = dict(zip(df.cell_id, df.productive))
    res = pd.DataFrame(
        columns=["productive", "non-productive", "total"],
        index=groups,
    )

    adata.obs[locus + "_productive"] = pd.Series(dict_df)
    for i in range(res.shape[0]):
        cell = res.index[i]
        res.loc[cell, "total"] = sum(adata.obs[groupby] == cell)
        if res.loc[cell, "total"] > 0:
            res.loc[cell, "productive"] = (
                sum(
                    adata.obs.loc[
                        adata.obs[groupby] == cell, locus + "_productive"
                    ].isin(["T"])
                )
                / res.loc[cell, "total"]
                * 100
            )
            res.loc[cell, "non-productive"] = (
                sum(
                    adata.obs.loc[
                        adata.obs[groupby] == cell, locus + "_productive"
                    ].isin(
                        [
                            "F",
                        ]
                    )
                )
                / res.loc[cell, "total"]
                * 100
            )
    res[groupby] = res.index
    res["productive+non-productive"] = res["productive"] + res["non-productive"]
    out = {"results": res, "locus": locus, "groupby": groupby}
    adata.uns["productive_ratio"] = out
    logg.info(
        " finished",
        time=start,
        deep=("Updated AnnData: \n" "   'uns', productive_ratio"),
    )


def vj_usage_pca(
    adata: AnnData,
    groupby: str,
    min_size: int = 20,
    mode: Literal["B", "abT", "gdT"] = "abT",
    transfer_mapping=None,
    n_comps: int = 30,
    groups: list[str] | None = None,
    allowed_chain_status: list[str] | None = [
        "Single pair",
        "Extra pair",
        "Extra pair-exception",
        "Orphan VDJ-exception",
    ],
    verbose=False,
    **kwargs,
) -> AnnData:
    """
    Extract productive V/J gene usage from single cell data and compute PCA.

    Parameters
    ----------
    adata : AnnData
        AnnData object holding the cell level metadata with Dandelion VDJ info transferred.
    groupby : str
        Column name in `adata.obs` to groupby as observations for PCA.
    min_size : int, optional
        Minimum cell size numbers to keep for computing the final matrix. Defaults to 20.
    mode : Literal["B", "abT", "gdT"], optional
        Mode for extract the V/J genes.
    transfer_mapping : None, optional
        If provided, the columns will be mapped to the output AnnData from the original AnnData.
    n_comps : int, optional
        Number of principal components to compute. Defaults to 30.
    groups : list[str] | None, optional
        If provided, only the following groups/categories will be used for computing the PCA.
    allowed_chain_status : list[str] | None, optional
        If provided, only the ones in this list are kept from the `chain_status` column.
        Defaults to ["Single pair", "Extra pair", "Extra pair-exception", "Orphan VDJ-exception"].
    verbose : bool, optional
        Whether to display progress
    **kwargs
        Additional keyword arguments passed to `scanpy.pp.pca`.

    Returns
    -------
    AnnData
        AnnData object with obs as groups and V/J genes as features.
    """
    start = logg.info("Computing PCA for V/J gene usage")
    if allowed_chain_status is not None:
        adata_ = adata[
            adata.obs["chain_status"].isin(allowed_chain_status)
        ].copy()

    if groups is not None:
        adata_ = adata_[adata_.obs[group].isin(groups)].copy()

    if "v_call_genotyped_VDJ" in adata_.obs:
        v_call = "v_call_genotyped_"
    else:
        v_call = "v_call_"

    # prep data
    adata_.obs[v_call + mode + "_VJ_main"] = [
        x.split("|")[0] for x in adata_.obs[v_call + mode + "_VJ"]
    ]
    adata_.obs["j_call_" + mode + "_VJ_main"] = [
        x.split("|")[0] for x in adata_.obs["j_call_" + mode + "_VJ"]
    ]
    adata_.obs[v_call + mode + "_VDJ_main"] = [
        x.split("|")[0] for x in adata_.obs[v_call + mode + "_VDJ"]
    ]
    adata_.obs["j_call_" + mode + "_VDJ_main"] = [
        x.split("|")[0] for x in adata_.obs["j_call_" + mode + "_VDJ"]
    ]

    df1 = pd.DataFrame(
        {
            groupby: Counter(adata_.obs[groupby]).keys(),
            "cellcount": Counter(adata_.obs[groupby]).values(),
        }
    )
    new_list = df1.loc[df1["cellcount"] >= min_size, groupby]

    vj_v_list = [
        x
        for x in list(set(adata_.obs[v_call + mode + "_VJ_main"]))
        if x not in ["None", "No_contig"]
    ]
    vj_j_list = [
        x
        for x in list(set(adata_.obs["j_call_" + mode + "_VJ_main"]))
        if x not in ["None", "No_contig"]
    ]
    vdj_v_list = [
        x
        for x in list(set(adata_.obs[v_call + mode + "_VDJ_main"]))
        if x not in ["None", "No_contig"]
    ]
    vdj_j_list = [
        x
        for x in list(set(adata_.obs["j_call_" + mode + "_VDJ_main"]))
        if x not in ["None", "No_contig"]
    ]

    new_list = df1.loc[
        df1["cellcount"] > min_size, groupby
    ]  # smp_celltype of at least 20 cells
    vdj_list = vj_v_list + vj_j_list + vdj_v_list + vdj_j_list

    vdj_df = pd.DataFrame(columns=vdj_list, index=new_list)
    # this is a bit slow
    for i in tqdm(
        range(vdj_df.shape[0]),
        desc="Tabulating V/J gene usage",
        bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
        disable=not verbose,
    ):
        cell = vdj_df.index[i]
        counter1 = Counter(
            adata_.obs.loc[adata_.obs[groupby] == cell, v_call + mode + "_VJ"]
        )
        for vj_v in vj_v_list:
            vdj_df.loc[cell, vj_v] = counter1[vj_v]

        counter2 = Counter(
            adata_.obs.loc[
                adata_.obs[groupby] == cell, "j_call_" + mode + "_VJ"
            ]
        )
        for vj_j in vj_j_list:
            vdj_df.loc[cell, vj_j] = counter2[vj_j]

        counter3 = Counter(
            adata_.obs.loc[adata_.obs[groupby] == cell, v_call + mode + "_VDJ"]
        )
        for vdj_v in vdj_v_list:
            vdj_df.loc[cell, vdj_v] = counter3[vdj_v]
        counter5 = Counter(
            adata_.obs.loc[
                adata_.obs[groupby] == cell, "j_call_" + mode + "_VDJ"
            ]
        )
        for vdj_j in vdj_j_list:
            vdj_df.loc[cell, vdj_j] = counter5[vdj_j]
        # normalise
        vdj_df.loc[cell, vdj_df.columns.isin(vj_v_list)] = (
            vdj_df.loc[cell, vdj_df.columns.isin(vj_v_list)]
            / np.sum(vdj_df.loc[cell, vdj_df.columns.isin(vj_v_list)])
            * 100
        )
        vdj_df.loc[cell, vdj_df.columns.isin(vj_j_list)] = (
            vdj_df.loc[cell, vdj_df.columns.isin(vj_j_list)]
            / np.sum(vdj_df.loc[cell, vdj_df.columns.isin(vj_j_list)])
            * 100
        )
        vdj_df.loc[cell, vdj_df.columns.isin(vdj_v_list)] = (
            vdj_df.loc[cell, vdj_df.columns.isin(vdj_v_list)]
            / np.sum(vdj_df.loc[cell, vdj_df.columns.isin(vdj_v_list)])
            * 100
        )
        vdj_df.loc[cell, vdj_df.columns.isin(vdj_j_list)] = (
            vdj_df.loc[cell, vdj_df.columns.isin(vdj_j_list)]
            / np.sum(vdj_df.loc[cell, vdj_df.columns.isin(vdj_j_list)])
            * 100
        )

    df2 = pd.DataFrame(index=vdj_df.index, columns=["cell_type"])
    df2["cell_type"] = list(vdj_df.index)
    vdj_df_adata = sc.AnnData(
        X=vdj_df.values, obs=df2, var=pd.DataFrame(index=vdj_df.columns)
    )
    vdj_df_adata.obs["cell_count"] = pd.Series(
        dict(zip(df1[groupby], df1["cellcount"]))
    )
    sc.pp.pca(
        vdj_df_adata, n_comps=n_comps, use_highly_variable=False, **kwargs
    )
    if transfer_mapping is not None:
        collapsed_obs = adata_.obs.drop_duplicates(subset=groupby)
        for to in transfer_mapping:
            transfer_dict = dict(zip(collapsed_obs[groupby], collapsed_obs[to]))
            vdj_df_adata.obs[to] = [
                transfer_dict[x] for x in vdj_df_adata.obs_names
            ]
    logg.info(
        " finished",
        time=start,
        deep=("Returned AnnData: \n" "   'obsm', X_pca for V/J gene usage"),
    )
    return vdj_df_adata


def group_sequences(
    input_vdj: pd.DataFrame,
    junction_key: str,
    recalculate_length: bool = True,
    by_alleles: bool = False,
    locus: Literal["ig", "tr-ab", "tr-gd"] = "ig",
):
    """
    Groups sequence IDs with the same V and J genes, and then splits the groups based on the lengths of the sequences.

    Parameters
    ----------
    input_vdj : pd.DataFrame
        The input data frame.
    junction_key : str
        The name of the column in the input data that contains the junction sequences.
    recalculate_length : bool, optional
        Whether to recalculate the lengths of the junction sequences.
    by_alleles : bool, optional
        Whether to group sequence IDs by their V and J gene alleles, rather than by only the V and J genes.
    locus : Literal["ig", "tr-ab", "tr-gd"], optional
        The locus of the input data. One of "ig", "tr-ab", or "tr-gd".

    Returns
    -------
    vj_len_grp : Tree
        A nested dictionary that groups sequence IDs by V and J gene, and then by sequence length.
    seq_grp : Tree
        A nested dictionary that groups sequence IDs by V and J gene, and then by sequence.

    Raises
    ------
    ValueError
        Raised when the sequence length column is not found in the input table.
    """
    locus_log1_dict = {"ig": "IGH", "tr-ab": "TRB", "tr-gd": "TRD"}
    # retrieve the J genes and J genes
    if not by_alleles:
        if VCALLG in input_vdj.columns:
            V = [re.sub(STRIPALLELENUM, "", str(v)) for v in input_vdj[VCALLG]]
        else:
            V = [re.sub(STRIPALLELENUM, "", str(v)) for v in input_vdj[VCALL]]
        J = [re.sub(STRIPALLELENUM, "", str(j)) for j in input_vdj[JCALL]]
    else:
        if VCALLG in input_vdj.columns:
            V = [str(v) for v in input_vdj[VCALLG]]
        else:
            V = [str(v) for v in input_vdj[VCALL]]
        J = [str(j) for j in input_vdj[JCALL]]
    # collapse the alleles to just genes
    V = [",".join(list(set(v.split(",")))) for v in V]
    J = [",".join(list(set(j.split(",")))) for j in J]
    seq = dict(zip(input_vdj.index, input_vdj[junction_key]))
    if recalculate_length:
        seq_length = [len(str(l)) for l in input_vdj[junction_key]]
    else:
        try:
            seq_length = [l for l in input_vdj[junction_key + "_length"]]
        except:
            raise ValueError(
                "{} not found in {} input table.".format(
                    junction_key + "_length", locus_log1_dict[locus]
                )
            )
    seq_length_dict = dict(zip(input_vdj.index, seq_length))
    # Create a dictionary and group sequence ids with same V and J genes
    V_J = dict(zip(input_vdj.index, zip(V, J)))
    vj_grp = defaultdict(list)
    for key, val in sorted(V_J.items()):
        vj_grp[val].append(key)
    # and now we split the groups based on lengths of the seqs
    vj_len_grp = Tree()
    seq_grp = Tree()
    for g in vj_grp:
        # first, obtain what's the unique lengths
        jlen = []
        for contig_id in vj_grp[g]:
            jlen.append(seq_length_dict[contig_id])
        setjlen = list(set(jlen))
        # then for each unique length, we add the contigs to a new tree if it matches the length
        # and also the actual seq sequences into a another one
        for s in setjlen:
            for contig_id in vj_grp[g]:
                jlen_ = seq_length_dict[contig_id]
                if jlen_ == s:
                    vj_len_grp[g][s][contig_id].value = 1
                    seq_grp[g][s][seq[contig_id]].value = 1
                    for c in [contig_id]:
                        vj_len_grp[g][s][c] = seq[c]
    return vj_len_grp, seq_grp


def group_pairwise_hamming_distance(
    clonotype_vj_len_group: Tree,
    clonotype_sequence_group: Tree,
    identity: float,
    locus: str,
    chain: Literal["VDJ", "VJ"],
    verbose: bool = True,
) -> Tree:
    """
    Group clonotypes by pairwise hamming distance and generate a nested dictionary of clonotype groups.

    Parameters
    ----------
    clonotype_vj_len_group : Tree
        A nested dictionary that groups sequence IDs by V and J gene, and then by sequence length.
    clonotype_sequence_group : Tree
        A nested dictionary that groups sequence IDs by V and J gene, and then by sequence.
    identity : float
        The identity threshold for grouping sequences.
    locus : str
        The locus of the chain.
    chain : Literal["VDJ", "VJ"]
        The chain type.
    verbose : bool, optional
        Whether to show the progress bar, by default True.

    Returns
    -------
    Tree
        A nested dictionary containing contigs assigned to clonotype groups.
    """
    clones = Tree()

    # for each seq group, calculate the hamming distance matrix
    for g in tqdm(
        clonotype_sequence_group,
        desc=f"Finding clones based on {locus} cell {chain} chains ".format(),
        bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
        disable=not verbose,
    ):
        for l in clonotype_sequence_group[g]:
            seq_ = list(clonotype_sequence_group[g][l])
            tdarray = np.array(seq_).reshape(-1, 1)
            d_mat = squareform(pdist(tdarray, lambda x, y: hamming(x[0], y[0])))
            # then calculate what the acceptable threshold is for each length of sequence
            tr = math.floor(int(l) * (1 - identity))
            # convert diagonal and upper triangle to zeroes
            d_mat = np.tril(d_mat)
            np.fill_diagonal(d_mat, 0)
            # get the coordinates/indices of seqs to match against the threshold later
            indices_temp = []
            indices = []
            indices_temp = [list(x) for x in np.tril_indices_from(d_mat)]
            indices = list(zip(indices_temp[0], indices_temp[1]))
            # if there's more than 1 contig, remove the diagonal
            if len(indices) > 1:
                for pairs in indices:
                    # remove diagonals
                    if pairs[0] == pairs[1]:
                        indices.remove(pairs)
            indices_j = []
            # use the coordinates/indices to retrieve the seq sequences
            for p in range(0, len(indices)):
                a1, b1 = indices[p]
                indices_j.append(seq_[a1])
                indices_j.append(seq_[b1])
            # convert the distance matrix to coordinate (source) and distance (target) and create it
            # as a dictionary
            source, target = d_mat.nonzero()
            source_target = list(zip(source.tolist(), target.tolist()))
            if len(source) == 0 & len(target) == 0:
                source_target = list([(0, 0)])
            dist = {}
            for st in source_target:
                dist.update({st: d_mat[st]})
            if d_mat.shape[0] > 1:
                seq_tmp_dict = clustering(dist, tr, seq_)
            else:
                seq_tmp_dict = {seq_[0]: tuple([seq_[0]])}
            # sort the list so that clones that are larger have a smaller number
            clones_tmp = sorted(
                list(set(seq_tmp_dict.values())),
                key=len,
                reverse=True,
            )
            for x in range(0, len(clones_tmp)):
                clones[g][l][x + 1] = clones_tmp[x]
    # now to retrieve the contig ids that are grouped together
    cid = Tree()
    for g in clones:
        for l in clones[g]:
            # retrieve the clone 'numbers'
            for c in clones[g][l]:
                grp_seq = clones[g][l][c]
                for key, value in clonotype_vj_len_group[g][l].items():
                    if value in grp_seq:
                        cid[g][l][c][key].value = 1

    return cid


def rename_clonotype_ids(
    clonotype_groups: Tree,
    cell_type: str = "",
    suffix: str = "",
    prefix: str = "",
) -> dict[str, str]:
    """
    Renames clonotype IDs to numerical barcode system.

    Parameters
    ----------
    clonotype_groups : Tree
        A nested dictionary that containing clonotype groups of contigs.
    cell_type : str, optional
        Cell type name to append to front, if necessary.
    suffix : str, optional
        Suffix to append to the end of the clonotype ID.
    prefix : str, optional
        Prefix to append to the beginning of the clonotype ID.
    Returns
    -------
    dict[str, str]
        A dictionary that maps sequence IDs to clonotype IDs.
    """
    clone_dict = {}
    first_key = []
    for k1 in clonotype_groups.keys():
        first_key.append(k1)
    first_key = list(set(first_key))
    first_key_dict = dict(zip(first_key, range(1, len(first_key) + 1)))
    for g in clonotype_groups:
        second_key = []
        for k2 in clonotype_groups[g].keys():
            second_key.append(k2)
        second_key = list(set(second_key))
        second_key_dict = dict(zip(second_key, range(1, len(second_key) + 1)))
        for l in clonotype_groups[g]:
            third_key = []
            for k3 in clonotype_groups[g][l].keys():
                third_key.append(k3)
            third_key = list(set(third_key))
            third_key_dict = dict(zip(third_key, range(1, len(third_key) + 1)))
            for key, value in dict(clonotype_groups[g][l]).items():
                for v in value:
                    if type(v) is int:
                        break
                    clone_dict[v] = (
                        cell_type
                        + prefix
                        + str(first_key_dict[g])
                        + "_"
                        + str(second_key_dict[l])
                        + "_"
                        + str(third_key_dict[key])
                        + suffix
                    )

    return clone_dict


def refine_clone_assignment(
    dat: pd.DataFrame,
    clone_key: str,
    clone_dict_vj: dict,
    verbose: bool = True,
) -> None:
    """
    Refines the clone assignment of VDJ sequences based on their VJ chain pairing.

    Parameters
    ----------
    dat : pd.DataFrame
        The input data frame containing the full VDJ/VJ data.
    clone_key : str
        The name of the column in the input data that contains the clone IDs.
    clone_dict_vj : dict
        A dictionary for the VJ chains for us to map sequence IDs to clone IDs.
    verbose : bool, optional
        Whether to print progress messages.
    """
    cellclonetree = Tree()
    seqcellclonetree = Tree()
    for c, s, z in zip(dat["cell_id"], dat["sequence_id"], dat[clone_key]):
        seqcellclonetree[c][s].value = 1
        if pd.notnull(z):
            cellclonetree[c][z].value = 1

    for c in cellclonetree:
        cellclonetree[c] = list(cellclonetree[c])

    fintree = Tree()
    for c in tqdm(
        cellclonetree,
        desc="Refining clone assignment based on VJ chain pairing ",
        bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
        disable=not verbose,
    ):
        suffix = [
            clone_dict_vj[x] for x in seqcellclonetree[c] if x in clone_dict_vj
        ]
        fintree[c] = []
        if len(suffix) > 0:
            for cl in cellclonetree[c]:
                if present(cl):
                    for s in suffix:
                        if check_same_celltype(cl, s):
                            fintree[c].append(
                                cl + "_" + "".join(s.split("_", 1)[1])
                            )
                        else:
                            fintree[c].append(cl + "|" + s)
        else:
            for cl in cellclonetree[c]:
                if present(cl):
                    fintree[c].append(cl)
        fintree[c] = "|".join(fintree[c])
    dat[clone_key] = [fintree[x] for x in dat["cell_id"]]
    for i, row in dat.iterrows():  # is this going to be slow...?
        if not present(row[clone_key]):
            if i in clone_dict_vj and i in dat.index:
                dat.at[i, clone_key] = clone_dict_vj[i]


def check_chains(dat_vdj: pd.DataFrame, dat_vj: pd.DataFrame) -> pd.DataFrame:
    """
    Generate a summary table of whether the chain is orphan or not.

    Parameters
    ----------
    dat_vdj : pd.DataFrame
        Input dataframe containing VDJ chains of a particular locus.
    dat_vj : pd.DataFrame
        Input dataframe containing VJ chains of a particular locus.

    Returns
    -------
    pd.DataFrame
        Output dataframe containing chain status.
    """

    vj_check = pd.crosstab(dat_vj.cell_id, dat_vj.locus).apply(sum, axis=1)
    vj_check[vj_check > 1] = 1
    vdj_check = pd.crosstab(dat_vdj.cell_id, dat_vdj.locus).apply(sum, axis=1)
    vdj_check[vdj_check > 1] = 1
    chain_check = pd.concat([vdj_check, vj_check], axis=1)
    chain_check.columns = ["VDJ", "VJ"]
    chain_check["Orphan VDJ"] = pd.notnull(chain_check["VDJ"]) & (
        pd.isnull(chain_check["VJ"])
    )
    chain_check["Orphan VJ"] = pd.notnull(chain_check["VJ"]) & (
        pd.isnull(chain_check["VDJ"])
    )
    chain_check["All VDJ"] = pd.notnull(chain_check["VDJ"])
    chain_check["All VJ"] = pd.notnull(chain_check["VJ"])
    chain_check["VDJ"] = (
        (~chain_check["Orphan VDJ"])
        & (~chain_check["Orphan VJ"])
        & (pd.notnull(chain_check["VDJ"]))
    )
    chain_check["VJ"] = (
        (~chain_check["Orphan VDJ"])
        & (~chain_check["Orphan VJ"])
        & (pd.notnull(chain_check["VJ"]))
    )
    return chain_check
