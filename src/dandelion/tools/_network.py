from __future__ import annotations

import multiprocessing
import re
import time

import networkx as nx
import numpy as np
import pandas as pd

from collections import defaultdict
from joblib import Parallel, delayed
from pathlib import Path
from polyleven import levenshtein
from scanpy import logging as logg
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from tqdm import tqdm
from typing import Callable, Literal, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


from dandelion.tools._tools import vdj_sample
from dandelion.tools._layout import generate_layout
from dandelion.utilities._core import Dandelion, Query
from dandelion.utilities._distances import (
    Metric,
    dist_func_long_sep,
    resolve_metric,
)
from dandelion.utilities._utilities import (
    present,
    sanitize_data,
    Tree,
    FALSES,
)


def generate_network(
    vdj_data: Dandelion,
    adata: AnnData | None = None,
    key: str | None = None,
    clone_key: str | None = None,
    min_size: int = 2,
    sample: int | None = None,
    force_replace: bool = False,
    verbose: bool = True,
    compute_graph: bool = True,
    compute_layout: bool = True,
    layout_method: Literal["mod_fr", "sfdp"] = "sfdp",
    expanded_only: bool = False,
    use_existing_graph: bool = True,
    num_cores: int = 1,
    sequential_chain: bool = False,
    distance_mode: Literal["clone", "full"] = "clone",
    dist_func: Callable | str | None = None,
    pad_to_max: bool = False,
    lazy: bool = False,
    zarr_path: Path | str | None = None,
    chunk_size: int | None = None,
    chunk_clone_limit: int | None = None,
    memory_limit_gb: float | None = None,
    memory_safety_fraction: float = 0.3,
    compress: bool = True,
    random_state: int | np.random.RandomState | None = None,
    **kwargs,
) -> Dandelion | tuple[Dandelion, AnnData]:
    """
    Generate a Levenshtein distance network based on VDJ and VJ sequences.

    The distance matrices are then combined into a singular matrix.

    Parameters
    ----------
    vdj_data : Dandelion
        Dandelion object.
    key : str | None, optional
        column name for distance calculations. None defaults to 'sequence_alignment_aa'.
    clone_key : str | None, optional
        column name to build network on.
    min_size : int, optional
        For visualization purposes, two graphs are created where one contains all cells and a trimmed second graph.
        This value specifies the minimum number of edges required otherwise node will be trimmed in the secondary graph.
    sample : int | None, optional
        If specified, cells will be randomly sampled to the integer provided. If the integer is larger than the number of cells,
        sampling with replacement is used and the same cell may appear multiple times with different sequence and cell ids. If None,
        no resampling is performed. A new Dandelion class will be returned.
    force_replace : bool, optional
        whether or not to sample with replacement when `sample` is smaller or equal to than the number of cells.
    verbose : bool, optional
        whether or not to print the progress bars.
    compute_graph : bool, optional
        whether or not to generate the graph after distance matrix calculation.
    compute_layout : bool, optional
        whether or not to generate the layout. May be time consuming if too many cells.
    layout_method : Literal["sfdp", "mod_fr"], optional
        accepts one of 'sfdp' or 'mod_fr'. 'sfdp' refers to `sfdp_layout` from `graph_tool` (C++ implementation; fast)
        whereas 'mod_fr' refers to modified Fruchterman-Reingold layout originally implemented in dandelion (python
        implementation; slow).
    expanded_only : bool, optional
        whether or not to only compute layout on expanded clonotypes.
    use_existing_graph : bool, optional
        whether or not to just compute the layout using the existing graph if it exists in the object.
    num_cores : int, optional
        number of cores to use for parallelizable steps. -1 uses all available cores.
    sequential_chain : bool, optional
        whether or not to use the original method for distance calculation method where each chain is calculated
        separately and sequentially added to the total distance matrix. This method is slower but would be more
        precise calculation. If False, concatenated sequences with a long separator are used for distance calculation.
        Ignored if lazy=True as the lazy method always uses the long separator approach. The long separator approach
        inserts a long string of consistent characters on a per-chain basis to ensure that distances between chains are large
        and do not interfere with intra-chain distances.
    distance_mode : Literal["clone", "full"], optional
        method to compute distance matrix. 'clone' refers to the original membership-based distance calculation where
        only distances within clones are calculated. Whereas 'full' computes the full pairwise distance matrix.
    dist_func : Callable | str | None, optional
        distance function to use. If None, `polyleven.levenshtein` is used. If a string is provided, it will use Bio.Align's
        substitution matrices (e.g., 'BLOSUM62', 'PAM250'). See `Bio.Align.substitution_matrices.load` for available options.
    blosom_matrix : str | None, optional
        If provided, dist_func is ignored and this substitution matrix from Bio.Align is used instead.
    pad_to_max : bool, optional
        whether or not to pad sequences to the maximum length in the dataset before distance calculation. This will
        allow for distance calculations that need sequences of the same length (e.g., Hamming distance). Note that this
        may increase memory usage and computation time.
    lazy: bool, optional
        If True, computation will be performed lazily using Dask/Zarr arrays. True will also return a Dask array view of the
        distance matrix stored on disk instead of a numpy array stored in memory.
    zarr_path: Path | str | None, optional
        Path to store Zarr array when using lazy mode. If None, "distance_matrix.zarr" will be created in the current working directory.
    chunk_size: int | None, optional
        Chunk size for distance matrix computation when using lazy mode. If None, chunk size is automatically computed
        based on available memory and number of cores. The automatic chunk size can be further adjusted using
        `memory_limit_gb` and `memory_safety_fraction` parameters.
    chunk_clone_limit: int | None, optional
        Maximum number of clones to process per chunk when using lazy mode and distance method = "clone". If None, chunk sizes will be
        automatically determined based on available memory and number of cores.
    memory_limit_gb: float | None, optional
        Memory limit per worker in GB for Dask. None defaults to all available memory/cores.
    memory_safety_fraction: float, optional
        Fraction of available memory to use. Defaults to 0.3 (i.e., 30% of available memory will be used for chunk size calculation).
    compress: bool, optional
        Whether to compress the Zarr array using Blosc with zstd.
    **kwargs
        additional kwargs passed to options specified in `networkx.drawing.layout.spring_layout` or
        `graph_tool.draw.sfdp_layout`.

    Returns
    -------
    Dandelion | tuple[Dandelion, AnnData]
        Dandelionbject with `.edges`, `.layout`, `.graph` initialized.

    Raises
    ------
    ValueError
        if any errors with dandelion input.
    """
    # normalize num_cores convention (-1 => use all CPUs)
    if num_cores == -1:
        num_cores = multiprocessing.cpu_count()
    num_cores = max(1, int(num_cores))

    clone_key = clone_key if clone_key is not None else "clone_id"
    dist_func = levenshtein if dist_func is None else dist_func
    metric = resolve_metric(dist_func)
    if not compute_graph:
        compute_layout = False

    if distance_mode == "clone" or compute_graph or compute_layout:
        if clone_key not in vdj_data.data:
            raise ValueError(
                "Data does not contain clone information. Please run ddl.tl.find_clones."
            )
    regenerate = True
    if vdj_data.graph is not None:
        if (min_size != 2) or (sample is not None):
            pass
        elif use_existing_graph:
            start = logg.info(
                "Generating network layout from pre-computed network"
            )
            if isinstance(vdj_data, Dandelion):
                regenerate = False
                g, g_, lyt, lyt_ = generate_layout(
                    vertices=None,
                    edges=None,
                    min_size=min_size,
                    weight=None,
                    verbose=verbose,
                    compute_layout=compute_layout,
                    layout_method=layout_method,
                    expanded_only=expanded_only,
                    graphs=(vdj_data.graph[0], vdj_data.graph[1]),
                    **kwargs,
                )

    if regenerate:
        start = logg.info("Generating network")

        key_ = key if key is not None else "sequence_alignment_aa"

        if key_ not in vdj_data.data:
            raise ValueError(f"key {key_} not found in data.")

        if lazy:
            from dandelion.tools._lazydistances import dask_safe_slice_square

        if sample is not None:
            vdj_data, adata = vdj_sample(
                vdj_data,
                adata,
                size=sample,
                force_replace=force_replace,
                random_state=random_state,
            )
            dat = vdj_data.data.copy()
        else:
            if "ambiguous" in vdj_data.data:
                dat = vdj_data.data[
                    vdj_data.data["ambiguous"].isin(FALSES)
                ].copy()
            else:
                dat = vdj_data.data.copy()

        dat_h = dat[dat["locus"].isin(["IGH", "TRB", "TRD"])].copy()
        dat_l = dat[dat["locus"].isin(["IGK", "IGL", "TRA", "TRG"])].copy()
        dat_ = pd.concat([dat_h, dat_l], ignore_index=True)
        dat_ = sanitize_data(dat_, ignore=clone_key)

        # retrieve sequence columns and clone info (unchanged)
        querier = Query(dat_, verbose=verbose)
        dat_seq = querier.retrieve(query=key_, retrieve_mode="split")
        # ensure that dat_seq matches order of vdj_data.metadata
        dat_seq = dat_seq.reindex(vdj_data.metadata.index)
        dat_seq.columns = [re.sub(key_ + "_", "", i) for i in dat_seq.columns]
        if distance_mode == "clone":
            dat_clone = querier.retrieve(
                query=clone_key, retrieve_mode="merge and unique only"
            )

            dat_clone = dat_clone[clone_key].str.split("|", expand=True)
            membership = Tree()
            for i, j in dat_clone.iterrows():
                jjj = [jj for jj in j if present(jj)]
                for ij in jjj:
                    membership[ij][i].value = 1
            membership = {i: list(j) for i, j in dict(membership).items()}

        # compute total_dist using chosen mode (original uses membership)
        logg.info(
            f"Calculating distance matrix {'lazily' if lazy else ''} with distance_mode = '{distance_mode}'\n"
        )
        if distance_mode == "clone":
            if lazy:
                from dandelion.tools._lazydistances import (
                    calculate_distance_matrix_zarr,
                )

                total_dist = calculate_distance_matrix_zarr(
                    dat_seq,
                    metric=metric,
                    pad_to_max=pad_to_max,
                    membership=membership,
                    zarr_path=zarr_path,
                    chunk_size=chunk_size,
                    max_clones_per_chunk=chunk_clone_limit,
                    num_cores=num_cores,
                    memory_limit_gb=memory_limit_gb,
                    memory_safety_fraction=memory_safety_fraction,
                    compress=compress,
                    lazy=lazy,
                    verbose=verbose,
                )
            else:
                if sequential_chain:
                    total_dist = calculate_distance_matrix_original(
                        dat_seq,
                        membership,
                        metric=metric,
                        verbose=verbose,
                    )
                else:
                    total_dist = calculate_distance_matrix_long(
                        dat_seq,
                        membership=membership,
                        metric=metric,
                        pad_to_max=pad_to_max,
                        num_cores=num_cores,
                        verbose=verbose,
                    )
        elif distance_mode == "full":
            if lazy:
                from dandelion.tools._lazydistances import (
                    calculate_distance_matrix_zarr,
                )

                total_dist = calculate_distance_matrix_zarr(
                    dat_seq,
                    metric=metric,
                    pad_to_max=pad_to_max,
                    membership=None,
                    zarr_path=zarr_path,
                    chunk_size=chunk_size,
                    num_cores=num_cores,
                    memory_limit_gb=memory_limit_gb,
                    memory_safety_fraction=memory_safety_fraction,
                    compress=compress,
                    lazy=lazy,
                    verbose=verbose,
                )
            else:
                if sequential_chain:
                    total_dist = calculate_distance_matrix_original_full(
                        dat_seq,
                        metric=metric,
                        num_cores=num_cores,
                        verbose=verbose,
                    )
                else:
                    total_dist = calculate_distance_matrix_long(
                        dat_seq,
                        membership=None,
                        metric=metric,
                        pad_to_max=pad_to_max,
                        num_cores=num_cores,
                        verbose=verbose,
                    )
        if compute_graph:

            # build cluster_dist as before (only include groups with >1 member)
            tmp_clusterdist2 = {
                c: membership[c] for c in membership if len(membership[c]) > 1
            }

            # extract per-clone distance matrices
            cluster_dist = {}
            for c_ in tqdm(
                tmp_clusterdist2,
                desc="Sorting into clusters ",
                disable=not verbose,
                bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
            ):
                idxs = [
                    i
                    for i in tmp_clusterdist2[c_]
                    if i in vdj_data.metadata.index
                ]
                if len(idxs) > 1:
                    pos = [dat_seq.index.get_loc(i) for i in idxs]
                    if lazy:
                        dist_mat_ = pd.DataFrame(
                            dask_safe_slice_square(total_dist, pos).compute(),
                            index=idxs,
                            columns=idxs,
                        )
                    else:
                        dist_mat_ = pd.DataFrame(
                            total_dist[np.ix_(pos, pos)],
                            index=idxs,
                            columns=idxs,
                        )
                    cluster_dist[c_] = dist_mat_
            # Minimum spanning trees (unchanged)
            # mst_tree = mst(cluster_dist, num_cores=num_cores, verbose=verbose)
            mst_tree = mst(cluster_dist, num_cores=1, verbose=verbose)
            # generate edge list
            edge_list = Tree()
            for c in tqdm(
                mst_tree,
                desc="Generating edge list ",
                disable=not verbose,
                bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
            ):
                edge_list[c] = nx.to_pandas_edgelist(mst_tree[c])
                if edge_list[c].shape[0] > 0:
                    edge_list[c]["weight"] = edge_list[c]["weight"] - 1
                    edge_list[c].loc[edge_list[c]["weight"] < 0, "weight"] = 0
            # try to combine edge lists (simple version)
            try:
                edge_list_final = pd.concat([edge_list[x] for x in edge_list])
                set_edge_list_index(edge_list_final)
            except Exception:
                edge_list_final = None

            # final layout + graph creation (unchanged)
            g, g_, lyt, lyt_ = generate_layout(
                vertices=vdj_data.metadata.index.tolist(),
                edges=edge_list_final,
                min_size=min_size,
                weight=None,
                verbose=verbose,
                compute_layout=compute_layout,
                layout_method=layout_method,
                expanded_only=expanded_only,
                **kwargs,
            )

    logg.info(
        " finished. Updated Dandelion object: \n",
        time=start,
        deep=(
            "   'layout', graph layout\n"
            if compute_layout
            else (
                ""
                "   'graph', network constructed from distance matrices of VDJ- and VJ- chains\n"
                if compute_graph
                else (
                    "" "   'distances', VDJ + VJ distance matrix\n"
                    if regenerate
                    else ""
                )
            )
        ),
    )

    # return or re-initialize vdj_data
    germline = getattr(vdj_data, "germline", None)
    if regenerate:
        distances = total_dist if lazy else csr_matrix(total_dist)
    else:
        distances = vdj_data.distances
    graph = (g, g_) if compute_graph else None
    layout = (lyt, lyt_) if compute_graph and compute_layout else None
    if sample is not None:
        out = Dandelion(
            data=dat_,
            metadata=vdj_data.metadata,
            clone_key=clone_key,
            layout=layout,
            graph=graph,
            distances=distances,
            germline=germline,
            verbose=False,
        )
        if adata is None:
            return out
        else:
            return out, adata
    else:
        vdj_data.__init__(
            data=vdj_data.data,
            metadata=vdj_data.metadata,
            clone_key=clone_key,
            layout=layout,
            graph=graph,
            distances=distances,
            germline=germline,
            initialize=False,
            verbose=False,
        )


def mst(
    mat: dict,
    num_cores: int | None = None,
    verbose: bool = True,
) -> Tree:
    """
    Construct minimum spanning tree based on supplied matrix in dictionary.

    Parameters
    ----------
    mat : dict
        Dictionary containing numpy ndarrays.
    num_cores: int, optional
        Number of cores to run this step. Parallelise using joblib if more than 1.
    verbose : bool, optional
        Whether or not to show logging information.

    Returns
    -------
    Tree
        Dandelion `Tree` object holding DataFrames of constructed minimum spanning trees.
    """
    mst_tree = Tree()

    if num_cores == 1:
        for c in tqdm(
            mat,
            desc="Calculating minimum spanning tree ",
            disable=not verbose,
            bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
        ):
            _, mst_tree[c] = process_mst_per_clonotype(mat=mat, c=c)
    else:
        results = Parallel(n_jobs=num_cores)(
            delayed(process_mst_per_clonotype)(mat, c)
            for c in tqdm(
                mat,
                desc=f"Calculating minimum spanning tree, parallelized across {num_cores} cores ",
                disable=not verbose,
                bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
            )
        )
        for result in results:
            mst_tree[result[0]] = result[1]
    return mst_tree


def process_mst_per_clonotype(
    mat: dict[str, pd.DataFrame], c: str
) -> tuple[str, nx.Graph]:
    """
    Function to calculate minimum spanning tree.

    Parameters
    ----------
    mat : dict[str, pd.DataFrame]
        Dictionary holding distance matrices.
    c : str
        Name of clonotype.

    Returns
    -------
    tuple[str, nx.Graph]
        Graph holding minimum spanning tree.
    """
    tmp = mat[c] + 1
    tmp[np.isnan(tmp)] = 0
    G = create_networkx_graph(tmp, drop_zero=True)
    return c, nx.minimum_spanning_tree(G)


def create_networkx_graph(
    adjacency: pd.DataFrame, drop_zero: bool = True
) -> nx.Graph:
    """
    Create a networkx graph from an adjacency matrix in chunks.

    Parameters
    ----------
    adjacency : pd.DataFrame
        Adjacency matrix.
    drop_zero : bool, optional
        Whether or not to drop edges with zero weight, by default True.

    Returns
    -------
    nx.Graph
        NetworkX graph.
    """
    G = nx.Graph()
    # add nodes
    nodes = list(adjacency.index)
    G.add_nodes_from(nodes)
    # convert adjacency matrix to edge list
    edge_list = adjacency_to_edge_list(adjacency, drop_zero=drop_zero)
    G.add_weighted_edges_from(ebunch_to_add=edge_list.values)
    return G


def set_edge_list_index(edge_list: pd.DataFrame) -> None:
    """
    Set the index of the edge list in-place.

    Parameters
    ----------
    edge_list : pd.DataFrame
        Edge list.
    """
    edge_list.index = [
        str(s) + "|" + str(t)
        for s, t in zip(edge_list["source"], edge_list["target"])
    ]


def adjacency_to_edge_list(
    adjacency: pd.DataFrame, drop_zero: bool = False, rename_index: bool = False
) -> pd.DataFrame:
    """
    Convert adjacency matrix to edge list that excludes self-loops.

    Parameters
    ----------
    adjacency : pd.DataFrame
        Adjacency matrix.
    drop_zero : bool, optional
        Whether or not to drop edges with zero weight, by default False.
    rename_index : bool, optional
        Whether or not to rename the index, by default False.

    Returns
    -------
    pd.DataFrame
        Edge list.
    """
    edge_list = adjacency.stack().reset_index()
    edge_list.columns = ["source", "target", "weight"]
    if rename_index:
        set_edge_list_index(edge_list)
    if drop_zero:
        edge_list = edge_list[edge_list["weight"] != 0]
    return edge_list


def clone_degree(vdj_data: Dandelion, weight: str | None = None) -> Dandelion:
    """
    Calculate node degree in BCR/TCR network.

    Parameters
    ----------
    vdj_data : Dandelion
        Dandelionect after `tl.generate_network` has been run.
    weight : str | None, optional
        Attribute name for retrieving edge weight in graph. None defaults to ignoring this. See `networkx.Graph.degree`.

    Raises
    ------
    AttributeError
        if graph not found.
    TypeError
        if input is not Dandelion class.
    """
    if isinstance(vdj_data, Dandelion):
        if vdj_data.graph is None:
            raise AttributeError(
                "Graph not found. Please run tl.generate_network."
            )
        else:
            G = vdj_data.graph[0]
            cd = pd.DataFrame.from_dict(G.degree(weight=weight))
            cd.set_index(0, inplace=True)
            vdj_data.metadata["clone_degree"] = pd.Series(cd[1])
    else:
        raise TypeError("Input object must be of {}".format(Dandelion))


def clone_centrality(vdj_data: Dandelion):
    """
    Calculate node closeness centrality in BCR/TCR network.

    Parameters
    ----------
    vdj_data : Dandelion
        Dandelion object after `tl.generate_network` has been run.

    Raises
    ------
    AttributeError
        if graph not found.
    TypeError
        if input is not Dandelion class.
    """
    if isinstance(vdj_data, Dandelion):
        if vdj_data.graph is None:
            raise AttributeError(
                "Graph not found. Please run tl.generate_network."
            )
        else:
            G = vdj_data.graph[0]
            cc = nx.closeness_centrality(G)
            cc = pd.DataFrame.from_dict(
                cc, orient="index", columns=["clone_centrality"]
            )
            vdj_data.metadata["clone_centrality"] = pd.Series(
                cc["clone_centrality"]
            )
    else:
        raise TypeError("Input object must be of {}".format(Dandelion))


def calculate_distance_matrix_original(
    dat_seq: pd.DataFrame,
    membership: dict,
    metric: Metric,
    verbose: bool = True,
) -> np.ndarray:
    """
    Re-implementation of original membership-based distance calculation.

    Parameters
    ----------
    dat_seq : pd.DataFrame
        DataFrame with sequence columns; index corresponds to cell IDs (or whatever unique ids you use).
    membership : dict
        Mapping from clone_id -> list of indices (these indices must be present in dat_seq.index).
    metric : Metric
        Distance metric to use.
    verbose : bool, optional
        Whether to show progress.

    Returns
    -------
    total_dist : np.ndarray
        Aggregated distance matrix across all columns; diagonal set to NaN by caller.
    """
    n = dat_seq.shape[0]
    index_list = list(dat_seq.index)
    # dmat_per_column will collect for each column a list of DataFrames (one per clone) to concat
    dmat_per_column = defaultdict(list)
    # iterate clones (membership) exactly like original
    for clone in tqdm(
        membership,
        desc="Calculating distances (original, by membership)",
        disable=not verbose,
        bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
    ):
        tmp = dat_seq.loc[membership[clone]]
        if tmp.shape[0] > 1:
            tmp = (
                tmp.replace("[.]", "", regex=True)
                .fillna("")
                .replace("None", "")
            )
            for col in tmp.columns:
                tdarray = np.array(tmp[col]).reshape(-1, 1)
                # keep the original NaN-check logic: return 0 when either is NaN
                d_mat_tmp = squareform(
                    pdist(
                        tdarray,
                        lambda x, y: (
                            metric.compute(x[0], y[0])
                            if (pd.notnull(x[0]) and pd.notnull(y[0]))
                            else 0
                        ),
                    )
                )
                df_block = pd.DataFrame(
                    d_mat_tmp, index=tmp.index, columns=tmp.index
                )
                dmat_per_column[col].append(df_block)
    # For each column, concat its blocks, resolve duplicates (sum), reindex to full index
    dist_matrices = []
    for col, blocks in dmat_per_column.items():
        if not blocks:
            continue
        full = pd.concat(blocks)
        # If duplicates occur at the index-level, group & sum them (same as original)
        if any(full.index.duplicated()):
            dup_indices = full.index[full.index.duplicated()]
            tmp1 = full.drop(dup_indices)
            tmp2 = full.loc[dup_indices]
            tmp2 = tmp2.groupby(level=0).apply(lambda df: df.sum(axis=0))
            full = pd.concat([tmp1, tmp2])
        # reindex to full matrix indices (missing rows/cols -> fill with zeros)
        full = full.reindex(index=index_list, columns=index_list).fillna(0.0)
        dist_matrices.append(full.values)
    if len(dist_matrices) == 0:
        total_dist = np.zeros((n, n))
    else:
        total_dist = np.sum(dist_matrices, axis=0)

    np.fill_diagonal(total_dist, np.nan)
    return total_dist


def calculate_distance_matrix_original_full(
    dat_seq: pd.DataFrame,
    metric: Metric,
    num_cores: int = 1,
    verbose: bool = True,
) -> np.ndarray:
    """
    Re-implementation of original membership-based distance calculation.

    Parameters
    ----------
    dat_seq : pd.DataFrame
        DataFrame with sequence columns; index corresponds to cell IDs (or whatever unique ids you use).
    metric : Metric
        Distance metric to use.
    num_cores : int, optional
        Number of cores to run this step. Parallelise using joblib if more than 1
    verbose : bool, optional
        Whether to show progress.

    Returns
    -------
    total_dist : np.ndarray
        Aggregated distance matrix across all columns; diagonal set to NaN by caller.
    """
    start_time = time.time()
    n = dat_seq.shape[0]
    total_dist = np.zeros((n, n), dtype=float)

    for col in dat_seq.columns:
        seq_series = (
            dat_seq[col]
            .replace("[.]", "", regex=True)
            .fillna("")
            .replace("None", "")
        )
        nonnull = seq_series.dropna()
        if nonnull.shape[0] <= 1:
            continue
        # seq_indices = list(nonnull.index)
        seqs = nonnull.to_numpy(dtype=object)
        m = len(seqs)
        if num_cores > 1:
            results = Parallel(n_jobs=num_cores)(
                delayed(_compute_row)(i, seqs, m, metric)
                for i in tqdm(
                    range(n), disable=not verbose, leave=False, position=0
                )
            )
            for i in tqdm(
                range(n), disable=not verbose, leave=False, position=0
            ):
                for j in range(i + 1, n):
                    total_dist[i][j] += results[i][j]
                    total_dist[j][i] += results[i][j]
        else:
            results = squareform(
                pdist(
                    seqs.reshape(-1, 1),
                    lambda x, y: (
                        metric.compute(x[0], y[0])
                        if (pd.notnull(x[0]) and pd.notnull(y[0]))
                        else 0
                    ),
                )
            )
            total_dist += results

    np.fill_diagonal(total_dist, np.nan)
    if verbose:
        end_time = time.time()
        logg.info(
            f"Distances calculated in {end_time - start_time:.2f} seconds"
        )
    return total_dist


def calculate_distance_matrix_long(
    dat_seq: pd.DataFrame,
    membership: dict | None,
    metric: Metric,
    pad_to_max: bool = False,
    sep: str = "#",
    num_cores: int = 1,
    verbose: bool = True,
) -> np.ndarray:
    """
    Re-implementation of original membership-based distance calculation but using concatenated sequences
    using a long separator.

    Parameters
    ----------
    dat_seq : pd.DataFrame
        DataFrame with sequence columns; index corresponds to cell IDs (or whatever unique ids you use).
    membership : dict | None
        Mapping from clone_id -> list of indices (these indices must be present in dat_seq.index).
        None indicates full pairwise distance calculation.
    metric : Metric
        Distance metric to use.
    pad_to_max : bool, optional
        whether or not to pad sequences to the maximum length in the dataset before distance calculation. This will
        allow for distance calculations that need sequences of the same length (e.g., Hamming distance). Note that this
        may increase memory usage and computation time.
    num_cores : int, optional
        Number of cores to run this step. Parallelise using joblib if more than 1.
    verbose : bool, optional
        Whether to show progress.

    Returns
    -------
    total_dist : np.ndarray (n x n)
        Aggregated distance matrix across all columns; diagonal set to NaN by caller.
    """
    start_time = time.time()
    # Step 1: clean sequences
    dat_seq_clean = (
        dat_seq.replace("[.]", "", regex=True).fillna("").replace("None", "")
    )
    # Step 2: initialize distance matrix
    n = dat_seq_clean.shape[0]
    index_list = list(dat_seq_clean.index)
    total_dist = np.zeros((n, n))
    if membership is None:
        # Step 3: compute full distance matrix at once
        seqs = dat_seq_clean.values
        if num_cores > 1:
            results = Parallel(n_jobs=num_cores)(
                delayed(_compute_row)(i, seqs, n, metric)
                for i in tqdm(
                    range(n), disable=not verbose, leave=False, position=0
                )
            )
            for i in tqdm(
                range(n), disable=not verbose, leave=False, position=0
            ):
                for j in range(i + 1, n):
                    total_dist[i][j] += results[i][j]
                    total_dist[j][i] += results[i][j]
        else:
            results = squareform(
                pdist(
                    seqs,
                    metric=lambda x, y: dist_func_long_sep(
                        x,
                        y,
                        metric=metric,
                        pad_to_max=pad_to_max,
                        sep=sep,
                    ),
                )
            )
            total_dist += results
    else:
        # Step 3: iterate over clone memberships
        for clone in tqdm(
            membership,
            desc="Calculating distances (original, by membership)",
            disable=not verbose,
            bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
        ):
            tmp = dat_seq_clean.loc[membership[clone]]
            if tmp.shape[0] > 1:
                seqs = tmp.values
                d_mat_tmp = squareform(
                    pdist(
                        seqs,
                        metric=lambda x, y: dist_func_long_sep(
                            x,
                            y,
                            metric=metric,
                            pad_to_max=pad_to_max,
                            sep=sep,
                        ),
                    )
                )
                total_dist[
                    np.ix_(
                        [index_list.index(i) for i in tmp.index],
                        [index_list.index(j) for j in tmp.index],
                    )
                ] += d_mat_tmp

    np.fill_diagonal(total_dist, np.nan)
    if verbose:
        end_time = time.time()
        logg.info(
            f"Distances calculated in {end_time - start_time:.2f} seconds"
        )
    return total_dist


def _compute_distance(
    i: int,
    j: int,
    sequences: np.ndarray,
    metric: Metric,
    sep: str = "#",
    pad_to_max: bool = False,
) -> tuple:
    """Helper to compute distance between two sequences at indices i and j."""
    return (
        dist_func_long_sep(
            sequences[i],
            sequences[j],
            metric=metric,
            sep=sep,
            pad_to_max=pad_to_max,
        )
        if (pd.notnull(sequences[i]) and pd.notnull(sequences[j]))
        else 0
    )


def _compute_row(
    i: int,
    sequences: np.ndarray,
    n: int,
    metric: Metric,
    sep: str = "#",
    pad_to_max: bool = False,
) -> list:
    """Helper to compute a row of the distance matrix."""
    row_distances = [
        (
            _compute_distance(
                i,
                j,
                sequences,
                metric=metric,
                sep=sep,
                pad_to_max=pad_to_max,
            )
            if j > i
            else 0
        )
        for j in range(n)
    ]
    return row_distances
