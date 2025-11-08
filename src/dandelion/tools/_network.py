import multiprocessing
import math
import re

import networkx as nx
import numpy as np
import pandas as pd

from anndata import AnnData
from collections import defaultdict
from joblib import Parallel, delayed
from polyleven import levenshtein
from scanpy import logging as logg
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from tqdm import tqdm
from typing import Literal

try:
    from networkx.utils import np_random_state as random_state
except:
    from networkx.utils import random_state

from dandelion.utilities._core import Dandelion, Query
from dandelion.utilities._utilities import present, sanitize_data, Tree, FALSES


def generate_network(
    vdj_data: Dandelion,
    adata: AnnData | None = None,
    key: str | None = None,
    clone_key: str | None = None,
    min_size: int = 2,
    sample: int | None = None,
    force_replace: bool = False,
    verbose: bool = True,
    compute_layout: bool = True,
    layout_method: Literal["sfdp", "mod_fr"] = "sfdp",
    expanded_only: bool = False,
    use_existing_graph: bool = True,
    num_cores: int = 1,
    distance_mode: Literal["original", "full"] = "original",
    chunk_size: int | None = None,
    memory_limit_gb: float | None = None,
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
    compute_layout : bool, optional
        whether or not to generate the layout. May be time consuming if too many cells.
    layout_method : Literal["sfdp", "mod_fr"], optional
        accepts one of 'sfdp' or 'mod_fr'. 'sfdp' refers to `sfdp_layout` from `graph_tool` (C++ implementation; fast)
        whereas 'mod_fr' refers to modified Fruchterman-Reingold layout originally implemented in dandelion (python
        implementation; slow).
    expanded_only : bool, optional
        whether or not to only compute layout on expanded clonotypes.
    use_existing_graph : bool, optional
        whether or not to just compute the layout using the existing graph if it exists in the Dandelion object.
    num_cores : int, optional
        number of cores to use for parallelizable steps. -1 uses all available cores.
    distance_mode : Literal["original", "full"], optional
        method to compute distance matrix. 'original' refers to the original membership-based distance calculation
        whereas 'full' computes the full pairwise distance matrix using Dask delayed blocks.
    chunk_size : int | None, optional
        number of sequences to process as chunks for blockwise computation when using `distance_mode='full'`. If None,
        an automatic chunk size will be determined based on available memory and number of cores.
    memory_limit_gb: float | None, optional
        memory limit per worker in GB when using `distance_mode='full'`. If None, dask will control this automatically.
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
    if clone_key not in vdj_data.data:
        raise ValueError(
            "Data does not contain clone information. Please run find_clones."
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
                g, g_, lyt, lyt_ = _generate_layout(
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
        if key is None:
            key_ = "sequence_alignment_aa"
        else:
            key_ = key

        if key_ not in vdj_data.data:
            raise ValueError(f"key {key_} not found in data.")

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
        dat_seq.columns = [re.sub(key_ + "_", "", i) for i in dat_seq.columns]
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
        logg.info(f"Calculating distance matrix with mode = '{distance_mode}'")
        if distance_mode == "original":
            total_dist = calculate_distance_matrix_original(
                dat_seq, membership, verbose=verbose
            )
        elif distance_mode == "full":
            total_dist = calculate_distance_matrix_full(
                dat_seq,
                chunk_size=chunk_size,
                num_cores=num_cores,
                memory_limit_gb=memory_limit_gb,
                verbose=verbose,
            )
        else:
            raise ValueError(f"Unknown distance_mode: {distance_mode}")

        # follow the existing pipeline: make tmp dataframe, extract per-clone submatrices, mst, etc.
        df = pd.DataFrame(
            total_dist, index=dat_seq.index, columns=dat_seq.index
        )
        # reorder df's index and columns to match vdj_data.metadata index
        df = df.reindex(
            index=vdj_data.metadata.index, columns=vdj_data.metadata.index
        )

        # build cluster_dist as before (only include groups with >1 member)
        tmp_clusterdist2 = {}
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
            idxs = list(tmp_clusterdist2[c_])
            dist_mat_ = df.loc[idxs, idxs]
            s1, s2 = dist_mat_.shape
            if s1 > 1 and s2 > 1:
                cluster_dist[c_] = dist_mat_
        # Minimum spanning trees (unchanged)
        mst_tree = mst(cluster_dist, num_cores=num_cores, verbose=verbose)
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

        vertice_list = list(df.index)

        # final layout + graph creation (unchanged)
        g, g_, lyt, lyt_ = _generate_layout(
            vertices=vertice_list,
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
        " finished",
        time=start,
        deep=(
            "Updated Dandelion object: \n"
            "   'data', contig-indexed AIRR table\n"
            "   'metadata', cell-indexed observations table\n"
            "   'layout', graph layout\n"
            "   'graph', network constructed from distance matrices of VDJ- and VJ- chains"
        ),
    )

    # return or re-initialize vdj_data
    germline_ = getattr(vdj_data, "germline", None)
    if sample is not None:
        if (lyt and lyt_) is not None:
            out = Dandelion(
                data=dat_,
                metadata=vdj_data.metadata,
                clone_key=clone_key,
                layout=(lyt, lyt_),
                graph=(g, g_),
                germline=germline_,
                verbose=False,
                distances=(
                    csr_matrix(df.values) if regenerate else vdj_data.distances
                ),
            )
        else:
            out = Dandelion(
                data=dat_,
                metadata=vdj_data.metadata,
                clone_key=clone_key,
                graph=(g, g_),
                germline=germline_,
                verbose=False,
                distances=(
                    csr_matrix(df.values) if regenerate else vdj_data.distances
                ),
            )
        if adata is None:
            return out
        else:
            return out, adata
    else:
        if (lyt and lyt_) is not None:
            vdj_data.__init__(
                data=vdj_data.data,
                metadata=vdj_data.metadata,
                clone_key=clone_key,
                layout=(lyt, lyt_),
                graph=(g, g_),
                germline=germline_,
                initialize=False,
                verbose=False,
                distances=(
                    csr_matrix(df.values) if regenerate else vdj_data.distances
                ),
            )
        else:
            vdj_data.__init__(
                data=vdj_data.data,
                metadata=vdj_data.metadata,
                clone_key=clone_key,
                layout=None,
                graph=(g, g_),
                germline=germline_,
                initialize=False,
                verbose=False,
                distances=(
                    csr_matrix(df.values) if regenerate else vdj_data.distances
                ),
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


def calculate_distance_matrix_original(
    dat_seq: pd.DataFrame, membership: dict, verbose: bool = True
) -> np.ndarray:
    """
    Re-implementation of original membership-based distance calculation.

    Parameters
    ----------
    dat_seq : pd.DataFrame
        DataFrame with sequence columns; index corresponds to cell IDs (or whatever unique ids you use).
    membership : dict
        Mapping from clone_id -> list of indices (these indices must be present in dat_seq.index).
    verbose : bool
        Whether to show progress.

    Returns
    -------
    total_dist : np.ndarray (n x n)
        Aggregated distance matrix across all columns; diagonal set to NaN by caller.
    """
    n = dat_seq.shape[0]
    index_list = list(dat_seq.index)
    idx_to_pos = {idx: i for i, idx in enumerate(index_list)}

    # dmat_per_column will collect for each column a list of DataFrames (one per clone) to concat
    dmat_per_column = defaultdict(list)

    # iterate clones (membership) exactly like original
    iterator = membership
    if verbose:
        iterator = tqdm(
            membership,
            desc="Calculating distances (original, by membership)",
            disable=not verbose,
            bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
        )

    for t in iterator:
        tmp = dat_seq.loc[membership[t]]
        if tmp.shape[0] > 1:
            tmp = tmp.replace("[.]", "", regex=True)
            for col in tmp.columns:
                tdarray = np.array(tmp[col]).reshape(-1, 1)
                # keep the original NaN-check logic: return 0 when either is NaN
                d_mat_tmp = squareform(
                    pdist(
                        tdarray,
                        lambda x, y: (
                            levenshtein(x[0], y[0])
                            if (x[0] == x[0] and y[0] == y[0])
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


def _compute_block(chunk1, chunk2):
    """Helper used by Dask to compute distances between two 1-D arrays of strings."""
    return np.array([[levenshtein(a, b) for b in chunk2] for a in chunk1])


def calculate_distance_matrix_full(
    dat_seq: pd.DataFrame,
    chunk_size: int | None = None,
    num_cores: int = 1,
    memory_limit_gb: float | None = None,
    verbose: bool = True,
) -> np.ndarray:
    """
    Compute full pairwise distance matrix using Dask delayed blocks.

    Parameters
    ----------
    dat_seq: pd.DataFrame
        DataFrame of sequences.
    chunk_size: int | None, optional
        number of sequences to process as chunks for blockwise computation when using `distance_mode='full'`. If None,
        an automatic chunk size will be determined based on available memory and number of cores.
    num_cores: int, optional
        number of cores to use.
    memory_limit_gb: float | None, optional
        memory limit per worker in GB. If None, dask will control this automatically.
    verbose : bool, optional
        whether or not to show progress.

    Returns
    -------
    np.ndarray
        Aggregated distance matrix; diagonal set to NaN.
    """
    try:
        import dask
        import psutil

        from dask.distributed import Client, progress
        from dask.diagnostics import ProgressBar

        def _auto_chunk_size(
            n: int,
            num_cores: int,
            memory_limit_gb: float | None = None,
            safety_fraction=0.7,
        ) -> tuple[int, int]:
            """Compute dynamic chunk size to stay within memory budget per worker."""
            available_mem = psutil.virtual_memory().available / (1024**3)
            if memory_limit_gb is None:
                memory_limit_gb = available_mem * safety_fraction / num_cores
            mem_per_core = min(memory_limit_gb, available_mem / num_cores)
            # each element is 8 bytes; solve m^2 * 8 / 1024^3 ≈ mem_per_core
            chunk_size = int(math.sqrt((mem_per_core * (1024**3)) / 8))
            n_chunks = max(1, math.ceil(n / chunk_size))
            return chunk_size, n_chunks

    except ImportError:
        raise ImportError(
            "Please install dask, distributed and psutil: pip install dask[complete]"
        )
    n = dat_seq.shape[0]
    if chunk_size is None:
        chunk_size, n_chunks = _auto_chunk_size(
            n, num_cores=num_cores, memory_limit_gb=memory_limit_gb
        )
        if verbose:
            logg.info(f"Auto chunk size ≈ {chunk_size} sequences per block")
    n_chunks_eff = max(1, math.ceil(n / chunk_size))

    index_list = list(dat_seq.index)
    idx_to_pos = {idx: i for i, idx in enumerate(index_list)}

    client = None
    if num_cores > 1:
        if memory_limit_gb is None:
            client = Client(
                n_workers=num_cores,
                threads_per_worker=1,
                processes=False,
            )
        else:
            client = Client(
                n_workers=num_cores,
                threads_per_worker=1,
                processes=False,
                memory_limit=f"{memory_limit_gb}GB",
            )

    total_dist = np.zeros((n, n), dtype=float)

    for col in dat_seq.columns:
        seq_series = dat_seq[col].replace("[.]", "", regex=True)
        nonnull = seq_series.dropna()
        if nonnull.shape[0] <= 1:
            continue

        seq_indices = list(nonnull.index)
        seqs = nonnull.to_numpy(dtype=object)

        n_chunks_eff = max(1, min(n_chunks, len(seqs)))
        chunks = np.array_split(seqs, n_chunks_eff)
        chunk_sizes = [len(c) for c in chunks]
        cum = np.cumsum([0] + chunk_sizes)

        delayed_blocks = []
        block_positions = []

        for i in range(len(chunks)):
            for j in range(i, len(chunks)):
                if chunk_sizes[i] == 0 or chunk_sizes[j] == 0:
                    continue
                delayed_blocks.append(
                    dask.delayed(_compute_block)(chunks[i], chunks[j])
                )
                block_positions.append((i, j))

        # Compute blocks in parallel and show distributed progress
        if client is not None:
            futures = client.compute(delayed_blocks)
            if verbose:
                progress(futures)
            results = client.gather(futures)
        else:
            # fallback to single-threaded computation
            from dask import compute

            if verbose:
                from dask.diagnostics import ProgressBar

                with ProgressBar():
                    results = compute(*delayed_blocks, scheduler="threads")
            else:
                results = compute(*delayed_blocks, scheduler="threads")

        # Assemble temporary matrix
        m = len(seqs)
        tmp_block_mat = np.zeros((m, m), dtype=float)
        for (i, j), block in zip(block_positions, results):
            start_i, end_i = cum[i], cum[i + 1]
            start_j, end_j = cum[j], cum[j + 1]
            tmp_block_mat[start_i:end_i, start_j:end_j] = block
            if i != j:
                tmp_block_mat[start_j:end_j, start_i:end_i] = block.T

        pos_list = [idx_to_pos[idx] for idx in seq_indices]
        total_dist[np.ix_(pos_list, pos_list)] += tmp_block_mat

    if client is not None:
        client.close()

    np.fill_diagonal(total_dist, np.nan)
    return total_dist


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


def _generate_layout(
    vertices: list | None = None,
    edges: pd.DataFrame | None = None,
    min_size: int = 2,
    weight: str | None = None,
    verbose: bool = True,
    compute_layout: bool = True,
    layout_method: Literal["sfdp", "mod_fr"] = "sfdp",
    expanded_only: bool = False,
    graphs: tuple[nx.Graph, nx.Graph] = None,
    **kwargs,
) -> tuple[nx.Graph, nx.Graph, dict, dict]:
    """Generate layout.

    Parameters
    ----------
    vertices : list
        list of vertices
    edges : pd.DataFrame, optional
        edge list in a pandas data frame.
    min_size : int, optional
        minimum clone size.
    weight : str | None, optional
        name of weight column.
    verbose : bool, optional
        whether or not to print status
    compute_layout : bool, optional
        whether or not to compute layout.
    layout_method : Literal["sfdp", "mod_fr"], optional
        layout method.
    expanded_only : bool, optional
        whether or not to only compute layout on expanded clones.
    graphs: tuple[nx.Graph, nx.Graph], optional
        tuple of graphs.
    dist_mat : pd.DataFrame | None, optional
        distance matrix.
    **kwargs
        passed to fruchterman_reingold_layout.

    Returns
    -------
    tuple[nx.Graph, nx.Graph, dict, dict]
        graphs and layout positions.
    """
    if graphs is None:
        if vertices is not None:
            G = nx.Graph()
            G.add_nodes_from(vertices)
            if edges is not None:
                G.add_weighted_edges_from(
                    [
                        (x, y, z)
                        for x, y, z in zip(
                            edges["source"], edges["target"], edges["weight"]
                        )
                    ]
                )
        G_ = G.copy()
    else:
        G = graphs[0]
        G_ = graphs[1]
    if min_size == 2:
        if edges is not None:
            G_.remove_nodes_from(nx.isolates(G))
        else:
            pass
    elif min_size > 2:
        if edges is not None:
            for component in list(nx.connected_components(G_)):
                if len(component) < min_size:
                    for node in component:
                        G_.remove_node(node)
        else:
            pass
    if compute_layout:
        if layout_method == "mod_fr":
            if not expanded_only:
                if verbose:
                    print("Computing network layout")
                pos = _fruchterman_reingold_layout(G, weight=weight, **kwargs)
            else:
                pos = None
            if verbose:
                print("Computing expanded network layout")
            pos_ = _fruchterman_reingold_layout(G_, weight=weight, **kwargs)
        elif layout_method == "sfdp":
            try:
                from graph_tool.all import sfdp_layout
            except ImportError:
                print(
                    "To benefit from faster layout computation, please install graph-tool: "
                    "conda install -c conda-forge graph-tool"
                )
                nographtool = True
            if "nographtool" in locals():
                if not expanded_only:
                    if verbose:
                        print("Computing network layout")
                    pos = _fruchterman_reingold_layout(
                        G, weight=weight, **kwargs
                    )
                else:
                    pos = None
                if verbose:
                    print("Computing expanded network layout")
                pos_ = _fruchterman_reingold_layout(G_, weight=weight, **kwargs)
            else:
                if not expanded_only:
                    gtg = nx2gt(G)
                    if verbose:
                        print("Computing network layout")
                    posx = sfdp_layout(gtg, **kwargs)
                    pos = dict(
                        zip(list(gtg.vertex_properties["id"]), list(posx))
                    )
                else:
                    pos = None
                gtg_ = nx2gt(G_)
                if verbose:
                    print("Computing expanded network layout")
                posx_ = sfdp_layout(gtg_, **kwargs)
                pos_ = dict(
                    zip(list(gtg_.vertex_properties["id"]), list(posx_))
                )
        if pos is None:
            G = G_
            pos = pos_

        return (G, G_, pos, pos_)
    else:
        return (G, G_, None, None)


# when dealing with a lot of unconnected vertices, the pieces fly out to infinity and the original fr layout can't be
# used
# work around from https://stackoverflow.com/questions/14283341/how-to-increase-node-spacing-for-networkx-spring-layout
# code chunk from networkx's layout.py https://github.com/networkx/networkx/blob/master/networkx/drawing/layout.py


def _process_params(G, center, dim):
    """Some boilerplate code."""
    if not isinstance(G, nx.Graph):
        empty_graph = nx.Graph()
        empty_graph.add_nodes_from(G)
        G = empty_graph

    if center is None:
        center = np.zeros(dim)
    else:
        center = np.asarray(center)

    if len(center) != dim:
        msg = "length of center coordinates must match dimension of layout"
        raise ValueError(msg)

    return G, center


def _fruchterman_reingold_layout(
    G,
    k=None,
    pos=None,
    fixed=None,
    iterations=50,
    threshold=1e-4,
    weight="weight",
    scale=1,
    center=None,
    dim=2,
    seed=None,
    **kwargs,
):
    """
    Position nodes using Fruchterman-Reingold force-directed algorithm.

    The algorithm simulates a force-directed representation of the network
    treating edges as springs holding nodes close, while treating nodes
    as repelling objects, sometimes called an anti-gravity force.
    Simulation continues until the positions are close to an equilibrium.
    There are some hard-coded values: minimal distance between
    nodes (0.01) and "temperature" of 0.1 to ensure nodes don't fly away.
    During the simulation, `k` helps determine the distance between nodes,
    though `scale` and `center` determine the size and place after
    rescaling occurs at the end of the simulation.
    Fixing some nodes doesn't allow them to move in the simulation.
    It also turns off the rescaling feature at the simulation's end.
    In addition, setting `scale` to `None` turns off rescaling.

    Parameters
    ----------
    G : NetworkX graph or list of nodes
        A position will be assigned to every node in G.
    k : float (default=None)
        Optimal distance between nodes.  If None the distance is set to
        1/sqrt(n) where n is the number of nodes.  Increase this value
        to move nodes farther apart.
    pos : dict or None  Optional (default=None)
        Initial positions for nodes as a dictionary with node as keys
        and values as a coordinate list or tuple.  If None, then use
        random initial positions.
    fixed : list or None  Optional (default=None)
        Nodes to keep fixed at initial position.
        ValueError raised if `fixed` specified and `pos` not.
    iterations : int  Optional (default=50)
        Maximum number of iterations taken
    threshold: float Optional (default = 1e-4)
        Threshold for relative error in node position changes.
        The iteration stops if the error is below this threshold.
    weight : string or None   Optional (default='weight')
        The edge attribute that holds the numerical value used for
        the edge weight.  If None, then all edge weights are 1.
    scale : number or None (default: 1)
        Scale factor for positions. Not used unless `fixed is None`.
        If scale is None, no rescaling is performed.
    center : array-like or None
        Coordinate pair around which to center the layout.
        Not used unless `fixed is None`.
    dim : int
        Dimension of layout.
    seed : int, RandomState instance or None  Optional (default=None)
        Set the random state for deterministic node layouts.
        If int, `seed` is the seed used by the random number generator,
        if numpy.random.RandomState instance, `seed` is the random
        number generator,
        if None, the random number generator is the RandomState instance used
        by numpy.random.

    Returns
    -------
    pos : dict
        A dictionary of positions keyed by node

    Examples
    --------
    >>> G = nx.path_graph(4)
    >>> pos = nx.spring_layout(G)
    # The same using longer but equivalent function name
    >>> pos = nx.fruchterman_reingold_layout(G)
    """
    G, center = _process_params(G, center, dim)

    if fixed is not None:
        if pos is None:
            raise ValueError("nodes are fixed without positions given")
        for node in fixed:
            if node not in pos:
                raise ValueError("nodes are fixed without positions given")
        nfixed = {node: i for i, node in enumerate(G)}
        fixed = np.asarray([nfixed[node] for node in fixed])

    if pos is not None:
        # Determine size of existing domain to adjust initial positions
        dom_size = max(coord for pos_tup in pos.values() for coord in pos_tup)
        if dom_size == 0:
            dom_size = 1
        pos_arr = seed.rand(len(G), dim) * dom_size + center

        for i, n in enumerate(G):
            if n in pos:
                pos_arr[i] = np.asarray(pos[n])
    else:
        pos_arr = None
        dom_size = 1

    if len(G) == 0:
        return {}
    if len(G) == 1:
        return {nx.utils.arbitrary_element(G.nodes()): center}

    try:
        # Sparse matrix
        if len(G) < 500:  # sparse solver for large graphs
            raise ValueError
        if int(nx.__version__[0]) > 2:
            A = nx.to_scipy_sparse_array(G, weight=weight, dtype="f")
        else:
            A = nx.to_scipy_sparse_matrix(G, weight=weight, dtype="f")
        if k is None and fixed is not None:
            # We must adjust k by domain size for layouts not near 1x1
            nnodes, _ = A.shape
            k = dom_size / np.sqrt(nnodes)
        pos = _sparse_fruchterman_reingold(
            A, k, pos_arr, fixed, iterations, threshold, dim, seed
        )
    except ValueError:
        A = nx.to_numpy_array(G, weight=weight)
        if k is None and fixed is not None:
            # We must adjust k by domain size for layouts not near 1x1
            nnodes, _ = A.shape
            k = dom_size / np.sqrt(nnodes)
        pos = _fruchterman_reingold(
            A, k, pos_arr, fixed, iterations, threshold, dim, seed
        )
    if fixed is None and scale is not None:
        pos = _rescale_layout(pos, scale=scale) + center
    pos = dict(zip(G, pos))
    return pos


@random_state(7)
def _fruchterman_reingold(
    A,
    k=None,
    pos=None,
    fixed=None,
    iterations=50,
    threshold=1e-4,
    dim=2,
    seed=None,
    **kwargs,
):
    """Fruchterman Reingold algorithm."""
    # Position nodes in adjacency matrix A using Fruchterman-Reingold
    # Entry point for NetworkX graph is fruchterman_reingold_layout()
    import numpy as np

    try:
        nnodes, _ = A.shape
    except AttributeError as e:
        msg = "fruchterman_reingold() takes an adjacency matrix as input"
        raise nx.NetworkXError(msg) from e

    if pos is None:
        # random initial positions
        pos = np.asarray(seed.rand(nnodes, dim), dtype=A.dtype)
    else:
        # make sure positions are of same type as matrix
        pos = pos.astype(A.dtype)

    # optimal distance between nodes
    if k is None:
        k = np.sqrt(1.0 / nnodes)
    # the initial "temperature"  is about .1 of domain area (=1x1)
    # this is the largest step allowed in the dynamics.
    # We need to calculate this in case our fixed positions force our domain
    # to be much bigger than 1x1
    t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1])) * 0.1
    # simple cooling scheme.
    # linearly step down by dt on each iteration so last iteration is size dt.
    dt = t / float(iterations + 1)
    delta = np.zeros((pos.shape[0], pos.shape[0], pos.shape[1]), dtype=A.dtype)
    # the inscrutable (but fast) version
    # this is still O(V^2)
    # could use multilevel methods to speed this up significantly
    for iteration in range(iterations):
        # matrix of difference between points
        delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
        # distance between points
        distance = np.linalg.norm(delta, axis=-1)
        # enforce minimum distance of 0.01
        np.clip(distance, 0.001, None, out=distance)
        # displacement "force"
        displacement = np.einsum(
            "ijk,ij->ik", delta, (k * k / distance**2 - A * distance / k)
        )
        displacement = displacement - pos / (k * np.sqrt(nnodes))
        # update positions
        length = np.linalg.norm(displacement, axis=-1)
        length = np.where(length < 0.01, 0.1, length)
        delta_pos = np.einsum("ij,i->ij", displacement, t / length)
        if fixed is not None:
            # don't change positions of fixed nodes
            delta_pos[fixed] = 0.0
        pos += delta_pos
        # cool temperature
        t -= dt
        err = np.linalg.norm(delta_pos) / nnodes
        if err < threshold:
            break
    return pos


@random_state(7)
def _sparse_fruchterman_reingold(
    A,
    k=None,
    pos=None,
    fixed=None,
    iterations=50,
    threshold=1e-4,
    dim=2,
    seed=None,
    **kwargs,
):
    """Sparse Fruchterman Reingold algorithm."""
    # Position nodes in adjacency matrix A using Fruchterman-Reingold
    # Entry point for NetworkX graph is fruchterman_reingold_layout()
    # Sparse version
    import numpy as np

    try:
        nnodes, _ = A.shape
    except AttributeError as e:
        msg = "fruchterman_reingold() takes an adjacency matrix as input"
        raise nx.NetworkXError(msg) from e
    try:
        from scipy.sparse import coo_matrix
    except ImportError as e:
        msg = "_sparse_fruchterman_reingold() scipy numpy: http://scipy.org/ "
        raise ImportError(msg) from e
    # make sure we have a LIst of Lists representation
    try:
        A = A.tolil()
    except AttributeError:
        A = (coo_matrix(A)).tolil()

    if pos is None:
        # random initial positions
        pos = np.asarray(seed.rand(nnodes, dim), dtype=A.dtype)
    else:
        # make sure positions are of same type as matrix
        pos = pos.astype(A.dtype)

    # no fixed nodes
    if fixed is None:
        fixed = []

    # optimal distance between nodes
    if k is None:
        k = np.sqrt(1.0 / nnodes)
    # the initial "temperature"  is about .1 of domain area (=1x1)
    # this is the largest step allowed in the dynamics.
    t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1])) * 0.1
    # simple cooling scheme.
    # linearly step down by dt on each iteration so last iteration is size dt.
    dt = t / float(iterations + 1)

    displacement = np.zeros((dim, nnodes))
    for iteration in range(iterations):
        displacement *= 0
        # loop over rows
        for i in range(A.shape[0]):
            if i in fixed:
                continue
            # difference between this row's node position and all others
            delta = (pos[i] - pos).T
            # distance between points
            distance = np.sqrt((delta**2).sum(axis=0))
            # enforce minimum distance of 0.01
            distance = np.where(distance < 0.01, 0.01, distance)
            # the adjacency matrix row
            Ai = np.asarray(A.getrowview(i).toarray())
            # displacement "force"
            displacement[:, i] += (
                delta * (k * k / distance**2 - Ai * distance / k)
            ).sum(axis=1)
        displacement = displacement - pos / (k * np.sqrt(nnodes))
        # update positions
        length = np.sqrt((displacement**2).sum(axis=0))
        length = np.where(length < 0.01, 0.1, length)
        delta_pos = (displacement * t / length).T
        pos += delta_pos
        # cool temperature
        t -= dt
        err = np.linalg.norm(delta_pos) / nnodes
        if err < threshold:
            break
    return pos


def _rescale_layout(pos, scale=1):
    """
    Return scaled position array to (-scale, scale) in all axes.

    The function acts on NumPy arrays which hold position information.
    Each position is one row of the array. The dimension of the space
    equals the number of columns. Each coordinate in one column.
    To rescale, the mean (center) is subtracted from each axis separately.
    Then all values are scaled so that the largest magnitude value
    from all axes equals `scale` (thus, the aspect ratio is preserved).
    The resulting NumPy Array is returned (order of rows unchanged).

    Parameters
    ----------
    pos : numpy array
        positions to be scaled. Each row is a position.
    scale : number (default: 1)
        The size of the resulting extent in all directions.

    Returns
    -------
    pos : numpy array
        scaled positions. Each row is a position.
    """
    # Find max length over all dimensions
    lim = 0  # max coordinate for all axes
    for i in range(pos.shape[1]):
        pos[:, i] -= pos[:, i].mean()
        lim = max(abs(pos[:, i]).max(), lim)
    # rescale to (-scale, scale) in all directions, preserves aspect
    if lim > 0:
        for i in range(pos.shape[1]):
            pos[:, i] *= scale / lim
    return pos


def extract_edge_weights(
    vdj_data: Dandelion, expanded_only: bool = False
) -> list:
    """
    Retrieve edge weights (BCR levenshtein distance) from graph.

    Parameters
    ----------
    vdj_data : Dandelion
        Dandelion object after `tl.generate_network` has been run.
    expanded_only : bool, optional
        whether to retrieve the edge weights from the expanded only graph or entire graph.

    Returns
    -------
    list
        list of edge weights.
    """
    if expanded_only:
        try:
            edges, weights = zip(
                *nx.get_edge_attributes(vdj_data.graph[1], "weight").items()
            )
        except ValueError as e:
            print(
                "{} i.e. the graph does not contain edges. Therefore, edge weights not returned.".format(
                    e
                )
            )
    else:
        try:
            edges, weights = zip(
                *nx.get_edge_attributes(vdj_data.graph[0], "weight").items()
            )
        except ValueError as e:
            print(
                "{} i.e. the graph does not contain edges. Therefore, edge weights not returned.".format(
                    e
                )
            )
    if "weights" in locals():
        return weights


# from https://bbengfort.github.io/2016/06/graph-tool-from-networkx/
def nx2gt(nxG):
    """Convert a networkx graph to a graph-tool graph."""
    try:
        import graph_tool as gt
    except ImportError:
        raise ImportError(
            "Please install graph-tool: conda install -c conda-forge graph-tool"
        )
    # Phase 0: Create a directed or undirected graph-tool Graph
    gtG = gt.Graph(directed=nxG.is_directed())

    # Add the Graph properties as "internal properties"
    for key, value in nxG.graph.items():
        # Convert the value and key into a type for graph-tool
        tname, value, key = get_prop_type(value, key)

        prop = gtG.new_graph_property(tname)  # Create the PropertyMap
        gtG.graph_properties[key] = prop  # Set the PropertyMap
        gtG.graph_properties[key] = value  # Set the actual value

    # Phase 1: Add the vertex and edge property maps
    # Go through all nodes and edges and add seen properties
    # Add the node properties first
    nprops = set()  # cache keys to only add properties once
    for node, data in list(nxG.nodes(data=True)):
        # Go through all the properties if not seen and add them.
        for key, val in data.items():
            if key in nprops:
                continue  # Skip properties already added

            # Convert the value and key into a type for graph-tool
            tname, _, key = get_prop_type(val, key)

            prop = gtG.new_vertex_property(tname)  # Create the PropertyMap
            gtG.vertex_properties[key] = prop  # Set the PropertyMap

            # Add the key to the already seen properties
            nprops.add(key)

    # Also add the node id: in NetworkX a node can be any hashable type, but
    # in graph-tool node are defined as indices. So we capture any strings
    # in a special PropertyMap called 'id' -- modify as needed!
    gtG.vertex_properties["id"] = gtG.new_vertex_property("string")

    # Add the edge properties second
    eprops = set()  # cache keys to only add properties once
    for src, dst, data in list(nxG.edges(data=True)):
        # Go through all the edge properties if not seen and add them.
        for key, val in data.items():
            if key in eprops:
                continue  # Skip properties already added

            # Convert the value and key into a type for graph-tool
            tname, _, key = get_prop_type(val, key)

            prop = gtG.new_edge_property(tname)  # Create the PropertyMap
            gtG.edge_properties[key] = prop  # Set the PropertyMap

            # Add the key to the already seen properties
            eprops.add(key)

    # Phase 2: Actually add all the nodes and vertices with their properties
    # Add the nodes
    vertices = {}  # vertex mapping for tracking edges later
    for node, data in list(nxG.nodes(data=True)):
        # Create the vertex and annotate for our edges later
        v = gtG.add_vertex()
        vertices[node] = v

        # Set the vertex properties, not forgetting the id property
        data["id"] = str(node)
        for key, value in data.items():
            gtG.vp[key][v] = value  # vp is short for vertex_properties

    # Add the edges
    for src, dst, data in list(nxG.edges(data=True)):
        # Look up the vertex structs from our vertices mapping and add edge.
        e = gtG.add_edge(vertices[src], vertices[dst])

        # Add the edge properties
        for key, value in data.items():
            gtG.ep[key][e] = value  # ep is short for edge_properties

    # Done, finally!
    return gtG


def get_prop_type(value, key=None):
    """
    Perform typing and value conversion for the graph_tool PropertyMap class.

    If a key is provided, it also ensures the key is in a format that can be
    used with the PropertyMap. Returns a tuple, (type name, value, key)
    """
    # Deal with the value
    if isinstance(value, bool):
        tname = "bool"

    elif isinstance(value, int):
        tname = "float"
        value = float(value)

    elif isinstance(value, float):
        tname = "float"

    elif isinstance(value, dict):
        tname = "object"

    else:
        tname = "string"
        value = str(value)

    return tname, value, key
