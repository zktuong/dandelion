from __future__ import annotations
import math
import zarr
import zarr.storage

import dask.array as da
import numpy as np
import pandas as pd

from dask import compute
from dask.distributed import Client, progress, Lock
from pathlib import Path
from scanpy import logging as logg
from typing import Callable
from zarr.codecs import BloscCodec

from dandelion.utilities._utilities import dist_func_long_sep


def calculate_distance_matrix_zarr(
    dat_seq: pd.DataFrame,
    func: Callable,
    pad_to_max: bool = False,
    membership: dict | None = None,
    zarr_path: Path | str | None = None,
    chunk_size: int | None = None,
    num_cores: int = 1,
    memory_limit_gb: float | None = None,
    memory_safety_fraction: float = 0.3,
    compress: bool = True,
    lazy: bool = False,
    verbose: bool = True,
) -> np.ndarray | da.Array:
    """
    Compute full pairwise distance matrix writing directly to Zarr.

    This approach:
    - Uses Dask to parallelize distance computation
    - Writes results directly to compressed Zarr storage AS blocks complete
    - Never materializes full matrix in memory
    - Suitable for 100k-10M sequences with compression

    The output matrix will have the same row/column order as dat_seq.index.
    Reorder dat_seq before calling this function if you need a specific order.

    Parameters
    ----------
    dat_seq : pd.DataFrame
        DataFrame of sequences (rows=samples, columns=different sequence types).
        Row order determines the order of the distance matrix.
    func : Callable
        Function to compute distance between two sequences.
    pad_to_max : bool, optional
        Whether to pad sequences to maximum length before distance calculation.
    zarr_path : str, optional
        Path to save Zarr array (directory). If None, uses 'distance_matrix.zarr'.
    chunk_size : int | None, optional
        Number of sequences per computation chunk. If None, auto-computed.
    num_cores : int
        Number of cores/workers to use
    memory_limit_gb : float | None, optional
        Memory limit per worker in GB
    memory_safety_fraction : float
        Fraction of available memory to use (conservative for Dask overhead)
    chunk_size : int | None, optional
        Chunk size for Zarr storage (can differ from computation chunk_size).
        Larger = better compression, smaller = better random access.
        If None, uses same as chunk_size.
    compress : bool, optional
        Whether to compress the Zarr array using Blosc with zstd.
    lazy: bool, optional
        If True, computation will be performed lazily using Dask/Zarr arrays. True will also return a Dask array view of the
        distance matrix stored on disk instead of a numpy array stored in memory.
    verbose : bool
        Whether to show progress

    Returns
    -------
    np.ndarray | da.Array
        Distance matrix in requested format
    """
    # Try import the dask dependencies otherwise raise error and tell user to install dask[complete]
    try:
        import dask.array as da
    except ImportError as e:
        raise ImportError(
            "Please install dask with the complete set of dependencies: pip install 'dask[complete]'"
        )
    if zarr_path is None:
        zarr_path = "distance_matrix.zarr"

    n = dat_seq.shape[0]

    # Determine computation chunk size
    if chunk_size is None:
        chunk_size, _ = _auto_chunk_size(
            n,
            num_cores=num_cores,
            memory_limit_gb=memory_limit_gb,
            safety_fraction=memory_safety_fraction,
        )
        if verbose:
            logg.info(
                f"Auto computation chunk size: {chunk_size} sequences per block"
            )
            logg.info(
                f"Estimated memory per block: {(chunk_size**2 * 8) / (1024**3):.2f} GB"
            )
    else:
        if verbose:
            logg.info(f"Computation chunk size: {chunk_size}")

    if verbose:
        logg.info(f"Zarr storage chunk size: {chunk_size}")
        logg.info(f"Output size: {(n * n * 8) / (1024**3):.2f} GB uncompressed")

    # Setup compressor
    comp = (
        BloscCodec(cname="zstd", clevel=3, shuffle="bitshuffle")
        if compress
        else None
    )

    # Create Zarr array
    store = zarr.storage.LocalStore(zarr_path)
    root = zarr.open_group(store=store, mode="w")

    if comp is not None:
        z_array = root.create_array(
            "distance_matrix",
            shape=(n, n),
            chunks=(chunk_size, chunk_size),
            dtype="float64",
            compressors=[comp],
            fill_value=0.0,
        )
    else:
        z_array = root.create_array(
            "distance_matrix",
            shape=(n, n),
            chunks=(chunk_size, chunk_size),
            dtype="float64",
            fill_value=0.0,
        )

    # zarr_lock = Lock("zarr_distance_matrix_lock") if num_cores > 1 else None
    zarr_lock_name = "zarr_distance_matrix_lock" if num_cores > 1 else None

    # Store metadata
    z_array.attrs["n_sequences"] = n
    z_array.attrs["chunk_size"] = chunk_size
    z_array.attrs["index"] = list(dat_seq.index)
    z_array.attrs["columns"] = list(dat_seq.columns)

    if verbose:
        logg.info(f"\nCreated Zarr array at: {zarr_path}")
        logg.info(f"Compressor: {comp}")

    # Setup Dask client
    client = _setup_dask_client(num_cores, memory_limit_gb, verbose)

    try:
        # Compute distances and write blocks as they complete
        _ = _compute_multicol_distances_streaming(
            dat_seq=dat_seq,
            func=func,
            pad_to_max=pad_to_max,
            chunk_size=chunk_size,
            z_array=z_array,
            lock_name=zarr_lock_name,
            client=client,
            membership=membership,
            serialise=False if num_cores > 1 else True,
            verbose=verbose,
        )

        # Set diagonal to NaN
        for i in range(0, n, chunk_size):
            end = min(i + chunk_size, n)
            diag_block = z_array[i:end, i:end]
            np.fill_diagonal(diag_block, np.nan)
            z_array[i:end, i:end] = diag_block

        if verbose:
            logg.info(f"\n{'='*60}\n")
            logg.info("✓ Distance matrix complete!")
            logg.info("\nFinal array info:")
            logg.info(f"  Shape: {z_array.shape}")
            logg.info(f"  Dtype: {z_array.dtype}")
            logg.info(f"  Chunks: {z_array.chunks}")
            logg.info(
                f"  Size on disk: {z_array.nbytes / (1024**3):.2f} GB uncompressed"
            )

            try:
                actual_size = sum(
                    f.stat().st_size
                    for f in Path(zarr_path).rglob("*")
                    if f.is_file()
                )
                logg.info(
                    f"  Actual size: {actual_size / (1024**3):.2f} GB (with compression)"
                )
                compression_ratio = (
                    z_array.nbytes / actual_size if actual_size > 0 else 0
                )
                logg.info(f"  Compression ratio: {compression_ratio:.1f}x")
            except Exception:
                pass

        # Return based on return_type
        if lazy:
            return da.from_zarr(z_array)
        else:
            return z_array[:]

    finally:
        if client is not None:
            client.close()


def _compute_multicol_distances_streaming(
    dat_seq: pd.DataFrame,
    func: Callable,
    pad_to_max: bool,
    chunk_size: int,
    z_array: zarr.Array,
    lock_name: str | None,
    client: Client | None = None,
    membership: dict | None = None,
    serialise: bool = False,
    verbose: bool = True,
):
    """
    Compute distance matrix using concatenation across all columns,
    writing blocks directly to Zarr as they complete.
    """
    try:
        import dask
    except ImportError:
        raise ImportError(
            "Please install dask with the complete set of dependencies: pip install 'dask[complete]'"
        )
    # Clean all sequences
    dat_seq_clean = (
        dat_seq.replace("[.]", "", regex=True).fillna("").replace("None", "")
    )

    seqs = dat_seq_clean.to_numpy(dtype=object)
    m = len(seqs)

    if m <= 1:
        return None

    delayed_blocks = []
    if membership is not None:
        index_list = list(dat_seq_clean.index)
        index_map = {idx: i for i, idx in enumerate(index_list)}  # fast lookup
        delayed_blocks = []
        # reduce to only clonotypes with >1 member
        membership = {k: v for k, v in membership.items() if len(v) > 1}
        if verbose:
            logg.info(
                f"Computing distances within {len(membership)} clonotypes with more than 1 member..."
            )
        for _, members in membership.items():
            # # check if len(members) > chunk_size to decide whether to further chunk
            # if len(members) > chunk_size:
            #     # split members into even subchunks
            #     n_subchunks = max(1, math.ceil(len(members) / chunk_size))
            #     member_chunks = np.array_split(members, n_subchunks)
            #     chunk_sizes = [len(c) for c in member_chunks]
            #     # iterate i/j in member_chunks:
            #     for i, member_chunk_i in enumerate(member_chunks):
            #         for j, member_chunk_j in enumerate(
            #             # so that off-diagonal blocks (i<j) are computed only once
            #             # is_diagonal should take care of transposing it later.
            #             member_chunks[i:],
            #             start=i,
            #         ):
            #             tmp_i = dat_seq_clean.loc[member_chunk_i]
            #             tmp_j = dat_seq_clean.loc[member_chunk_j]
            #             seqs_i = tmp_i.to_numpy(dtype=object)
            #             seqs_j = tmp_j.to_numpy(dtype=object)
            #             idxs_i = [index_map[m] for m in member_chunk_i]
            #             idxs_j = [index_map[m] for m in member_chunk_j]
            #             is_diagonal = i == j
            #             if serialise:
            #                 delayed_blocks.append(
            #                     _compute_and_write_membership(
            #                         chunk_i=seqs_i,
            #                         chunk_j=seqs_j,
            #                         func=func,
            #                         pad_to_max=pad_to_max,
            #                         z_array=z_array,
            #                         idxs_i=idxs_i,
            #                         idxs_j=idxs_j,
            #                         is_diagonal=is_diagonal,
            #                         lock_name=lock_name,
            #                     )
            #                 )
            #             else:
            #                 delayed_blocks.append(
            #                     dask.delayed(_compute_and_write_membership)(
            #                         chunk_i=seqs_i,
            #                         chunk_j=seqs_j,
            #                         func=func,
            #                         pad_to_max=pad_to_max,
            #                         z_array=z_array,
            #                         idxs_i=idxs_i,
            #                         idxs_j=idxs_j,
            #                         is_diagonal=is_diagonal,
            #                         lock_name=lock_name,
            #                     )
            #                 )
            # else:
            #     tmp = dat_seq_clean.loc[members]
            #     seqs_tmp = tmp.to_numpy(dtype=object)
            #     idxs = [index_map[m] for m in members]
            #     if serialise:
            #         delayed_blocks.append(
            #             _compute_and_write_membership(
            #                 chunk_i=seqs_tmp,
            #                 chunk_j=seqs_tmp,
            #                 func=func,
            #                 pad_to_max=pad_to_max,
            #                 z_array=z_array,
            #                 idxs_i=idxs,
            #                 idxs_j=idxs,
            #                 is_diagonal=True,  # diagonal (within group)
            #                 lock_name=lock_name,
            #             )
            #         )
            #     else:
            #         delayed_blocks.append(
            #             dask.delayed(_compute_and_write_membership)(
            #                 chunk_i=seqs_tmp,
            #                 chunk_j=seqs_tmp,
            #                 func=func,
            #                 pad_to_max=pad_to_max,
            #                 z_array=z_array,
            #                 idxs_i=idxs,
            #                 idxs_j=idxs,
            #                 is_diagonal=True,  # diagonal (within group)
            #                 lock_name=lock_name,
            #             )
            #         )
            if serialise:
                delayed_blocks.append(
                    _compute_and_write_membership(
                        members=members,
                        chunk_size=chunk_size,
                        seqs_df=dat_seq_clean.loc[members],
                        index_map=index_map,
                        func=func,
                        pad_to_max=pad_to_max,
                        z_array=z_array,
                        lock_name=lock_name,
                    )
                )
            else:
                delayed_blocks.append(
                    dask.delayed(_compute_and_write_membership)(
                        members=members,
                        chunk_size=chunk_size,
                        seqs_df=dat_seq_clean.loc[members],
                        index_map=index_map,
                        func=func,
                        pad_to_max=pad_to_max,
                        z_array=z_array,
                        lock_name=lock_name,
                    )
                )
    else:
        if verbose:
            logg.info(
                f"Processing {m} sequences across {dat_seq_clean.shape[1]} columns..."
            )
        # Determine number of chunks
        n_chunks = max(1, math.ceil(m / chunk_size))
        chunks_list = np.array_split(seqs, n_chunks)
        chunk_sizes = [len(c) for c in chunks_list]
        cum_sizes = np.cumsum([0] + chunk_sizes)
        for i in range(len(chunks_list)):
            for j in range(i, len(chunks_list)):
                if chunk_sizes[i] == 0 or chunk_sizes[j] == 0:
                    continue
                start_i, end_i = cum_sizes[i], cum_sizes[i + 1]
                start_j, end_j = cum_sizes[j], cum_sizes[j + 1]
                is_diagonal = i == j
                if serialise:
                    delayed_blocks.append(
                        _compute_and_write_block(
                            chunk_i=chunks_list[i],
                            chunk_j=chunks_list[j],
                            func=func,
                            pad_to_max=pad_to_max,
                            z_array=z_array,
                            start_i=start_i,
                            end_i=end_i,
                            start_j=start_j,
                            end_j=end_j,
                            is_diagonal=is_diagonal,
                            lock_name=lock_name,
                        )
                    )
                else:
                    delayed_blocks.append(
                        dask.delayed(_compute_and_write_block)(
                            chunk_i=chunks_list[i],
                            chunk_j=chunks_list[j],
                            func=func,
                            pad_to_max=pad_to_max,
                            z_array=z_array,
                            start_i=start_i,
                            end_i=end_i,
                            start_j=start_j,
                            end_j=end_j,
                            is_diagonal=is_diagonal,
                            lock_name=lock_name,
                        )
                    )

    if verbose:
        logg.info(f"Computing and writing {len(delayed_blocks)} blocks...")
    # Compute blocks - they write to Zarr as they complete
    if client is not None:
        futures = client.compute(delayed_blocks)
        if verbose:
            progress(futures)
        client.gather(futures)
    else:
        if verbose:
            from dask.diagnostics import ProgressBar

            with ProgressBar():
                compute(*delayed_blocks, scheduler="threads")
        else:
            compute(*delayed_blocks, scheduler="threads")

    return True


def _compute_block_multicol(
    seqs_i: np.ndarray,
    seqs_j: np.ndarray,
    func: Callable,
    pad_to_max: bool,
):
    """
    Compute pairwise distances between two sequence chunks across multiple columns.
    seqs_i and seqs_j are 2D arrays: shape (n_rows, n_cols)
    """
    n_i, n_j = len(seqs_i), len(seqs_j)
    dist_block = np.zeros((n_i, n_j), dtype=float)

    for i, row_i in enumerate(seqs_i):
        for j, row_j in enumerate(seqs_j):
            dist_block[i, j] = dist_func_long_sep(
                row_i, row_j, func=func, pad_to_max=pad_to_max
            )

    return dist_block


def dask_safe_slice_square(arr: da.Array, pos: list) -> da.Array:
    """Return arr[pos, pos], Dask-safe."""
    if isinstance(arr, da.Array):
        return arr[pos, :][:, pos]
    else:
        return arr[np.ix_(pos, pos)]


def _compute_and_write_block(
    chunk_i: np.ndarray,
    chunk_j: np.ndarray,
    func: Callable,
    pad_to_max: bool,
    z_array: zarr.Array,
    start_i: int,
    end_i: int,
    start_j: int,
    end_j: int,
    is_diagonal: bool,
    lock_name: str | None,
):
    """
    Compute a single block and write it directly to Zarr.
    """
    lock = Lock(lock_name) if lock_name is not None else None
    block = _compute_block_multicol(chunk_i, chunk_j, func, pad_to_max)

    if lock is not None:
        with lock:
            # Write to Zarr immediately
            z_array[start_i:end_i, start_j:end_j] += block

            # If not diagonal block, write transpose too
            if not is_diagonal:
                z_array[start_j:end_j, start_i:end_i] += block.T
    else:
        # Write to Zarr immediately
        z_array[start_i:end_i, start_j:end_j] += block

        # If not diagonal block, write transpose too
        if not is_diagonal:
            z_array[start_j:end_j, start_i:end_i] += block.T

    return True  # Just return success flag


# def _compute_and_write_membership(
#     chunk_i: np.ndarray,
#     chunk_j: np.ndarray,
#     func: Callable,
#     pad_to_max: bool,
#     z_array: zarr.Array,
#     idxs_i: list[int],
#     idxs_j: list[int],
#     is_diagonal: bool,
#     lock: Lock | None,
# ):
#     """
#     Compute a block for specific row/column indices (e.g., membership subsets)
#     and write it directly to Zarr.
#     """
#     # Compute the block using the distance function
#     block = _compute_block_multicol(chunk_i, chunk_j, func, pad_to_max)
#     if lock is not None:
#         with lock:
#             # Use np.ix_ to assign to the correct indices in the zarr array
#             z_array[np.ix_(idxs_i, idxs_j)] = block

#             # If not diagonal, write transpose to the symmetric location
#             if not is_diagonal:
#                 z_array[np.ix_(idxs_j, idxs_i)] = block.T
#     else:
#         z_array[np.ix_(idxs_i, idxs_j)] = block
#         if not is_diagonal:
#             z_array[np.ix_(idxs_j, idxs_i)] = block.T

#     return True


def _compute_and_write_membership(
    members: list,
    chunk_size: int,
    seq_df: pd.DataFrame,
    index_map: dict,
    func: Callable,
    pad_to_max: bool,
    z_array: zarr.Array,
    lock_name: str | None,
):
    """
    Extract sequences for members and delegate to computation function
    that handles chunking internally.
    """
    # Extract sequences and indices for all members
    # tmp = dat_seq_clean.loc[members]

    seqs = seq_df.to_numpy(dtype=object)
    idxs = [index_map[m] for m in members]

    # Delegate to computation function that handles chunking
    _compute_and_write_with_chunking(
        seqs=seqs,
        idxs=idxs,
        chunk_size=chunk_size,
        func=func,
        pad_to_max=pad_to_max,
        z_array=z_array,
        lock_name=lock_name,
    )

    return True


def _compute_and_write_with_chunking(
    seqs: np.ndarray,
    idxs: list[int],
    chunk_size: int,
    func: Callable,
    pad_to_max: bool,
    z_array: zarr.Array,
    lock_name: str | None,
):
    """
    Compute blocks with internal chunking if needed and write to Zarr.
    """
    lock = Lock(lock_name) if lock_name is not None else None
    # Check if chunking is needed
    if len(seqs) > chunk_size:
        # Split into subchunks
        n_subchunks = max(1, math.ceil(len(seqs) / chunk_size))
        seq_chunks = np.array_split(seqs, n_subchunks)
        idx_chunks = np.array_split(idxs, n_subchunks)

        # Iterate over chunk pairs (upper triangle only)
        for i in range(len(seq_chunks)):
            for j in range(i, len(seq_chunks)):
                seqs_i = seq_chunks[i]
                seqs_j = seq_chunks[j]
                idxs_i = idx_chunks[i]
                idxs_j = idx_chunks[j]
                is_diagonal = i == j

                # Compute the block
                block = _compute_block_multicol(
                    seqs_i, seqs_j, func, pad_to_max
                )

                # Write to Zarr
                if lock is not None:
                    with lock:
                        z_array[np.ix_(idxs_i, idxs_j)] = block
                        if not is_diagonal:
                            z_array[np.ix_(idxs_j, idxs_i)] = block.T
                else:
                    z_array[np.ix_(idxs_i, idxs_j)] = block
                    if not is_diagonal:
                        z_array[np.ix_(idxs_j, idxs_i)] = block.T
    else:
        # No chunking needed
        block = _compute_block_multicol(seqs, seqs, func, pad_to_max)

        if lock is not None:
            with lock:
                z_array[np.ix_(idxs, idxs)] = block
        else:
            z_array[np.ix_(idxs, idxs)] = block


# def _process_chunk_pair(
#     member_chunk_i: np.ndarray,
#     member_chunk_j: np.ndarray,
#     dat_seq_clean: pd.DataFrame,
#     index_map: dict,
#     func: Callable,
#     pad_to_max: bool,
#     z_array: zarr.Array,
#     is_diagonal: bool,
#     lock: Lock | None,
# ):
#     """
#     Process a single pair of member chunks and write to Zarr.
#     """
#     # Extract sequences and indices
#     tmp_i = dat_seq_clean.loc[member_chunk_i]
#     tmp_j = dat_seq_clean.loc[member_chunk_j]
#     seqs_i = tmp_i.to_numpy(dtype=object)
#     seqs_j = tmp_j.to_numpy(dtype=object)
#     idxs_i = [index_map[m] for m in member_chunk_i]
#     idxs_j = [index_map[m] for m in member_chunk_j]

#     # Compute the block
#     block = _compute_block_multicol(seqs_i, seqs_j, func, pad_to_max)

#     # Write to Zarr
#     if lock is not None:
#         with lock:
#             z_array[np.ix_(idxs_i, idxs_j)] = block
#             if not is_diagonal:
#                 z_array[np.ix_(idxs_j, idxs_i)] = block.T
#     else:
#         z_array[np.ix_(idxs_i, idxs_j)] = block
#         if not is_diagonal:
#             z_array[np.ix_(idxs_j, idxs_i)] = block.T


def _auto_chunk_size(
    n: int,
    num_cores: int,
    memory_limit_gb: float | None = None,
    safety_fraction: float = 0.3,
) -> tuple[int, int]:
    """
    Compute dynamic chunk size to stay within memory budget.

    Parameters
    ----------
    n : int
        Total number of sequences
    num_cores : int
        Number of cores/workers
    memory_limit_gb : float, optional
        Memory limit per worker in GB
    safety_fraction : float
        Fraction of available memory to use (conservative for Dask overhead)

    Returns
    -------
    chunk_size : int
        Recommended chunk size
    n_chunks : int
        Number of chunks along one dimension
    """
    try:
        import psutil
    except ImportError:
        raise ImportError(
            "Please install psutil to enable automatic chunk size calculation: pip install psutil"
        )
    available_mem = psutil.virtual_memory().available / (1024**3)

    if memory_limit_gb is None:
        memory_limit_gb = available_mem * safety_fraction / num_cores

    mem_per_core = min(memory_limit_gb, available_mem / num_cores)

    # Each chunk block is chunk_size^2 * 8 bytes
    # Keep blocks small enough that 2-3 can fit in memory per worker
    chunk_size = int(math.sqrt((mem_per_core * (1024**3)) / 8 / 3))

    # Ensure minimum chunk size
    chunk_size = max(100, min(chunk_size, n))

    n_chunks = max(1, math.ceil(n / chunk_size))

    return chunk_size, n_chunks


def _setup_dask_client(
    num_cores: int,
    memory_limit_gb: float | None = None,
    verbose: bool = True,
) -> Client | None:
    """
    Setup Dask distributed client.

    Parameters
    ----------
    num_cores : int
        Number of workers
    memory_limit_gb : float, optional
        Memory limit per worker
    verbose : bool, optional
        Whether to logg.info client dashboard link.

    Returns
    -------
    Client or None
        Dask client if num_cores > 1, else None
    """
    try:
        from dask.distributed import Client
    except ImportError:
        raise ImportError(
            "Please install dask distributed to enable parallel processing: pip install dask distributed"
        )
    if num_cores <= 1:
        return None

    client_kwargs = {
        "n_workers": num_cores,
        "threads_per_worker": 1,
        "processes": True,  # Critical for memory isolation
    }

    if memory_limit_gb is not None:
        client_kwargs["memory_limit"] = f"{memory_limit_gb}GB"

    client = Client(**client_kwargs)

    if verbose:
        logg.info(f"Dask client started: {client.dashboard_link}")

    return client
