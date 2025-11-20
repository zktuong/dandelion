from __future__ import annotations
import math
import zarr
import zarr.storage

import dask.array as da
import numpy as np
import pandas as pd

from dask import compute
from dask.diagnostics import ProgressBar
from dask.distributed import Client, progress, Lock
from pathlib import Path
from scanpy import logging as logg
from zarr.codecs import BloscCodec

from dandelion.utilities._distances import dist_func_long_sep, Metric


def calculate_distance_matrix_zarr(
    dat_seq: pd.DataFrame,
    metric: Metric,
    pad_to_max: bool = False,
    membership: dict | None = None,
    zarr_path: Path | str | None = None,
    chunk_size: int | None = None,
    max_clones_per_chunk: int | None = None,
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
    metric : Metric
        Distance metric to use.
    pad_to_max : bool, optional
        Whether to pad sequences to maximum length before distance calculation.
    zarr_path : str, optional
        Path to save Zarr array (directory). If None, uses 'distance_matrix.zarr'.
    chunk_size : int | None, optional
        Number of sequences per computation chunk. If None, auto-computed.
    max_clones_per_chunk : int | None, optional
        Maximum number of clones per computation chunk. If None, no limit is applied
        and instead if will be based on memory limits.
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
    # reduce to only clonotypes with >1 member
    if membership is not None:
        membership = {k: v for k, v in membership.items() if len(v) > 1}
        if len(membership) < 1:
            # there's nothing to compute. so just return
            logg.info(
                msg=(
                    "No clonotypes with more than one member found. "
                    "Try using distance_mode = 'full' instead.\n"
                    "Skipping distance computation.\n"
                )
            )
            return
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

    zarr_lock = Lock("zarr_distance_matrix_lock") if num_cores > 1 else None

    # Store metadata
    z_array.attrs["n_sequences"] = n
    z_array.attrs["chunk_size"] = chunk_size
    z_array.attrs["index"] = list(dat_seq.index)
    z_array.attrs["columns"] = list(dat_seq.columns)

    if verbose:
        logg.info(f"\nCreated Zarr array at: {zarr_path}")
        logg.info(f"Compressor: {comp}\n")

    # Setup Dask client
    client = _setup_dask_client(num_cores, memory_limit_gb)

    try:
        # Compute distances and write blocks as they complete
        _compute_multicol_distances_streaming(
            dat_seq=dat_seq,
            metric=metric,
            pad_to_max=pad_to_max,
            z_array=z_array,
            chunk_size=chunk_size,
            num_cores=num_cores,
            chunk_batch_limit=max_clones_per_chunk,
            client=client,
            membership=membership,
            lock=zarr_lock,
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
    metric: Metric,
    pad_to_max: bool,
    z_array: zarr.Array,
    chunk_size: int,
    num_cores: int,
    chunk_batch_limit: int | None = None,
    client: Client | None = None,
    membership: dict | None = None,
    lock: Lock | None = None,
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

    seqs_np = dat_seq_clean.to_numpy(dtype=object)
    seqs = da.from_delayed(
        dask.delayed(lambda x: x)(seqs_np), shape=seqs_np.shape, dtype=object
    )
    m = len(seqs)
    if m <= 1:
        return None

    delayed_blocks = []
    logg.info(
        msg="Setting up Dask scheduler and task graph...",
    )
    if membership is not None:
        # Build index map
        index_list = list(dat_seq_clean.index)
        index_map = {idx: i for i, idx in enumerate(index_list)}
        # Build reordering (flatten membership)
        group_idx_flat = np.array(
            [index_map[cell] for group in membership.values() for cell in group]
        )

        group_lengths = np.fromiter(
            (len(v) for v in membership.values()), dtype=int
        )
        boundaries = np.cumsum(group_lengths[:-1])

        # Reorder seqs contiguously - convert result to dask array
        seqs_reordered_np = seqs_np[group_idx_flat]
        # Now each membership group is a contiguous slice (a view)
        seqs_m = np.split(seqs_reordered_np, boundaries)
        idx_m = np.split(group_idx_flat, boundaries)

        # seqs_reordered = [da.from_array(s, chunks=s.shape) for s in seqs_m]
        # group_idx_flat_da = [da.from_array(i, chunks=i.shape) for i in idx_m]
        # if client is not None:
        #     seqs_reordered = client.persist(seqs_reordered)
        #     group_idx_flat_da = client.persist(group_idx_flat_da)

        batched_idx_m = []
        batched_seqs_m = []

        current_batch_idx = []
        current_batch_seqs = []
        current_batch_size = 0  # total number of sequences in current batch
        current_batch_clones = 0  # total number of clones in current batch

        # Calculate soft batch limit if not provided
        if chunk_batch_limit is None and num_cores > 1:
            total_clonotypes = len(membership)
            # Target: distribute clonotypes across cores with some overhead for load balancing
            # Use 2x cores to allow for better parallelization
            target_batches = num_cores * 2
            soft_chunk_batch_limit = max(1, total_clonotypes // target_batches)
        else:
            soft_chunk_batch_limit = None

        for idx, seq in zip(idx_m, seqs_m):
            n = len(idx)
            # Determine if the current batch should be finalized
            finalize = False
            if (
                chunk_batch_limit is not None
                and current_batch_clones >= chunk_batch_limit
            ):
                finalize = True
            elif (
                soft_chunk_batch_limit is not None
                and current_batch_clones >= soft_chunk_batch_limit
            ):
                # Use soft limit if no hard limit and batch is getting large
                finalize = True
            elif (
                chunk_batch_limit is None
                and current_batch_size + n > chunk_size
                and current_batch_size > 0
            ):
                finalize = True

            if finalize:
                # finalize current batch
                batched_idx_m.append(current_batch_idx)
                batched_seqs_m.append(current_batch_seqs)
                current_batch_idx, current_batch_seqs = [], []
                current_batch_size, current_batch_clones = 0, 0

            # Add current clonotype to batch
            current_batch_idx.append(idx)
            current_batch_seqs.append(seq)
            current_batch_size += n
            current_batch_clones += 1

        # Add the last batch
        if current_batch_idx:
            batched_idx_m.append(current_batch_idx)
            batched_seqs_m.append(current_batch_seqs)

        # average number of clonotypes per batch
        avg_clonotypes_per_batch = sum(len(b) for b in batched_seqs_m) / len(
            batched_idx_m
        )
        # # lock per group
        # membership_locks = (
        #     [Lock(f"membership-{i}") for i in range(len(batched_idx_m))]
        #     if num_cores > 1
        #     else None
        # )
        logg.info(
            msg=f"Created {len(batched_seqs_m)} chunks for distance computation of ~{math.ceil(avg_clonotypes_per_batch)} clonotypes per batch...",
        )
        for idx, seq in zip(batched_idx_m, batched_seqs_m):
            delayed_blocks.append(
                dask.delayed(_process_batch)(
                    batch_seqs_m=seq,
                    batch_idx_m=idx,
                    metric=metric,
                    pad_to_max=pad_to_max,
                    z_array=z_array,
                    lock=lock,
                )
            )
    else:
        # Determine number of chunks
        n_chunks = max(1, math.ceil(m / chunk_size))
        if n_chunks < num_cores:
            n_chunks = num_cores
        # Work with numpy for splitting, then convert chunks to dask arrays
        seqs_computed = seqs.compute() if hasattr(seqs, "compute") else seqs_np
        chunks_list = np.array_split(seqs_computed, n_chunks)
        # Convert each chunk to a dask array if client is provided
        if client is not None:
            chunks_list = [
                client.persist(
                    da.from_delayed(
                        dask.delayed(lambda x: x)(chunk),
                        shape=chunk.shape,
                        dtype=object,
                    )
                )
                for chunk in chunks_list
            ]

        chunk_sizes = [len(c) for c in chunks_list]
        cum_sizes = np.cumsum([0] + chunk_sizes)
        for i in range(len(chunks_list)):
            for j in range(i, len(chunks_list)):
                if chunk_sizes[i] == 0 or chunk_sizes[j] == 0:
                    continue
                start_i, end_i = cum_sizes[i], cum_sizes[i + 1]
                start_j, end_j = cum_sizes[j], cum_sizes[j + 1]
                is_diagonal = i == j
                delayed_blocks.append(
                    dask.delayed(_compute_and_write_block)(
                        chunk_i=chunks_list[i],
                        chunk_j=chunks_list[j],
                        metric=metric,
                        pad_to_max=pad_to_max,
                        z_array=z_array,
                        start_i=start_i,
                        end_i=end_i,
                        start_j=start_j,
                        end_j=end_j,
                        is_diagonal=is_diagonal,
                        lock=lock,
                    )
                )

    logg.info(
        f"Starting computation of {len(delayed_blocks)} chunks...",
    )
    # Compute blocks - they write to Zarr as they complete
    if client is not None:
        futures = client.compute(delayed_blocks)
        progress(futures)
        client.gather(futures)
    else:
        with ProgressBar():
            compute(*delayed_blocks, scheduler="threads")

    return True


def _compute_block_multicol(
    seqs_i: np.ndarray,
    seqs_j: np.ndarray,
    metric: Metric,
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
                row_i, row_j, metric=metric, pad_to_max=pad_to_max
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
    metric: Metric,
    pad_to_max: bool,
    z_array: zarr.Array,
    start_i: int,
    end_i: int,
    start_j: int,
    end_j: int,
    is_diagonal: bool,
    lock: Lock | None,
):
    """
    Compute a single block and write it directly to Zarr.
    """
    block = _compute_block_multicol(chunk_i, chunk_j, metric, pad_to_max)

    if lock is not None:
        with lock:
            # Write to Zarr immediately
            z_array[start_i:end_i, start_j:end_j] = block

            # If not diagonal block, write transpose too
            if not is_diagonal:
                z_array[start_j:end_j, start_i:end_i] = block.T
    else:
        # Write to Zarr immediately
        z_array[start_i:end_i, start_j:end_j] = block

        # If not diagonal block, write transpose too
        if not is_diagonal:
            z_array[start_j:end_j, start_i:end_i] = block.T
    return True  # Just return success flag


def _process_batch(
    batch_seqs_m: list[np.ndarray],
    batch_idx_m: list[list[int]],
    metric: Metric,
    pad_to_max: bool,
    z_array: zarr.Array,
    lock: Lock | None,  # Single lock for this batch
):
    seqs_cat = np.concatenate(batch_seqs_m, axis=0)
    idxs_cat = np.concatenate(batch_idx_m, axis=0)

    lengths = [len(s) for s in batch_seqs_m]
    boundaries = np.cumsum(lengths)
    boundaries = np.insert(boundaries, 0, 0)

    _compute_and_write_membership_cat(
        seqs_cat=seqs_cat,
        idxs_cat=idxs_cat,
        boundaries=boundaries,
        metric=metric,
        pad_to_max=pad_to_max,
        z_array=z_array,
        lock=lock,
    )


def _compute_and_write_membership_cat(
    seqs_cat: np.ndarray,
    idxs_cat: np.ndarray,
    boundaries: np.ndarray,
    metric: Metric,
    pad_to_max: bool,
    z_array: zarr.Array,
    lock: list[Lock] | None,
):
    """
    Compute distances within each membership group only.
    """
    # Compute full block once
    block = _compute_block_multicol(seqs_cat, seqs_cat, metric, pad_to_max)

    n_groups = len(boundaries) - 1

    # Only process diagonal blocks (within each group)
    for g in range(n_groups):
        s, e = boundaries[g], boundaries[g + 1]
        idx = idxs_cat[s:e]

        # Extract subblock for this group only
        subblock = block[s:e, s:e]

        # Write to zarr
        if lock is not None:
            with lock:
                z_array[np.ix_(idx, idx)] = subblock
        else:
            z_array[np.ix_(idx, idx)] = subblock


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
    num_cores: int, memory_limit_gb: float | None = None
) -> Client | None:
    """
    Setup Dask distributed client.

    Parameters
    ----------
    num_cores : int
        Number of workers
    memory_limit_gb : float, optional
        Memory limit per worker

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
        "threads_per_worker": 1,  # for simplicity and to avoid GIL issues
        "processes": True,  # Critical for memory isolation
    }

    if memory_limit_gb is not None:
        client_kwargs["memory_limit"] = f"{memory_limit_gb}GB"

    client = Client(**client_kwargs)

    logg.info(f"Dask client started: {client.dashboard_link}")

    return client
