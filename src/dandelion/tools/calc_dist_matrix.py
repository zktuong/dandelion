import numpy as np
from scipy.spatial.distance import pdist, squareform
from polyleven import levenshtein
from tqdm import tqdm
from dask import delayed, compute
import dask


def compute_chunk_distances(chunk):
    """Compute pairwise distances for a chunk of sequences."""
    return squareform(
        pdist(chunk.reshape(-1, 1), lambda x, y: levenshtein(x[0], y[0]))
    )


def parallel_levenshtein_matrix(sequences, n_chunks=400):
    """
    Compute Levenshtein distance matrix in parallel using Dask.

    Args:
        sequences: Array of sequences to compare
        n_chunks: Number of chunks for parallel processing

    Returns:
        Numpy array of pairwise distances
    """
    dask.config.set(scheduler="threads")

    # Split data into chunks
    chunks = np.array_split(sequences, n_chunks)

    # Create delayed computations
    delayed_results = [
        delayed(compute_chunk_distances)(chunk) for chunk in chunks
    ]

    # Execute computations in parallel
    results = compute(*delayed_results)

    # Combine results into full matrix
    n = len(sequences)
    dist_matrix = np.zeros((n, n))

    start_idx = 0
    for res in results:
        end_idx = start_idx + res.shape[0]
        dist_matrix[start_idx:end_idx, start_idx:end_idx] = res
        start_idx = end_idx

    return dist_matrix


def calculate_distance_matrix(dat_seq):
    """
    Calculate combined distance matrix for all sequence columns.

    Args:
        dat_seq: DataFrame containing sequence columns, (passed from the network function)

    Returns:
        Combined distance matrix with diagonal set to NaN
    """
    distance_matrices = {}

    for col in tqdm(
        dat_seq.columns,
        desc="Calculating distance matrices",
        bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
    ):
        # Clean sequences and convert to array
        clean_sequences = dat_seq[col].str.replace("[.]", "", regex=True)
        seq_array = np.array(clean_sequences.dropna().tolist())

        # Compute distance matrix
        distance_matrices[col] = parallel_levenshtein_matrix(seq_array)

    # Combine all distance matrices
    valid_matrices = [
        m for m in distance_matrices.values() if isinstance(m, np.ndarray)
    ]
    total_dist = np.sum(valid_matrices, axis=0)

    # Set diagonal to NaN
    np.fill_diagonal(total_dist, np.nan)

    return total_dist
