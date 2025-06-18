import numpy as np


## from skbio==0.5.6
def validate_counts_vector(
    counts: np.array, suppress_cast: bool = False
) -> np.array:
    """Validate and convert input to an acceptable counts vector type.
    Note: may not always return a copy of `counts`!
    """
    counts = np.asarray(counts)
    if not np.all(np.isreal(counts)):
        raise ValueError("Counts vector must contain real-valued entries.")
    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")
    elif (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts
