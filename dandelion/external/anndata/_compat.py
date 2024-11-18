import numpy as np
import pandas as pd

from scipy.sparse import spmatrix, sparray
from typing import Sequence

# from https://github.com/scverse/anndata/blob/53537b5219ff82cbdee96b7733172fb114e80ca8/src/anndata/compat/__init__.py#L48-L49
Index1D = slice | int | str | np.int64 | np.ndarray
Index = Index1D | tuple[Index1D, Index1D] | spmatrix | sparray


# from https://github.com/scverse/anndata/blob/53537b5219ff82cbdee96b7733172fb114e80ca8/src/anndata/_core/index.py#L43-L111
def _normalize_index(
    indexer: (
        slice
        | np.integer
        | int
        | str
        | Sequence[bool | int | np.integer]
        | np.ndarray
        | pd.Index
    ),
    index: pd.Index,
) -> slice | int | np.ndarray:  # ndarray of int or bool
    """normalize index from anndata."""
    if not isinstance(index, pd.RangeIndex):
        msg = "Don’t call _normalize_index with non-categorical/string names"
        assert index.dtype != float, msg
        assert index.dtype != int, msg

    # the following is insanely slow for sequences,
    # we replaced it using pandas below
    def name_idx(i):
        """name index."""
        if isinstance(i, str):
            i = index.get_loc(i)
        return i

    if isinstance(indexer, slice):
        start = name_idx(indexer.start)
        stop = name_idx(indexer.stop)
        # string slices can only be inclusive, so +1 in that case
        if isinstance(indexer.stop, str):
            stop = None if stop is None else stop + 1
        step = indexer.step
        return slice(start, stop, step)
    elif isinstance(indexer, (np.integer, int)):
        return indexer
    elif isinstance(indexer, str):
        return index.get_loc(indexer)  # int
    elif isinstance(
        indexer, (Sequence, np.ndarray, pd.Index, spmatrix, np.matrix, sparray)
    ):
        if hasattr(indexer, "shape") and (
            (indexer.shape == (index.shape[0], 1))
            or (indexer.shape == (1, index.shape[0]))
        ):
            if isinstance(indexer, (spmatrix, sparray)):
                indexer = indexer.toarray()
            indexer = np.ravel(indexer)
        if not isinstance(indexer, (np.ndarray, pd.Index)):
            indexer = np.array(indexer)
            if len(indexer) == 0:
                indexer = indexer.astype(int)
        if issubclass(indexer.dtype.type, (np.integer, np.floating)):
            return indexer  # Might not work for range indexes
        elif issubclass(indexer.dtype.type, np.bool_):
            if indexer.shape != index.shape:
                raise IndexError(
                    f"Boolean index does not match AnnData’s shape along this "
                    f"dimension. Boolean index has shape {indexer.shape} while "
                    f"AnnData index has shape {index.shape}."
                )
            return indexer
        else:  # indexer should be string array
            positions = index.get_indexer(indexer)
            if np.any(positions < 0):
                not_found = indexer[positions < 0]
                raise KeyError(
                    f"Values {list(not_found)}, from {list(indexer)}, "
                    "are not valid obs/ var names or indices."
                )
            return positions  # np.ndarray[int]
    else:
        raise IndexError(f"Unknown indexer {indexer!r} of type {type(indexer)}")


def unpack_index(index: Index) -> tuple[Index1D, Index1D]:
    """unpack index from anndata."""
    if not isinstance(index, tuple):
        return index, slice(None)
    elif len(index) == 2:
        return index
    elif len(index) == 1:
        return index[0], slice(None)
    else:
        raise IndexError("invalid number of indices")
