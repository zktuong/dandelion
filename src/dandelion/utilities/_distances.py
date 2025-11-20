from __future__ import annotations

from Bio.Align import substitution_matrices as _submat
from polyleven import levenshtein
from typing import Callable, Protocol, runtime_checkable


def dist_func_long_sep(
    x: list[str],
    y: list[str],
    metric: Metric,
    pad_to_max: bool = False,
    sep: str = "#",
) -> float:
    """
    Concatenate column-wise with long separators and apply the metric.
    """
    # If the metric is a substitution matrix, skip sep padding
    if isinstance(metric, SubstitutionMatrixMetric):
        # Just concatenate without adding separators
        s1 = "".join(x)
        s2 = "".join(y)
    elif pad_to_max:
        max_len = [max(len(a), len(b)) for a, b in zip(x, y)]
        s1_parts, s2_parts = [], []

        for s1, s2, le in zip(x, y, max_len):
            # Pad each element
            s1_parts.append(s1.ljust(le + 1, sep))
            s2_parts.append(s2.ljust(le + 1, sep))

        # Join with per-column separators
        s1_result, s2_result = [], []
        for i, (s1_part, s2_part, le) in enumerate(
            zip(s1_parts, s2_parts, max_len)
        ):
            s1_result.append(s1_part)
            s2_result.append(s2_part)

            # Add separator between columns (but not after the last one)
            if i < len(max_len) - 1:
                col_sep = sep * (le + 2)  # Longer than padded length (le + 1)
                s1_result.append(col_sep)
                s2_result.append(col_sep)

        s1 = "".join(s1_result)
        s2 = "".join(s2_result)
    else:
        # Dynamically choose separator length: longer than the max column
        try:
            max_x = max(len(s) for s in x)
        except ValueError:
            max_x = 0
        try:
            max_y = max(len(s) for s in y)
        except ValueError:
            max_y = 0
        max_len = max(max_x, max_y)
        long_sep = sep * (max_len + 1)
        s1 = long_sep.join(x)
        s2 = long_sep.join(y)

    return metric.compute(s1, s2)


@runtime_checkable
class Metric(Protocol):
    """
    Protocol for metric objects used to compute distance between two strings.

    Implementors must provide `compute(s1: str, s2: str) -> float`.
    """

    def compute(self, s1: str, s2: str) -> float: ...


class CallableMetric:
    """Wraps a callable f(s1, s2) into a Metric object."""

    def __init__(self, func: Callable[[str, str], float]):
        self.func = func

    def compute(self, s1: str, s2: str) -> float:
        return float(self.func(s1, s2))


class LevenshteinMetric:
    """Metric that computes edit distance using the polyleven implementation."""

    def compute(self, s1: str, s2: str) -> float:
        return float(levenshtein(s1, s2))


class SubstitutionMatrixMetric:
    """
    Metric that uses a substitution matrix (e.g. BLOSUM62).

    Distance is computed as:
        distance = max_possible_score - actual_score
    Optionally normalized to [0, 1].

    Parameters
    ----------
    matrix_name : str
        Name of the Biopython substitution matrix (e.g., "BLOSUM62").
    gap_penalty : float, default=-4.0
        Penalty for gaps or characters not in the matrix.
    """

    def __init__(
        self,
        matrix_name: str,
        gap_penalty: float = -4.0,
    ):
        self.gap_penalty = float(gap_penalty)

        # Load Biopython substitution matrix
        try:
            bio_mat = _submat.load(matrix_name)
        except Exception as e:
            raise ValueError(
                f"Could not load substitution matrix '{matrix_name}'. "
                "Please ensure it is a valid matrix name supported by Biopython."
            ) from e

        # Convert to dict for fast lookup
        self.matrix: dict[tuple[str, str], float] = {}
        for a in bio_mat.alphabet:
            for b in bio_mat.alphabet:
                v = float(bio_mat[a, b])
                self.matrix[(a, b)] = v
                self.matrix[(b, a)] = v  # symmetric access

    def compute(self, s1: str, s2: str) -> float:
        """Compute distance between two sequences based on substitution matrix."""

        min_len = min(len(s1), len(s2))
        score = 0.0

        # Sum substitution matrix scores for aligned positions
        for i in range(min_len):
            a, b = s1[i], s2[i]
            score += self.matrix.get((a, b), self.gap_penalty)

        # Apply gap penalties for unequal lengths
        len_diff = abs(len(s1) - len(s2))
        score += len_diff * self.gap_penalty

        # Compute max possible score for s1 (perfect self-alignment)
        max_score = sum(self.matrix.get((c, c), 0.0) for c in s1)
        max_score += len_diff * self.gap_penalty  # consider gaps

        # Distance = max_score - actual_score, clamped at 0
        distance = max(0.0, max_score - score)

        return distance


# -------------------------
# Metric resolver
# -------------------------
def resolve_metric(
    dist_func: Callable[[str, str], float] | str | Metric | None,
) -> Metric:
    """
    Convert user-supplied dist_func into a Metric object.

    Parameters
    ----------
    dist_func : Callable[[str, str], float] | str | Metric | None
        Distance function specification.
    Returns
    -------
    Metric
        Resolved Metric object.
    """
    if dist_func is None:
        return LevenshteinMetric()

    if isinstance(dist_func, Metric):
        return dist_func

    if callable(dist_func):
        return CallableMetric(dist_func)

    if isinstance(dist_func, str):
        return SubstitutionMatrixMetric(dist_func)

    raise TypeError(
        "dist_func must be None, callable, Metric, or substitution matrix name string."
    )
