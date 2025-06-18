# Lifted from skibio==0.5.6
# because of issue with having skbio as a dependency

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from dandelion.external.skbio._utils import validate_counts_vector


def osd(counts):
    """Calculate observed OTUs, singles, and doubles.
    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    Returns
    -------
    osd : tuple
        Observed OTUs, singles, and doubles.
    See Also
    --------
    observed_otus
    singles
    doubles
    Notes
    -----
    This is a convenience function used by many of the other measures that rely
    on these three measures.
    """
    counts = validate_counts_vector(counts)
    return observed_otus(counts), singles(counts), doubles(counts)


def singles(counts):
    """Calculate number of single occurrences (singletons).
    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    Returns
    -------
    int
        Singleton count.
    """
    counts = validate_counts_vector(counts)
    return (counts == 1).sum()


def doubles(counts):
    """Calculate number of double occurrences (doubletons).
    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    Returns
    -------
    int
        Doubleton count.
    """
    counts = validate_counts_vector(counts)
    return (counts == 2).sum()


def observed_otus(counts):
    """Calculate the number of distinct OTUs.
    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    Returns
    -------
    int
        Distinct OTU count.
    """
    counts = validate_counts_vector(counts)
    return (counts != 0).sum()


def chao1(counts, bias_corrected=True):
    r"""Calculate chao1 richness estimator.
    Uses the bias-corrected version unless `bias_corrected` is ``False`` *and*
    there are both singletons and doubletons.
    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    bias_corrected : bool, optional
        Indicates whether or not to use the bias-corrected version of the
        equation. If ``False`` *and* there are both singletons and doubletons,
        the uncorrected version will be used. The biased-corrected version will
        be used otherwise.
    Returns
    -------
    double
        Computed chao1 richness estimator.
    See Also
    --------
    chao1_ci
    Notes
    -----
    The uncorrected version is based on Equation 6 in [1]_:
    .. math::
       chao1=S_{obs}+\frac{F_1^2}{2F_2}
    where :math:`F_1` and :math:`F_2` are the count of singletons and
    doubletons, respectively.
    The bias-corrected version is defined as
    .. math::
       chao1=S_{obs}+\frac{F_1(F_1-1)}{2(F_2+1)}
    References
    ----------
    .. [1] Chao, A. 1984. Non-parametric estimation of the number of classes in
       a population. Scandinavian Journal of Statistics 11, 265-270.
    """
    counts = validate_counts_vector(counts)
    o, s, d = osd(counts)

    if not bias_corrected and s and d:
        return o + s**2 / (d * 2)
    else:
        return o + s * (s - 1) / (2 * (d + 1))
