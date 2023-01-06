#!/usr/bin/env python
"""shannon module."""
# Lifted from skibio==0.5.6
# because of issue with having skbio as a dependency

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import numpy as np
from dandelion.utilities._utilities import _validate_counts_vector


def shannon(counts, base=2):
    r"""Calculate Shannon entropy of counts, default in bits.
    Shannon-Wiener diversity index is defined as:
    .. math::
       H = -\sum_{i=1}^s\left(p_i\log_2 p_i\right)
    where :math:`s` is the number of OTUs and :math:`p_i` is the proportion of
    the community represented by OTU :math:`i`.
    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    base : scalar, optional
        Logarithm base to use in the calculations.
    Returns
    -------
    double
        Shannon diversity index H.
    Notes
    -----
    The implementation here is based on the description given in the SDR-IV
    online manual [1]_ except that the default logarithm base used here is 2
    instead of :math:`e`.
    References
    ----------
    .. [1] http://www.pisces-conservation.com/sdrhelp/index.html
    """
    counts = _validate_counts_vector(counts)
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)
