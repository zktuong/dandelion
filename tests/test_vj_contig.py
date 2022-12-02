#!/usr/bin/env python
import pytest
from dandelion.preprocessing._preprocessing import check_productive_vj


def test_vj_contig():
    test = {"test_1": 3, "test_2": 3, "test_3": 3}
    keep, extra, ambiguous = check_productive_vj(test)
    assert len(keep) == 3
    assert len(extra) == 0
    assert len(ambiguous) == 0
