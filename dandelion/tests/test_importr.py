#!/usr/bin/env python


def test_importrpy2():

    from rpy2.robjects.packages import importr
    from rpy2.rinterface import NULL
    from rpy2.robjects import pandas2ri
    sh = importr('shazam')

    assert sh.__module__ == 'rpy2.robjects.packages'
