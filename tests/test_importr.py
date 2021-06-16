#!/usr/bin/env python

def test_importrpy2():
    import rpy2
    from rpy2.robjects.packages import importr
    sh = importr('shazam')
    assert sh.__module__ == 'rpy2.robjects.packages'
