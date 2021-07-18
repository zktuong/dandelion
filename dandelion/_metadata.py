#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2021-04-03 16:46:13
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-07-18 17:14:24

# try:
# 	from .version import __version__
# except ImportError:
# 	__version__ = ''

from setuptools_scm import get_version

__version__ = get_version(root="..", relative_to=__file__).split('+')[0]
__author__ = "Zewen Kelvin Tuong"
__email__ = "kt16@sanger.ac.uk"
__url__ = "https://www.github.com/zktuong/dandelion/"
__docs__ = "https://sc-dandelion.readthedocs.io/"
__classifiers__ = [
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: R",
    "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
]
