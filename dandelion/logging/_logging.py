#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2021-02-06 13:18:58
# @Last Modified by:   Kelvin
# @Last Modified time: 2022-11-21 17:17:10
"""logging module."""
from typing import List

modules = [
    "dandelion",
    "pandas",
    "numpy",
    "matplotlib",
    "networkx",
    "scipy",
    "distance",
    "polyleven",
]


# borrowed from scanpy's logging module
def _versions_dependencies(dependencies: List[str]):
    """Version dependencies."""
    for mod in dependencies:
        mod_name, dist_name = mod if isinstance(mod, tuple) else (mod, mod)
        try:
            imp = __import__(mod_name)
            yield dist_name, imp.__version__
        except (ImportError, AttributeError):
            pass


def print_versions(dependencies: List[str] = modules):
    """
    Versions that are essential for dandelion's operation.
    """
    print(
        " ".join(
            f"{mod}=={ver}" for mod, ver in _versions_dependencies(dependencies)
        )
    )


def print_header(dependencies: List[str] = modules):
    """
    Versions that are essential for dandelion's operation.
    """
    print(
        " ".join(
            f"{mod}=={ver}" for mod, ver in _versions_dependencies(dependencies)
        )
    )
