#!/opt/conda/envs/sc-dandelion-container/bin/python
"""dandelion container test script"""
import os
import pandas as pd
import dandelion as ddl
from subprocess import run


def test_callscript():
    """Test script to run preprocessing."""
    p = run(
        ["python", "/share/dandelion_preprocess.py", "-h"],
        capture_output=True,
        encoding="utf8",
    )
    assert p.returncode == 0
    assert p.stdout != ""


def test_container():
    """Test script to run container."""
    os.system(
        "cd /tests; python /share/dandelion_preprocess.py --meta test.csv;"
    )
    dat = pd.read_csv(
        "/tests/sample_test_10x/dandelion/filtered_contig_dandelion.tsv",
        sep="\t",
    )
    assert not dat["c_call"].empty
    assert not dat["v_call_genotyped"].empty
    assert not dat["mu_count"].empty
    assert not dat["mu_freq"].empty
    vdj = None
    try:
        vdj = ddl.Dandelion(dat)
    except:
        pass
    assert vdj is not None
