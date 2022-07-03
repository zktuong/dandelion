#!/usr/bin/env python
"""test tools clones."""
import dandelion as ddl
import pytest

from pathlib import Path


@pytest.mark.usefixtures("airr_generic")
def test_find_clones_other_options(airr_generic):
    """Test find clones."""
    vdj = ddl.pp.check_contigs(airr_generic, productive_only=False)
    with pytest.raises(ValueError):
        vdj = ddl.tl.find_clones(vdj, recalculate_length=False)
    vdj.data["junction_aa_length"] = 10
    ddl.tl.find_clones(vdj, recalculate_length=False)
    assert not vdj.data.clone_id.empty
    assert not vdj.metadata.clone_id.empty


@pytest.mark.usefixtures("create_testfolder", "airr_generic")
def test_find_clones_file(create_testfolder, airr_generic):
    """Test find clones from file."""
    in_file = str(create_testfolder) + "/test_airr.tsv"
    out_file = str(create_testfolder) + "/test_airr_clone.tsv"
    airr_generic.to_csv(in_file, sep="\t", index=False)
    ddl.tl.find_clones(in_file)
    assert Path(out_file) in list(create_testfolder.iterdir())
