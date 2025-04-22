import pytest
from ihm_vis import IHM_Builder
from pathlib import Path
import pandas as pd
import tempfile
import os

##################################
# Tests for file parsing
##################################

VALID_LOCAL_FILES = ["test_files/8zz1.cif", "test_data/8zzs.cif"]
INVALID_LOCAL_FILES = ["test_files/invalid.cif", "cif.txt", "test_data/8zz1"]

VALID_URLS = ["https://pdb-ihm.org/cif/9a3v.cif", "https://files.rcsb.org/view/9EV8.cif", "https://files.rcsb.org/view/8ZZ1.cif"]
INVALID_URLS = ["https://pdb-ihm.org/cif/9a3v", "https://files.rcsb.org/view/9EV8", "http://invalid-url.org/file.cif"]

@pytest.mark.parametrize("file_path", VALID_LOCAL_FILES)
def test_valid_local_file(file_path):
    builder = IHM_Builder(source=file_path)
    assert builder.source_type == "file"
    assert builder.cif is not None

@pytest.mark.parametrize("file_path", INVALID_LOCAL_FILES)
def test_invalid_local_file(file_path):
    with pytest.raises(ValueError):
        IHM_Builder(source=file_path)

@pytest.mark.parametrize("url", VALID_URLS)
def test_valid_url(url):
    builder = IHM_Builder(source=url)
    assert builder.source_type == "url"
    assert builder.cif is not None

@pytest.mark.parametrize("url", INVALID_URLS)
def test_invalid_url(url):
    with pytest.raises(ValueError):
        IHM_Builder(source=url)

@pytest.mark.parametrize("source", VALID_LOCAL_FILES + VALID_URLS)
def test_all_valid_sources(source):
    """Test for autodetecting of source"""
    builder = IHM_Builder(source=source)
    assert builder.cif is not None

@pytest.mark.paramtrize("url", VALID_URLS)
def test_null_coords():
    builder = IHM_Builder(source=url)
    df = builder.get_cross_links()
    assert not df['atom_id_1_coords'].isnull().any(), "atom_id_1_coords has null values"
    assert not df['atom_id_2_coords'].isnull().any(), "atom_id_2_coords has null values"

# Add test for empty dataframe when no restraints to visualize