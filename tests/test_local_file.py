import pytest
from pathlib import Path
import requests

from ihm_vis.utils.local_file import LocalFile

TEST_ROOT = (Path(".") / "tests").resolve()

def test_serve_local_file():

    test_cif_file = str(TEST_ROOT / "test_data" / "9a3v.cif")
    lf = LocalFile(test_cif_file)

    # Should only be present when using the serve context manager
    with pytest.raises(requests.ConnectionError) as e:
        before = requests.get(lf.url)

    # Should find file
    with lf.serve():
        during = requests.get(lf.url)

    # Should not find file now that context closed
    with pytest.raises(requests.ConnectionError) as e:
        after = requests.get(lf.url)

    assert during.status_code == 200, "Should find file during"

    with open(test_cif_file, "r") as f:
        expected = f.read()

    assert during.text == expected

