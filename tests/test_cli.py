from pathlib import Path
import shutil
import pytest
from ihm_vis import style

# Some tests will modify the user's style
# force a reload of the module to 
# ensure previous tests dont affect
# others
@pytest.fixture(scope="function", autouse=True)
def reset_user_style():
    style.reset_style()
    yield


TEST_ROOT = (Path(".") / "tests").resolve()

@pytest.mark.parametrize(
        "test_name,cmd,expected_output_file_name,files_to_copy",

        [
            # Test1 - simple example
            ("test1",
             "visualize_restraints https://pdb-ihm.org/cif/9a3v.cif -v",
             "9a3v.mvsj",
             []),


            # Test2 - all-in-one example from Readme
            ("test2",
             "visualize_restraints https://pdb-ihm.org/cif/9a3v.cif -v -f across_chains -s violated_and_compliant -c my_style.yaml",
             "9a3v.mvsj",
             ["my_style.yaml"]),


            # Test3 - local file
            ("test3",
             "visualize_restraints 9a3v.cif -v -t fromlocal -o myoutput",
             "myoutput.mvsj",
             ["9a3v.cif"]),

        ]
)
def test_cli(script_runner, tmp_path, monkeypatch,
             test_name, cmd, expected_output_file_name, files_to_copy):

    monkeypatch.chdir(tmp_path)

    for f in files_to_copy:
        shutil.copy(TEST_ROOT / "test_data" / f, tmp_path / f)

    ret = script_runner.run(cmd.split())
    assert ret.success

    actual_output_file = tmp_path / expected_output_file_name
    expected_output_file = TEST_ROOT / "test_data" / f"{test_name}_expected_output"

    with open(actual_output_file, "r") as af:
        with open(expected_output_file, "r") as ef:

            for i, (actual, expected) in enumerate(zip(af.readlines(), ef.readlines())):
                if "timestamp" in expected:
                    continue

                if "version" in expected:
                    continue

                assert actual == expected, f"line {i} of {actual_output_file.resolve()} vs {expected_output_file.resolve()} doesnt match: {expected}, {actual}"


    expected_stdout_file = TEST_ROOT / "test_data" / f"{test_name}_expected_stdout"
    with open(expected_stdout_file, "r") as f:
        expected_stdout = f.read()

    assert ret.stdout == expected_stdout



