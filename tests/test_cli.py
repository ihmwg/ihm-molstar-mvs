from pathlib import Path
import pytest
import filecmp

ROOT = (Path(".") / "tests" / "test_data").resolve()

@pytest.mark.parametrize(
        "test_name,cmd,expected_output_file_name",

        [
            # Test1
            ("test1",
             "visualize_restraints https://pdb-ihm.org/cif/9a3v.cif -v",
             "9a3v.mvsj")

        ]
)
def test_cli(script_runner, tmp_path, monkeypatch,
             test_name, cmd, expected_output_file_name):

    monkeypatch.chdir(tmp_path)

    ret = script_runner.run(cmd.split())
    assert ret.success

    actual_output_file = tmp_path / expected_output_file_name
    expected_output_file = ROOT / f"{test_name}_expected_output"

    with open(actual_output_file, "r") as af:
        with open(expected_output_file, "r") as ef:

            for i, (actual, expected) in enumerate(zip(af.readlines(), ef.readlines())):
                if "timestamp" in expected:
                    continue

                if "version" in expected:
                    continue

                assert actual == expected, f"line {i} of expected doesnt match with actual: {expected}, {actual}"


    expected_stdout_file = ROOT / f"{test_name}_expected_stdout"
    with open(expected_stdout_file, "r") as f:
        expected_stdout = f.read()

    assert ret.stdout == expected_stdout



