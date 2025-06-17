import os
import sys
import subprocess
import re
import pytest

OPENSN = os.environ.get("OPENSN")

# Test cases
TEST_CASES = [
    (
        os.path.join("HEU_MET_FAST_003", "HEU_MET_FAST_003.py"),  # Directory name, input file name
        {"k-eigenvalue:": 1.001655},                              # key, value to test or None
        12,                                                       # Number of MPI ranks
        1.0e-6                                                    # Tolerance or None
    ),
    (
        os.path.join("OpenSn_Logo_CAD", "opensn.py"),
        None,
        8,
        None,
    ),
    (
        os.path.join("Urban_Source", "urban_source.py"),
        None,
        96,
        None,
    ),
    (
        os.path.join("Six_1g_spherical_benchmarks", "Problem_1.py"),
        {"end-radius:": 1.00, "avg-value:": 8.212354},
        12,
        1.0e-6,
    ),
    (
        os.path.join("Six_1g_spherical_benchmarks", "Problem_2.py"),
        {"end-radius:": 1.00, "avg-value:": 8.212354},
        12,
        1.0e-6,
    ),
    (
        os.path.join("Six_1g_spherical_benchmarks", "Problem_3.py"),
        {"end-radius:": 1.00, "avg-value:": 0.105997},
        12,
        1.0e-6,
    ),
    (
        os.path.join("Six_1g_spherical_benchmarks", "Problem_4.py"),
        {"end-radius:": 1.00, "avg-value:": 12.144526},
        12,
        1.0e-6,
    ),
    (
        os.path.join("Six_1g_spherical_benchmarks", "Problem_5.py"),
        {"end-radius:": 1.00, "avg-value:": 11.660929},
        12,
        1.0e-6,
    ),
    (
        os.path.join("Six_1g_spherical_benchmarks", "Problem_6.py"),
        {"end-radius:": 1.00, "avg-value:": 0.033569},
        12,
        1.0e-6,
    ),
]


@pytest.mark.parametrize("input_file, expected, mpi_ranks, rel_tol", TEST_CASES)
def test_application_output(input_file, expected, mpi_ranks, rel_tol):
    # Ensure OPENSN is set and executable
    assert OPENSN, "OPENSN environment variable must be set"
    assert os.path.isfile(OPENSN) and os.access(OPENSN, os.X_OK), \
        f"OPENSN '{OPENSN}' not found or not executable"

    # Input files and working directory
    abs_input = os.path.abspath(input_file)
    assert os.path.isfile(abs_input), f"Input file not found: {abs_input}"
    cwd = os.path.dirname(abs_input)
    inp = os.path.basename(abs_input)

    # Validate mpi_ranks
    assert isinstance(mpi_ranks, int) and mpi_ranks > 0, \
        f"mpi_ranks must be a positive integer, got {mpi_ranks}"

    # Run opensn
    cmd = ["mpirun", "-n", str(mpi_ranks), OPENSN, "-i", inp]
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"Application exited with code {result.returncode}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}\n"
    )

    if expected is None or rel_tol is None:
        return

    # Validate rel_tol
    assert isinstance(rel_tol, float) and rel_tol > 0.0, \
        f"rel_tol must be a positive float, got {rel_tol}"

    # Combine outputs and process line-by-line attempting to match key/value pairs
    lines = (result.stdout + result.stderr).splitlines()
    matched = False
    for line in lines:
        tokens = [tok for tok in re.split(r"[ ,]+", line.strip()) if tok]
        ok = True
        for key, exp_val in expected.items():
            try:
                idx = tokens.index(key)
            except ValueError:
                ok = False
                break
            if idx + 1 >= len(tokens):
                ok = False
                break
            next_tok = tokens[idx + 1]
            if isinstance(exp_val, (int, float)):
                try:
                    val = float(next_tok)
                except ValueError:
                    ok = False
                    break
                if not (val == pytest.approx(exp_val, rel=rel_tol)):
                    ok = False
                    break
            else:
                if next_tok != exp_val:
                    ok = False
                    break
        if ok:
            matched = True
            break
    if not matched:
        full_output = result.stdout + " " + result.stderr
        pytest.fail(f"No line matches expected key/value pairs: {expected}\n {full_output}")


if __name__ == "__main__":
    if not OPENSN or not os.path.isfile(OPENSN) or not os.access(OPENSN, os.X_OK):
        print(f"Error: OPENSN '{OPENSN}' not found or not executable", file=sys.stderr)
        sys.exit(1)
    ret = pytest.main([os.path.abspath(__file__)])
    sys.exit(0 if ret == 0 else 1)
