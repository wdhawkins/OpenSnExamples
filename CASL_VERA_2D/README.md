# CASL VERA Problem 2 Benchmarks

This directory contains OpenSn inputs and support files for the VERA Core Physics
Benchmark Progression Problem 2 two-dimensional fuel-lattice cases.

The reusable benchmark assets live with this example suite: meshes are in
`assets/mesh/casl`, and cross-section libraries are in `assets/xs/casl`.

The checked-in meshes were generated with the included `spydermesh.py` utility.
Use `generate_lattice_meshes.py` to regenerate the full benchmark lattice meshes:

```bash
python generate_lattice_meshes.py
```

By default, this regenerates cases 2A through 2Q into `assets/mesh/casl`. To
regenerate selected cases or write to another directory, pass case names and
`--output-dir`:

```bash
python generate_lattice_meshes.py 2A 2B --output-dir /tmp/casl_meshes
```

## References

- VERA benchmark specification:
  [CASL-U-2012-0131-004, VERA Core Physics Benchmark Progression Problem Specifications](https://corephysics.com/docs/CASL-U-2012-0131-004.pdf)
- Original OpenSn CASL benchmark repository:
  [aly-tamu/VERA-Examples](https://github.com/aly-tamu/VERA-Examples/)

## Problem 2 Reference Eigenvalues

Quick reference for Table P2-5, "Problem 2 Reference Solution Eigenvalue
Results", from the VERA benchmark specification.

The OpenSn regression inputs in this directory are intended to track the VERA
reference eigenvalues to within a few hundred pcm, not to reproduce the
published Monte Carlo values exactly.

| Case | Description | Fuel Temperature | Moderator Density | Reference k-effective |
| --- | --- | ---: | ---: | ---: |
| 2A | No Poisons | 565 K | 0.743 g/cc | 1.182175 +/- 0.000017 |
| 2B | No Poisons | 600 K | 0.661 g/cc | 1.183360 +/- 0.000024 |
| 2C | No Poisons | 900 K | 0.661 g/cc | 1.173751 +/- 0.000023 |
| 2D | No Poisons | 1200 K | 0.661 g/cc | 1.165591 +/- 0.000023 |
| 2E | 12 Pyrex | 600 K | 0.743 g/cc | 1.069627 +/- 0.000024 |
| 2F | 24 Pyrex | 600 K | 0.743 g/cc | 0.976018 +/- 0.000026 |
| 2G | 24 AIC | 600 K | 0.743 g/cc | 0.847695 +/- 0.000025 |
| 2H | 24 B4C | 600 K | 0.743 g/cc | 0.788221 +/- 0.000025 |
| 2I | Instrument Thimble | 600 K | 0.743 g/cc | 1.179916 +/- 0.000024 |
| 2J | Instrument + 24 Pyrex | 600 K | 0.743 g/cc | 0.975193 +/- 0.000025 |
| 2K | Zoned + 24 Pyrex | 600 K | 0.743 g/cc | 1.020063 +/- 0.000025 |
| 2L | 80 IFBA | 600 K | 0.743 g/cc | 1.018915 +/- 0.000024 |
| 2M | 128 IFBA | 600 K | 0.743 g/cc | 0.938796 +/- 0.000025 |
| 2N | 104 IFBA + 20 WABA | 600 K | 0.743 g/cc | 0.869615 +/- 0.000025 |
| 2O | 12 Gadolinia | 600 K | 0.743 g/cc | 1.047729 +/- 0.000024 |
| 2P | 24 Gadolinia | 600 K | 0.743 g/cc | 0.927410 +/- 0.000024 |
| 2Q | Zircaloy Spacer Grid | 565 K | 0.743 g/cc | 1.171940 +/- 0.000016 |
