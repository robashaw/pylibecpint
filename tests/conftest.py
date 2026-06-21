"""Test utilities matching libecpint's testutil.hpp verification logic."""

from pathlib import Path

import numpy as np


REFERENCE_DIR = Path(__file__).parent / "reference"


def load_reference(filename: str) -> np.ndarray:
    """Load reference values from a .output file."""
    path = REFERENCE_DIR / filename
    values = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line:
            values.append(float(line))
    return np.array(values)


def check_against_reference(
    results: np.ndarray,
    filename: str,
    thresh: float = 1e-5,
    precision: float = 1e-10,
    maxrel: float = 1e-3,
) -> None:
    """Compare results to reference file, matching testutil.hpp logic.

    Raises AssertionError with diagnostics on failure.
    """
    benchmark = load_reference(filename)
    results_flat = np.asarray(results).ravel()

    assert len(results_flat) == len(benchmark), (
        f"Size mismatch: got {len(results_flat)}, expected {len(benchmark)}"
    )

    maxval = np.max(np.abs(benchmark))
    sig_floor = max(precision, 1e-6 * maxval)
    abs_tol = max(precision, 1e-6 * maxval)

    significant = np.abs(benchmark) > sig_floor
    abserrors = np.abs(benchmark - results_flat)

    if np.any(significant):
        relerrors_sig = abserrors[significant] / np.abs(benchmark[significant])
        avg_relerror = np.mean(relerrors_sig)
        worst_relerror = np.max(relerrors_sig)
    else:
        avg_relerror = 0.0
        worst_relerror = 0.0

    insignificant = ~significant
    if np.any(insignificant):
        worst_abs = np.max(abserrors[insignificant])
    else:
        worst_abs = 0.0

    assert avg_relerror <= thresh, (
        f"Average relative error {avg_relerror:.2e} exceeds threshold {thresh:.2e}"
    )
    assert worst_relerror <= maxrel, (
        f"Worst relative error {worst_relerror:.2e} exceeds maximum {maxrel:.2e}"
    )
    assert worst_abs <= abs_tol, (
        f"Worst absolute error on near-zero element {worst_abs:.2e} exceeds tolerance {abs_tol:.2e}"
    )
