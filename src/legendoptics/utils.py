from __future__ import annotations

import numpy as np
import pint
import scipy.interpolate
from importlib_resources import files
from numpy.typing import NDArray
from pint import Quantity

u = pint.get_application_registry()


def readdatafile(filename: str) -> tuple[NDArray, NDArray]:
    """Read ``(x, y)`` data points from `filename` with units.

    Accepted file format ::

        # unit1 unit2
        0.23453 2.3456
        0.49678 3.6841
        ...

    Units in the header must be parseable as :mod:`pint` units.
    """
    x = []
    y = []
    lines = files("legendoptics.data").joinpath(filename).read_text().split("\n")

    # parse header
    header = lines[0].lstrip()
    if header[0] != "#":
        raise RuntimeError(
            "input data file does not seem to contain header with (pint) units"
        )

    units = header.lstrip("#").split()

    lineno = 0
    for line in lines[1:-1]:
        lineno += 1
        if not line:
            continue

        val = line.split()
        if len(val) < 2:
            raise RuntimeError(f"could not parse line {lineno}: '{line}'")

        x.append(float(val[0]))
        y.append(float(val[1]))

    return (x * u[units[0]], y * u[units[1]])


class InterpolatingGraph:
    """Linear interpolation between data points, similar to Geant4 default interpolation.

    The data points are given as two 1-dimensional NDArrays with units.
    """

    def __init__(self, idx: Quantity[NDArray], vals: Quantity[NDArray]):
        self.idx = idx
        self.vals = vals
        self.d_min = min(idx)
        self.d_max = max(idx)
        self.n = len(idx)
        assert min(vals).m >= 0  # We only want positive values in the spectra
        fn = scipy.interpolate.interp1d(idx.m, vals.m)
        self.fn = lambda l: u.Quantity(fn(l.to(self.idx.u).m), self.vals.u)

    def __call__(self, pts: Quantity[float | NDArray]) -> Quantity[float | NDArray]:
        # return first/last value if pts out of defined range
        if isinstance(pts, np.ndarray):
            return np.piecewise(
                pts,
                [
                    pts < self.d_min,
                    ((pts >= self.d_min) & (pts <= self.d_max)),
                    pts > self.d_max,
                ],
                [self.vals.iloc[0], self.fn, self.vals.iloc[-1]],
            )
        if pts < self.d_min:
            return self.vals.iloc[0]
        if pts > self.d_max:
            return self.vals.iloc[-1]
        return self.fn(pts)
