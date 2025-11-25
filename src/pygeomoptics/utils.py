from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pint
import scipy.interpolate
from importlib_resources import files
from pint import Quantity

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def readdatafile(
    filename: str, pkg: str = "pygeomoptics.data", ncols: int = 2
) -> tuple[Quantity, Quantity]:
    """Read ``(x, y)`` data points from `filename` with units.

    Accepted file format ::

        # unit1 unit2
        0.23453 2.3456
        0.49678 3.6841
        ...

    After the first line, comments are also allowed after a ``#`` character.

    Units in the header must be parseable as :mod:`pint` units.

    Parameters
    ----------
    filename
        (relative) file name of the data file, including extension.
    pkg
        python package name used to access data files. Only needs to be set to access
        data files in other packages.
    ncols
        number of columns to read.
    """
    cols = [[] for _ in range(ncols)]
    path = files(pkg).joinpath(filename) if pkg is not None else Path(filename)
    lines = path.read_text().split("\n")
    lines = [line.strip() for line in lines]

    # parse header
    header = lines[0]
    if header[0] != "#" or len(header.lstrip("#").split()) != ncols:
        msg = "input data file does not seem to contain header with (pint) units"
        raise RuntimeError(msg)

    units = header.lstrip("#").split()

    lineno = 0
    for line in lines[1:-1]:
        lineno += 1
        if not line or line[0] == "#":  # skip empty lines and lines with comments.
            continue

        val = line.split()
        if len(val) != ncols:
            msg = f"could not parse line {lineno}: '{line}'"
            raise RuntimeError(msg)

        for i in range(ncols):
            cols[i].append(float(val[i]))

    return tuple(cols[i] * u(units[i]) for i in range(ncols))


class InterpolatingGraph:
    """Linear interpolation between data points, similar to Geant4 default interpolation.

    The data points are given as two 1-dimensional NDArrays with units.
    """

    def __init__(
        self,
        idx: Quantity,
        vals: Quantity,
        min_idx: Quantity | None = None,
        max_idx: Quantity | None = None,
        zero_outside: bool = False,
    ):
        if not isinstance(idx, Quantity) or not isinstance(vals, Quantity):
            msg = "only pint Quantities can be used as index or value"
            raise ValueError(msg)
        if len(idx.shape) != 1 or len(vals.shape) != 1:
            msg = "only 1-dimensional data can be used as index or value"
            raise ValueError(msg)
        if idx.shape != vals.shape:
            msg = "only and index and values of the same shape can be used"
            raise ValueError(msg)

        # Filter the supplied data points.
        f = np.full(idx.shape, True)
        if min_idx is not None:
            f = f & (idx >= min_idx)
        if max_idx is not None:
            f = f & (idx <= max_idx)
        idx = idx[f]
        vals = vals[f]

        self.idx = idx
        self.vals = vals
        self.d_min = min(idx)
        self.d_max = max(idx)
        self.n = len(idx)
        if min(vals).m < 0:
            msg = f"data to be interpolated can only contain positive values. min={min(vals).m}"
            raise ValueError(msg)

        self.fn = scipy.interpolate.interp1d(idx.m, vals.m)

        self.val_min = self.vals[0].m if not zero_outside else 0
        self.val_max = self.vals[-1].m if not zero_outside else 0

    def __call__(self, pts: Quantity) -> Quantity:
        # return first/last value if pts out of defined range
        if isinstance(pts.m, np.ndarray):
            p = pts.to(self.idx.u).m
            return Quantity(
                np.piecewise(
                    p,
                    [
                        p < self.d_min.m,
                        ((p >= self.d_min.m) & (p <= self.d_max.m)),
                        p > self.d_max.m,
                    ],
                    [self.val_min, self.fn, self.val_max],
                ),
                self.vals.u,
            )

        if pts < self.d_min:
            return self.val_min * self.vals.u
        if pts > self.d_max:
            return self.val_max * self.vals.u
        return self.fn(pts.to(self.idx.u).m) * self.vals.u


def g4gps_write_emission_spectrum(
    filename: str,
    output_macro: bool,
    λ_peak: Quantity,
    scint_em: Quantity,
    quantity_name: str,
) -> None:
    """Write a energy spectrum for use with G4GeneralParticleSource.

    It can be used like this in a Geant4 macro:

    .. code ::

        /gps/ene/type     Arb
        /gps/ene/diffspec true
        /gps/hist/type    arb
        /gps/hist/file    <filename>
        /gps/hist/inter   Lin
    """
    with u.context("sp"):
        pointwise = np.array([λ_peak.to("MeV").m, scint_em.m]).T

    # reorder the values to be in ascending energy order.
    sort = np.argsort(pointwise[:, 0])
    pointwise[:, 0] = pointwise[:, 0][sort]
    pointwise[:, 1] = pointwise[:, 1][sort]

    if pointwise.shape[0] > 1024:
        log.warning("G4GeneralParticleSource spectrum can only have 1024 bins.")

    if not output_macro:
        np.savetxt(filename, pointwise)
        return

    with Path.open(filename, "wt") as f:
        f.write(f"# {quantity_name} | pygeomoptics\n\n")
        f.write("/gps/ene/type     Arb\n")
        f.write("/gps/ene/diffspec true\n")
        f.write("/gps/hist/type    arb\n\n")
        for point in pointwise:
            f.write(f"/gps/hist/point   {point[0]} {point[1]}\n")
        f.write("\n/gps/hist/inter   Lin\n")
