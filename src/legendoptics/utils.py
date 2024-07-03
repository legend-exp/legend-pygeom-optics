from __future__ import annotations

import logging

import numpy as np
import pint
import scipy.interpolate
from importlib_resources import files
from pint import Quantity

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def readdatafile(filename: str) -> tuple[Quantity, Quantity]:
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
        msg = "input data file does not seem to contain header with (pint) units"
        raise RuntimeError(msg)

    units = header.lstrip("#").split()

    lineno = 0
    for line in lines[1:-1]:
        lineno += 1
        if not line:
            continue

        val = line.split()
        if len(val) < 2:
            msg = f"could not parse line {lineno}: '{line}'"
            raise RuntimeError(msg)

        x.append(float(val[0]))
        y.append(float(val[1]))

    return (x * u(units[0]), y * u(units[1]))


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
    ):
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
        assert min(vals).m >= 0  # We only want positive values in the spectra
        self.fn = scipy.interpolate.interp1d(idx.m, vals.m)

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
                    [self.vals[0].m, self.fn, self.vals[-1].m],
                ),
                self.vals.u,
            )

        if pts < self.d_min:
            return self.vals.iloc[0]
        if pts > self.d_max:
            return self.vals.iloc[-1]
        return self.fn(pts.to(self.idx.u).m) * self.vals.u


def g4gps_write_emission_spectrum(
    filename: str, λ_peak: Quantity, scint_em: Quantity
) -> None:
    """Write a energy spectrum for use with G4GeneralParticleSource

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

    if pointwise.shape[0] > 1024:
        log.warning("G4GeneralParticleSource spectrum can only have 1024 bins.")

    np.savetxt(filename, pointwise)
