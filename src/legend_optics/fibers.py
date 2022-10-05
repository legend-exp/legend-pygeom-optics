"""
Scinillating fibers BCF91-A from Saint Gobain.

.. [SaintGobainDataSheet] https://www.crystals.saint-gobain.com/sites/hps-mac3-cma-crystals/files/2021-11/Fiber-Product-Sheet.pdf
"""
from __future__ import annotations

import logging

import pint
from numpy.typing import NDArray
from pint import Quantity

from legend_optics.utils import InterpolatingGraph, readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def fiber_cladding2_refractive_index() -> float:
    """Refractive index of second fiber cladding material [SaintGobainDataSheet]_."""
    return 1.42


def fiber_cladding1_refractive_index() -> float:
    """Refractive index of first fiber cladding material [SaintGobainDataSheet]_."""
    return 1.49


def fiber_core_refractive_index() -> float:
    """Refractive index of fiber core material [SaintGobainDataSheet]_."""
    return 1.6


def fiber_wls_absorption(
    absAt400nm: Quantity[float] = 0.7 * u.mm,
) -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    r"""[SaintGobainDataSheet]_ reports the absorption spectrum for BCF-91A. Knowing that the fibers are 1mm thick one can
    extract the absorption length: starting from the trivial relation:

    :math:`1 - P(E) = \exp(-x/l(E))`

    where :math:`P(E)` is the probability (thus proportional to the absorption spectrum) for a photon
    travelling a distance :math:`x` to be absorbed in the material given the attenuation length
    :math:`l(E)`, one can extract :math:`l(E)` from :math:`P(E)`. By integrating over the thickness of
    the material :math:`L` one obtains:

    :math:`(1 - P(E)) \cdot L = l(E) \cdot (1 - \exp(-L/l(E)))`

    but the problem now is that :math:`l(E)` cannot be extracted analytically (inhomogeneus expression).
    Luigi wrote a Mathematica script that solves it numerically.
    Remeber that the units are arbitrary because the original absorption
    spectrum has arbitrary units.

    Measured an absorption length of 0.7 mm at 400 nm, the spectrum has been rescaled by
    that.
    """
    wvl, absorp = readdatafile("psfibers_wlsabslength.dat")  # arbitrary unit
    assert str(absorp.dimensionality) == "dimensionless"
    # scale factor for absorption lengths (abslength is 0.7mm at 400nm, see above)
    absorp *= absAt400nm / InterpolatingGraph(wvl, absorp)(400 * u.nm)
    return wvl, absorp


def fiber_wls_emission() -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    """[SaintGobainDataSheet]_ reports the emission spectrum for BCF-91A."""
    return readdatafile("psfibers_wlscomponent.dat")


def fiber_wls_timeconstant() -> float:
    """WLS time constant [SaintGobainDataSheet]_."""
    return 12 * u.ns


def fiber_absorption_length() -> float:
    """Absorption length of fiber [SaintGobainDataSheet]_. Note this is a macroscopical value for a 1 mm fiber."""
    return 3.5 * u.m
