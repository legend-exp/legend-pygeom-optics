from __future__ import annotations

import logging

import numpy as np
import pint
from numpy.typing import NDArray
from pint import Quantity

from legend_optics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def fiber_cl2_refractive_index() -> float:
    """Refractive index of second fiber cladding material [SaintGobainDataSheet]_."""
    return 1.42


def fiber_cl1_refractive_index() -> float:
    """Refractive index of first fiber cladding material [SaintGobainDataSheet]_."""
    return 1.49


def fiber_core_refractive_index() -> float:
    """Refractive index of fiber core material [SaintGobainDataSheet]_."""
    return 1.6


def fiber_wls_absorption(absAt400nm: Quantity[float] = 0.7*u.mm) -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    """
    [SaintGobainDataSheet]_ reports the absorption spectrum for BCF-91A. Knowing that the fibers are 1mm thick one can
    extract the absorption length: starting from the trivial relation:

    1 - P(E) = exp(-x/l(E))

    where P(E) is the probability (thus proportional to the absorption spectrum) for a photon
    travelling a distance x to be absorbed in the material given the attenuation length l(E), one
    can extract l(E) from P(E).  By integrating over the thickness of the material L one obtains:

    (1 - P(E)) * L = l(E) * (1 - exp(-L/l(E)))

    but the problem now is that l(E) cannot be extracted analytically (inhomogeneus expression).
    I wrote a Mathematica script that solves it numerically.
    Remeber that the units are arbitrary because the original absorption
    spectrum has arbitrary units.

    Measured an absorption length of 0.7 mm at 400 nm, the spectrum has been rescaled by
    that. Reference in Raphael Kneissl's bachelor thesis

    .. [SaintGobainDataSheet] https://www.crystals.saint-gobain.com/sites/hps-mac3-cma-crystals/files/2021-11/Fiber-Product-Sheet.pdf
    """
    wvl, absorp = readdatafile('psfibers_wlsabslength.dat') # arbitrary unit
    assert(absorp.check('dimensionless'))
    # scale factor for absorption lengths (abslength is 0.7mm at 400nm, see above)
    absorp *= absAt400nm / FibersAbsorptionGr(400)
    return wvl, absorp


def fiber_wls_emission() -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    """[SaintGobainDataSheet]_ reports the emission spectrum for BCF-91A."""
    return readdatafile('psfibers_wlscomponent.dat')


def fiber_wls_timeconstant() -> float:
    """WLS time constant [SaintGobainDataSheet]_."""
    return 12*u.ns


def fiber_absorption_length() -> float:
    """Absorption length of fiber [SaintGobainDataSheet]_. Note this is a macroscopical value for a 1 mm fiber."""
    return 3.5*u.m
