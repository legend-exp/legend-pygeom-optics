"""Copper."""

from __future__ import annotations

import logging

import pint
from pint import Quantity

from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def copper_reflectivity() -> tuple[Quantity, Quantity]:
    """Reflectivity of copper surfaces.

    Measurements from [Wegmann2017]_ (data points above 300 nm) and [Salamanna2022]_ (data points between 120 and 220 nm, averaged and smoothed).
    The interpolation between both domains is linear, but arbitrary.

    .. optics-plot::
    """
    return readdatafile("cu_reflectivity.dat")


def pyg4_copper_attach_reflectivity(mat, reg) -> None:
    """Attach the optical reflectivity to the given copper material instance.

    See Also
    --------
    .copper_reflectivity
    """
    λ, refl = copper_reflectivity()
    with u.context("sp"):
        mat.addVecPropertyPint("REFLECTIVITY", λ.to("eV"), refl)
