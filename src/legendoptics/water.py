"""
High puritiy water for LEGEND-200 watertank.

.. [Mason2016] John D. Mason, Michael T. Cone, and Edward S. Fry, “Ultraviolet (250-550 nm) absorption
    spectrum of pure water”. In: Appl. Opt. 55, 7163-7172 (2016). https://doi.org/10.1364/AO.55.007163
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity

from legendoptics import store

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def water_refractive_index() -> float:
    """Refractive index.

    .. optics-const::
    """
    return 1.33


@store.register_pluggable
def water_absorption() -> tuple[Quantity, Quantity]:
    """Ultra pure water absorption lengths, from [Mason2016]_.

    .. optics-plot::
    """

    wvl = [
        600,
        550,
        500,
        450,
        400,
        350,
        300,
        250,
        200,
        150,
        100,
    ] * u.nm

    abs = [
        10 * 1000,  # 10 m
        20 * 1000,  # 20 m
        50 * 1000,  # 50 m
        100 * 1000,  # 100 m
        100 * 1000,  # 100 m
        100 * 1000,  # 100 m
        90 * 1000,  # 90 m
        20 * 1000,  # 20 m
        1 * 1000,  # 1 m
        0.001,  # 0.001 mm
        0.0001,  # 0.0001 mm
    ] * u.mm

    assert abs.check("[length]")
    return wvl, abs


def pyg4_water_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given water material instance.

    See Also
    --------
    .water_refractive_index
    """
    λ = np.array([100, 600]) * u.nm
    r = [water_refractive_index()] * 2
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_water_attach_absorption(mat, reg) -> None:
    """Attach absorption to the given water material instance.

    See Also
    --------
    .water_absorption
    """
    λ, absorption = water_absorption()
    with u.context("sp"):
        mat.addVecPropertyPint("ABSLENGTH", λ.to("eV"), absorption)
