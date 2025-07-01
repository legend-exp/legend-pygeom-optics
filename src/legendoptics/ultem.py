"""Ultem (polyetherimide, PEI) structural components

.. [Zhang2020] X. Zhang et al. “Complex refractive indices measurements of polymers in visible and
    near-infrared bands”. In: Appl. Opt. 59, 2337-2344 (2020), https://doi.org/10.1364/AO.383831
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity
from scipy.signal import savgol_filter

from legendoptics import store
from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def ultem_refractive_index() -> tuple[Quantity, Quantity]:
    """Real refractive index reported for PEI in [Zhang2020]_, smoothed.

    The data points below 400 nm are an ad-hoc extension.

    .. optics-plot::
    """
    λ, r = readdatafile("ultem_rindex.dat")

    r_smoothed = savgol_filter(r.m, 30, 3, mode="mirror") * r.u
    r[λ > 500 * u.nm] = r_smoothed[λ > 500 * u.nm]

    λ = np.insert(λ, 0, [350, 375] * u.nm, axis=0)
    r = np.insert(r, 0, [1.684, 1.667], axis=0)

    r, λ = r[λ <= 650 * u.nm], λ[λ <= 650 * u.nm]
    return λ.to("nm"), r


@store.register_pluggable
def ultem_absorption() -> tuple[Quantity, Quantity]:
    """Absorption length, based on the complex refractive index reported for PEI in
    [Zhang2020]_.

    The data points below 400 nm are an ad-hoc extension.

    .. optics-plot::
    """
    λ, κ = readdatafile("ultem_rindex_imag.dat")
    κ, λ = κ[λ <= 650 * u.nm], λ[λ <= 650 * u.nm]

    λ = np.insert(λ, 0, [350, 375] * u.nm, axis=0)
    κ = np.insert(κ, 0, [κ[0], κ[0]], axis=0)

    α = λ.to("mm") / (4 * np.pi * κ)

    return λ.to("nm"), α


def pyg4_ultem_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given Ultem/PEI material instance.

    See Also
    --------
    .ultem_refractive_index
    """
    λ, r = ultem_refractive_index()
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_ultem_attach_absorption(mat, reg) -> None:
    """Attach the absorption length to the given Ultem/PEI material instance.

    See Also
    --------
    .ultem_absorption
    """
    λ, α = ultem_absorption()
    with u.context("sp"):
        mat.addVecPropertyPint("ABSLENGTH", λ.to("eV"), α)
