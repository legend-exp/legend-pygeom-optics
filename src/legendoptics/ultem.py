"""Ultem (polyetherimide, PEI) structural components

.. [Zhang2020] X. Zhang et al. “Complex refractive indices measurements of polymers in visible and
    near-infrared bands”. In: Appl. Opt. 59, 2337-2344 (2020), https://doi.org/10.1364/AO.383831
"""

from __future__ import annotations

import logging

import pint
from pint import Quantity
from scipy.signal import savgol_filter

from legendoptics import store
from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def ultem_refractive_index() -> tuple[Quantity, Quantity]:
    """Real refractive index reported in [Zhang2020]_, smoothed.

    .. optics-plot::
    """
    λ, r = readdatafile("ultem_rindex.dat")

    r_smoothed = savgol_filter(r.m, 30, 3, mode="mirror") * r.u
    r[λ > 500 * u.nm] = r_smoothed[λ > 500 * u.nm]
    return λ.to("nm"), r


def pyg4_ultem_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given Ultem/PEI material instance.

    See Also
    --------
    .ultem_refractive_index
    """
    λ, r = ultem_refractive_index()
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)
