"""
Tyvek reflector.
"""

from __future__ import annotations

import logging

import pint
from pint import Quantity

from legendoptics import store
from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def tyvek_reflectivity() -> tuple[Quantity, Quantity]:
    """Tyvek reflectivity from [Janacek2012]_.

       A little bit more conservative than in the paper
       (-1% reflectivity compared to the values in the paper).

    .. optics-plot::
    """
    return readdatafile("tyvek_reflectivity.dat")


def pyg4_tyvek_attach_reflectivity(mat, reg, reflectivity_scale: float = 1) -> None:
    """Attach the optical reflectivity to the given material instance.

    Parameters
    ----------
    reflectivity_scale
        Global scale for tyvek reflectivity.

    See Also
    --------
    .tyvek_reflectivity
    """
    λ, refl = tyvek_reflectivity()
    with u.context("sp"):
        mat.addVecPropertyPint("REFLECTIVITY", λ.to("eV"), reflectivity_scale * refl)
