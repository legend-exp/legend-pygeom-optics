"""
Tetratex reflector.

.. [Janecek2012] M. Janecek, "Reflectivity spectra for commonly used reflectors",
    https://www.osti.gov/servlets/purl/1184400
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity

from legendoptics import store
from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def tetratex_reflectivity() -> tuple[Quantity, Quantity]:
    """Tetratex reflectivity from [Janecek2012]_.

    He measures the reflectivity of 2 and 4 superimposed layers of 160um thick
    Tetratex. As our layer in GERDA/LEGEND is 254um thick I'm taking here his results
    for the two superimposed foils (= 320um). So, in reality, the reflectivity
    of our foil should be (negligibly) smaller.

    The measured spectrum only has data points above 250 nm. Below that, the reflectivity
    is expected to drop; however, the value is unknown ([Araujo2022]_ finds an 90 %-CL
    upper limit of ~ 17 %). We set an arbitrary reflectivity of 10 % here.

    .. optics-plot::
    """
    λ, r = readdatafile("tetratex_reflectivity.dat")
    # add the lower reflectivity for the VUV range.
    λ = np.insert(λ, 0, 150 * u.nm, axis=0)
    r = np.insert(r, 0, 0.1, axis=0)
    return λ, r


def pyg4_tetratex_attach_reflectivity(mat, reg, reflectivity_scale: float = 1) -> None:
    """Attach the optical reflectivity to the given tetratex material instance.

    Parameters
    ----------
    reflectivity_scale
        Global scale for tetratex reflectivity.

    See Also
    --------
    .tetratex_reflectivity
    """
    λ, refl = tetratex_reflectivity()
    with u.context("sp"):
        mat.addVecPropertyPint("REFLECTIVITY", λ.to("eV"), reflectivity_scale * refl)
