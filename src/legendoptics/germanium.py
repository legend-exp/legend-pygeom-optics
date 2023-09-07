"""Germanium.

.. [Wegmann2017] A. Wegmann “Characterization of the liquid argon veto of the GERDA
    experiment and its application for the measurement of the 76Ge half-life”
    (PhD thesis). https://www.mpi-hd.mpg.de/gerda/public/2017/phd2017-anneWegmann.pdf
.. [Salamanna2022] Measurements by E. Bernieri, G. Salamanna, D. Tagnani (Roma Tre & INFN),
    unpublished.
"""

from __future__ import annotations

import logging

import pint
from pint import Quantity

from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def germanium_reflectivity() -> tuple[Quantity, Quantity]:
    """Reflectivity of germanium surfaces.

    Measurements from [Wegmann2017]_ (with GERDA dead-layer Li-doped germanium, at room temperature; data-points above 300nm). Data points between 220 and 120 nm from [Salamanna2022]_ (averaged and smoothed).
    The interpolation between both domains is linear, but arbitrary.

    .. optics-plot::
    """
    return readdatafile("ge_reflectivity.dat")


def pyg4_germanium_attach_reflectivity(mat, reg) -> None:
    """Attach the optical reflectivity to the given germanium material instance.

    See Also
    --------
    .germanium_reflectivity
    """
    λ, refl = germanium_reflectivity()
    with u.context("sp"):
        mat.addVecPropertyPint("REFLECTIVITY", λ.to("eV"), refl)
