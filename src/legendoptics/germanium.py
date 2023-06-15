"""Germanium.

.. [Wegmann2017] A. Wegmann “Characterization of the liquid argon veto of the GERDA
    experiment and its application for the measurement of the 76Ge half-life”
    (PhD thesis). https://www.mpi-hd.mpg.de/gerda/public/2017/phd2017-anneWegmann.pdf
"""

from __future__ import annotations

import logging

import pint
from pint import Quantity

from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def germanium_reflectivity() -> tuple[Quantity, Quantity]:
    """Measurements from [Wegmann2017]_ (with GERDA dead-layer Li-doped germanium, at room temperature)."""
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
