"""Steel.

Optical properties of steel are expected to vary depending on surface finish, oxidation, ...
The properties here can only be a guess toward any given special material (i.e. the
GERDA/LEGEND cryostat) until we have a dedicated measurement available.

.. [STEEL1982] Optical constants and spectral selectivity of stainless steel and its oxides
    https://doi.org/10.1063/1.331503
.. [Soto-Oton2026] https://doi.org/10.1088/1748-0221/21/04/C04065
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity

from pygeomoptics import store

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def steel_reflectivity() -> tuple[Quantity, Quantity]:
    """Reflectivity of steel surfaces, modeled after [STEEL1982]_.

    For the value below 200 nm, we just assume a lower reflectivity. This seems to follow the
    general trend of metals. See also [Soto-Oton2026]_, which also measures the overall
    reflectivity over the whole spectrum to be lower. The 20% at 110 nm is a conservative (upper)
    guess.

    .. optics-plot::
    """
    λ = np.array([110, 200, 300, 400, 600, 800]) * u.nm
    refl = np.array([0.2, 0.35, 0.45, 0.55, 0.58, 0.60])
    return λ, refl


def pyg4_steel_attach_reflectivity(mat, reg) -> None:
    """Attach the optical reflectivity to the given steel material instance.

    See Also
    --------
    .steel_reflectivity
    """
    λ, refl = steel_reflectivity()
    with u.context("sp"):
        mat.addVecPropertyPint("REFLECTIVITY", λ.to("eV"), refl)
