"""Fused silica (SiO₂) (e.g., Suprasil, ...)

.. [Malitson1965] I.H. Malitson “Interspecimen Comparison of the Refractive Index of Fused Silica”.
    In: J. Opt. Soc. Am. 55, 1205-1209 (1965), https://doi.org/10.1364/JOSA.55.001205
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
def silica_refractive_index(λ: Quantity) -> Quantity:
    """Refractive index reported in [Malitson1965]_.

    .. optics-plot:: {'call_x': True, 'xlim': [200, 650]}
    """
    λ = λ.to("μm").m
    return np.sqrt(
        1
        + 0.6961663 * λ**2 / (λ**2 - 0.0684043**2)
        + 0.4079426 * λ**2 / (λ**2 - 0.1162414**2)
        + 0.8974794 * λ**2 / (λ**2 - 9.896161**2)
    )


def pyg4_silica_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given silica material instance.

    See Also
    --------
    .silica_refractive_index
    """
    from pygeomoptics.pyg4utils import pyg4_sample_λ

    λ = pyg4_sample_λ(200 * u.nm, 650 * u.nm)
    r = silica_refractive_index(λ)
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)
