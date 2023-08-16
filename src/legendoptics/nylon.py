"""
Nylon from Borexino tank, as described in [Agostini2018]_.

.. [Agostini2018] M. Agostini et al. “The Monte Carlo simulation of the Borexino detector”
    In: Astroparticle Physics 97 (2018). https://doi.org/10.1016/j.astropartphys.2017.10.003

.. [Benziger2007] J. Benziger et al. “The Nylon Scintillator Containment Vessels for the
    Borexino Solar Neutrino Experiment” In: International Journal of Modern Physics A, 29(16)
    (2014). https://doi.org/10.1016/j.nima.2007.08.176
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity

from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def nylon_refractive_index() -> float:
    """Refractive index in near-UV range, from [Benziger2007]_.

    .. optics-const::
    """
    return 1.53


def nylon_absorption() -> tuple[Quantity, Quantity]:
    """Values reported in [Agostini2018]_.

    .. optics-plot::
    """
    wvl, absorp = readdatafile("nylon_absorption.dat")
    assert absorp.check("[length]")
    return wvl, absorp


def pyg4_nylon_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given fiber core material instance.

    See Also
    --------
    .nylon_refractive_index
    """
    λ = np.array([650.0, 115.0]) * u.nm
    r = [nylon_refractive_index()] * 2
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_nylon_attach_absorption(mat, reg) -> None:
    """Attach absorption to the given nylon material instance.

    See Also
    --------
    .nylon_absorption
    """
    λ, absorption = nylon_absorption()
    with u.context("sp"):
        mat.addVecPropertyPint("ABSLENGTH", λ.to("eV"), absorption)
