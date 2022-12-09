"""
Nylon from Borexino tank, as described in [Agostini2018]_.

.. [Agostini2018] M. Agostini et al. “The Monte Carlo simulation of the Borexino detector”
    In: Astroparticle Physics 97 (2018). https://doi.org/10.1016/j.astropartphys.2017.10.003

.. [Benziger2007] J. Benziger et al. “The Nylon Scintillator Containment Vessels for the
    Borexino Solar Neutrino Experiment” In: nternational Journal of Modern Physics A, 29(16)
    (2014). https://doi.org/10.1016/j.nima.2007.08.176
"""

from __future__ import annotations

import logging

import pint
from numpy.typing import NDArray
from pint import Quantity

from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def nylon_refractive_index() -> float:
    """Refractive index in near-UV range, from [Benziger2007]_."""
    return 1.53


def nylon_absorption() -> tuple[Quantity, Quantity]:
    """Values reported in [Agostini2018]_."""
    wvl, absorp = readdatafile("nylon_absorption.dat")
    assert absorp.check("[length]")
    return wvl, absorp
