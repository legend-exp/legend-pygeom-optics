from __future__ import annotations

import logging

import numpy as np
import pint
from numpy.typing import NDArray
from pint import Quantity

from legend_optics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def tetratex_reflectivity() -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    """
    Tetratex reflectivity from [Janacek2012]_.

    He measures the reflectivity of 2 and 4 superimposed layers of 160um thick
    Tetratex. As our layer in GERDA/LEGEND is 254um thick I'm taking here his results
    for the two superimposed foils (= 320um). So, in reality, the reflectivity
    of our foil should be (negligibly) smaller.

    .. [Janacek2012] M. Janacek, "Reflectivity spectra for commonly used reflectors", https://www.osti.gov/servlets/purl/1184400
    """
    return readdatafile("tpb_wlsabslength.dat")
