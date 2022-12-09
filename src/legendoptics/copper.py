"""
Copper
"""
    
from __future__ import annotations

import logging

import pint
from numpy.typing import NDArray
from pint import Quantity

from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def copper_reflectivity() -> tuple[Quantity, Quantity]:
    """Measurements from [Wegmann2017]_."""
    return readdatafile("cu_reflectivity.dat")
