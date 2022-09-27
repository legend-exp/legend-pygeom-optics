"""
Tetraphenyl-Butadiene wavelength shifter

.. [Francini2013] "VUV-Vis optical characterization of Tetraphenyl-butadiene
    films on glass and specular reflector substrates from
    room to liquid Argon temperature", https://arxiv.org/abs/1304.6117
.. [Benson2017] "Measurements of the intrinsic quantum efficiency and absorption
    length of tetraphenyl butadiene thin films in the vacuum
    ultraviolet regime", https://arxiv.org/abs/1709.05002
.. [Araujo2022] "R&D of Wavelength-Shifting Reflectors and Characterization of
    the Quantum Efficiency of Tetraphenyl Butadiene and Polyethylene
    Naphthalate in Liquid Argon", https://arxiv.org/abs/2112.06675
.. [MolbaseTPB] http://www.molbase.com/en/overview_1450-63-1-moldata-77892.html
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from numpy.typing import NDArray
from pint import Quantity

from legend_optics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def tpb_quantum_efficiency() -> float:
    """Quantum efficiency

    * Current literature value of 0.85 from [Araujo2022]_ at LAr Temperature.
    * Other measurement from [Benson2017]_ reports ~0.6 at room temperature
    """
    return 0.85


def tpb_refractive_index() -> float:
    """Refractive index from [MolbaseTPB]_."""
    return 1.635


def tpb_wls_timeconstant() -> Quantity[float]:
    """Time constant: arbitrary small"""
    return 0.01*u.ns


def tpb_wls_emission() -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    """WLS Emission spectrum

    [Francini2013]_ measure the emission spectrum of TPB (~160 um thick layer) on VM2000
    at an excitation wavelength of 128nm and at 87K, so exactly in our experimental
    conditions. The major differences brougth by the LAr temperature are the vibronic
    structures that modify the shape of the spectrum.
    """
    return readdatafile("tpb_vm2000_wlscomponent.dat")


def tpb_wls_absorption() -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    """Values reported in [Benson2017]_ for TPB evaporated on utraviolet-transmitting acrylic substrate"""
    wvl, absorp = readdatafile("tpb_wlsabslength.dat")
    assert absorp.check('[length]')
    return wvl, absorp

