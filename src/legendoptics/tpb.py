"""
Tetraphenyl-Butadiene wavelength shifter.

.. [Francini2013] R. Francini et al. “VUV-Vis optical characterization of Tetraphenyl-butadiene films
    on glass and specular reflector substrates from room to liquid Argon temperature.”
    In: Journal of Instrumentation 8.09 (Sept. 2013). https://doi.org/10.1088/1748-0221/8/09/p09006.
.. [Benson2018] C. Benson et al. “Measurements of the intrinsic quantum eﬀiciency and absorption
    length of tetraphenyl butadiene thin films in the vacuum ultraviolet regime.”
    In: The European Physical Journal C 78.4 (Apr. 2018).
    https://doi.org/10.1140/epjc/s10052-018-5807-z
.. [Araujo2022] G. R. Araujo et al. “R&D of wavelength-shifting reflectors and characterization of
    the quantum eﬀiciency of tetraphenyl butadiene and polyethylene naphthalate in
    liquid argon.” In: The European Physical Journal C 82.5 (May 2022).
    https://doi.org/10.1140/epjc/s10052-022-10383-0
.. [MolbaseTPB] http://www.molbase.com/en/overview_1450-63-1-moldata-77892.html
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity

from legendoptics.utils import InterpolatingGraph, readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def tpb_quantum_efficiency() -> float:
    """Quantum efficiency.

    * Current literature value of 0.85 from [Araujo2022]_ at LAr Temperature.
    * Other measurement from [Benson2018]_ reports ~0.6 at room temperature

    .. optics-const::
    """
    return 0.85


def tpb_refractive_index() -> float:
    """Refractive index from [MolbaseTPB]_.

    .. optics-const::
    """
    return 1.635


def tpb_wls_timeconstant() -> Quantity:
    """Time constant: arbitrary small.

    .. optics-const::
    """
    return 0.01 * u.ns


def tpb_wls_emission() -> tuple[Quantity, Quantity]:
    """WLS Emission spectrum.

    [Francini2013]_ measure the emission spectrum of TPB (~160 um thick layer) on VM2000
    at an excitation wavelength of 128nm and at 87K, so exactly in our experimental
    conditions. The major differences brougth by the LAr temperature are the vibronic
    structures that modify the shape of the spectrum.

    .. optics-plot::
    """
    return readdatafile("tpb_vm2000_wlscomponent.dat")


def tpb_polystyrene_wls_emission() -> tuple[Quantity, Quantity]:
    """WLS Emission spectrum for TPB in a polystyrene matrix.

    [Francini2013]_ measure the emission spectrum of TPB in a polystyrene matrix
    at an excitation wavelength of 128nm and at 87K, so exactly in our experimental
    conditions. The major differences brougth by the embedding in the PS matrix is the
    shift of the main emission peak and a loss of vibronic structures.

    .. optics-plot::
    """
    return readdatafile("tpb_polystyrene_wlscomponent.dat")


def tpb_wls_absorption() -> tuple[Quantity, Quantity]:
    """Values reported in [Benson2018]_ for TPB evaporated on utraviolet-transmitting acrylic substrate.

    .. optics-plot:: {'yscale': 'log'}
    """
    wvl, absorp = readdatafile("tpb_wlsabslength.dat")
    assert absorp.check("[length]")
    return wvl, absorp


def pyg4_tpb_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given tpb material instance.

    See Also
    --------
    .tpb_refractive_index
    """
    λ = np.array([650.0, 115.0]) * u.nm
    r = [tpb_refractive_index()] * 2
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_tpb_attach_wls(
    mat,
    reg,
    quantum_efficiency: bool | float = True,
    emission_spectrum: str = "default",
) -> None:
    """Attach wavelength shifting properties to the given tpb material instance.

    Parameters
    ----------
    quantum_efficiency
        If `False`, disable attaching any photon number information. If `True`, use the values
        from .tpb_quantum_efficiency. If specified as a number, directly attach this number as
        mean number of emitted photons.
    emission_spectrum
        either `default` or `polystyrene_matrix`

    See Also
    --------
    .tpb_wls_absorption
    .tpb_wls_emission
    .tpb_polystyrene_wls_emission
    .tpb_wls_timeconstant
    .tpb_quantum_efficiency
    """
    from legendoptics.pyg4utils import pyg4_sample_λ

    if emission_spectrum not in ["default", "polystyrene_matrix"]:
        msg = "invalid parameter value of emission_spectrum"
        raise ValueError(msg)

    emission_fn = tpb_wls_emission
    if emission_spectrum == "polystyrene_matrix":
        emission_fn = tpb_polystyrene_wls_emission

    λ_full = pyg4_sample_λ(112 * u.nm, 650 * u.nm, 800)

    absorption = InterpolatingGraph(*tpb_wls_absorption())(λ_full)
    emission = InterpolatingGraph(*emission_fn())(λ_full)
    # make sure that the scintillation spectrum is zero at the boundaries.
    emission[0] = 0
    emission[-1] = 0

    with u.context("sp"):
        mat.addVecPropertyPint("WLSABSLENGTH", λ_full.to("eV"), absorption)
        mat.addVecPropertyPint("WLSCOMPONENT", λ_full.to("eV"), emission)

    mat.addConstPropertyPint("WLSTIMECONSTANT", tpb_wls_timeconstant())
    if quantum_efficiency is True:
        mat.addConstPropertyPint("WLSMEANNUMBERPHOTONS", tpb_quantum_efficiency())
    if isinstance(quantum_efficiency, float):
        mat.addConstPropertyPint("WLSMEANNUMBERPHOTONS", quantum_efficiency)
