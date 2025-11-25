"""Polyethylene Naphthalate (PEN) plastic scintillator and wavelength shifter.

.. [Manzanillas2022] L. Manzanillas et al. “Optical properties of low background PEN structural components
    for the LEGEND-200 experiment”. In: JINST 17 P09007 (2022), https://doi.org/10.1088/1748-0221/17/09/P09007
.. [Hong2017] N. Hong et al. “Mueller matrix characterization of flexible plastic substrates”.
    In: Applied Surface Science, 421:518-528 (2017), https://doi.org/10.1016/j.apsusc.2017.01.276
.. [Ouchi2006] I. Ouchi et al. “Features of Fluorescence Spectra of Polyethylene 2,6-Naphthalate Films”
    In: Journal of Applied Polymer Science, Vol. 105, 114-121 (2007), https://doi.org/10.1002/app.26085
.. [Hackett2024] B. Hackett et al. “Light response of poly(ethylene 2,6-naphthalate) to neutrons”. In:
    Nuclear Inst. and Methods in Physics Research, A (2024), https://doi.org/10.1016/j.nima.2024.169900
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity

from pygeomoptics import store
from pygeomoptics.scintillate import ScintConfig, ScintParticle
from pygeomoptics.utils import (
    InterpolatingGraph,
    g4gps_write_emission_spectrum,
    readdatafile,
)

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def pen_quantum_efficiency() -> float:
    """Quantum efficiency, from [Araujo2022]_ at LAr temperature.

    .. optics-const::
    """
    return 0.69


@store.register_pluggable
def pen_refractive_index() -> float:
    """Refractive index from [Hong2017]_.

    .. optics-const::
    """
    return 1.51


@store.register_pluggable
def pen_scint_timeconstant() -> Quantity:
    """Time constant, from [Manzanillas2022]_.

    .. optics-const::
    """
    return 25.3 * u.ns


@store.register_pluggable
def pen_scint_light_yield() -> Quantity:
    """PEN scintillation yield for electrons, from [Manzanillas2022]_.

    .. optics-const::
    """
    return 5440 / u.MeV


@store.register_pluggable
def pen_wls_emission() -> tuple[Quantity, Quantity]:
    """WLS Emission spectrum.

    [Leonhardt2024]_ measure the emission spectrum of a PEN sample made from the same pellets as
    in LEGEND-200 sample at an excitation wavelength of 128nm and at 87K, so exactly in our
    experimental conditions.

    .. note::
        Data points below 375 nm are an ad-hoc continuation of the measured data to zero, avoiding
        a steep step. They are not based on actual measurements.

    .. optics-plot::
    """
    return readdatafile("pen_wlscomponent.dat")


@store.register_pluggable
def pen_absorption() -> tuple[Quantity, Quantity]:
    """Bulk absorption reported in [Manzanillas2022]_.

    .. optics-plot::
    """
    wvl, absorp = readdatafile("pen_abslength.dat")
    assert absorp.check("[length]")
    return wvl, absorp


@store.register_pluggable
def pen_wls_absorption() -> tuple[Quantity, Quantity]:
    """WLS absorption of PEN.

    .. warning::
        There is no measurement of PEN similar to ours available, so this step-function is improvised.
        For geometries with thick PEN objects, the absorption length should not matter too much—in a
        certain range, all light will be absorbed anyway.

    The absorbing range and approximate magnitude have been extracted from [Ouchi2006]_, figure 1.

    .. optics-plot:: {'yscale': 'log'}
    """
    wvl = np.array([70, 71, 380, 381]) * u.nm
    absorp = np.array([1e3, 2e-8, 2e-8, 1e3]) * u.m  # 1e3 is "infinity"
    assert absorp.check("[length]")
    return wvl, absorp


@store.register_pluggable
def pen_scintillation_params() -> ScintConfig:
    """Get a :class:`ScintConfig` object for PEN.

    This implements the measured electron light yield. The light yield for other particle
    types (protons, alphas, ions) is derived from the Birk's constant in [Hackett2024]_
    assuming some common LET ranges. The quenching is not implemented in an energy-dependent way.

    .. optics-const::

    See Also
    --------
    .pen_scint_light_yield
    """
    # setup scintillation response just for electrons.
    return ScintConfig(
        flat_top=pen_scint_light_yield(),
        fano_factor=None,
        particles=[
            ScintParticle("electron", yield_factor=1, exc_ratio=None),
            # more particle types have to be added to simulate decays in or nearby PEN plates.
            # this is an ad-hoc guess from the Birk's constant in and common LET ranges from
            # ASTAR/PSTAR/ESTAR. Ions should have a larger LET than alphas.
            #
            # TODO: implementing a full energy-dependency with this "scintillation by particle
            # type" mechanism is impossible. Investigate the use of Birk's law in Geant4.
            ScintParticle("proton", yield_factor=0.5, exc_ratio=None),
            ScintParticle("alpha", yield_factor=0.1, exc_ratio=None),
            ScintParticle("ion", yield_factor=0.05, exc_ratio=None),
        ],
    )


def pyg4_pen_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given PEN material instance.

    See Also
    --------
    .pen_refractive_index
    """
    λ = np.array([650.0, 115.0]) * u.nm
    r = [pen_refractive_index()] * 2
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_pen_attach_attenuation(mat, reg) -> None:
    """Attach bulk absorption properties to the given PEN material instance.

    See Also
    --------
    .pen_absorption
    """
    with u.context("sp"):
        λ_abs, absorption = pen_absorption()
        # set absorption for lowest wavelength to "infinity".
        absorption[np.argmin(λ_abs)] = 1e3 * u.m
        mat.addVecPropertyPint("ABSLENGTH", λ_abs.to("eV"), absorption)


def pyg4_pen_attach_wls(mat, reg, quantum_efficiency: bool | float = True) -> None:
    """Attach wavelength shifting properties to the given PEN material instance.

    Parameters
    ----------
    quantum_efficiency
        If `False`, disable attaching any photon number information. If `True`, use the values
        from .pen_quantum_efficiency. If specified as a number, directly attach this number as
        mean number of emitted photons.

    See Also
    --------
    .pen_wls_absorption
    .pen_wls_emission
    .pen_scint_timeconstant
    .pen_quantum_efficiency
    """
    from pygeomoptics.pyg4utils import pyg4_sample_λ

    λ_abs, absorption = pen_wls_absorption()

    λ_scint = pyg4_sample_λ(350 * u.nm, 650 * u.nm, 800)  # sample more points for WLS.
    emission = InterpolatingGraph(*pen_wls_emission(), min_idx=350 * u.nm)(λ_scint)
    # make sure that the scintillation spectrum is zero at the boundaries.
    emission[0] = 0
    emission[-1] = 0

    with u.context("sp"):
        mat.addVecPropertyPint("WLSABSLENGTH", λ_abs.to("eV"), absorption)
        mat.addVecPropertyPint("WLSCOMPONENT", λ_scint.to("eV"), emission)

    mat.addConstPropertyPint("WLSTIMECONSTANT", pen_scint_timeconstant())
    if quantum_efficiency is True:
        mat.addConstPropertyPint("WLSMEANNUMBERPHOTONS", pen_quantum_efficiency())
    elif isinstance(quantum_efficiency, float):
        mat.addConstPropertyPint("WLSMEANNUMBERPHOTONS", quantum_efficiency)


def pyg4_pen_attach_scintillation(mat, reg) -> None:
    """Attach Geant4 properties for PEN scintillation response to the given material instance.

    .. note:: This currently only adds scintillation for energy deposited by electrons.

    See Also
    --------
    .pen_scint_light_yield
    .pen_wls_emission
    .pen_scint_timeconstant
    """
    from pygeomoptics.pyg4utils import pyg4_def_scint_by_particle_type, pyg4_sample_λ

    # sample the measured emission spectrum.
    λ_scint = pyg4_sample_λ(350 * u.nm, 650 * u.nm, 200)
    scint_em = InterpolatingGraph(*pen_wls_emission(), min_idx=350 * u.nm)(λ_scint)
    # make sure that the scintillation spectrum is zero at the boundaries.
    scint_em[0] = 0
    scint_em[-1] = 0

    with u.context("sp"):
        mat.addVecPropertyPint("SCINTILLATIONCOMPONENT1", λ_scint.to("eV"), scint_em)

    mat.addConstPropertyPint("SCINTILLATIONTIMECONSTANT1", pen_scint_timeconstant())

    # We do not know the PEN fano factor. Geant4 calculates σ = RESOLUTIONSCALE × √mean,
    # so we set it to 1; otherwise Geant4 will crash.
    mat.addConstPropertyPint("RESOLUTIONSCALE", 1)

    pyg4_def_scint_by_particle_type(mat, pen_scintillation_params())


def g4gps_pen_emissions_spectrum(filename: str, output_macro: bool) -> None:
    """Write a PEN emission energy spectrum for G4GeneralParticleSource.

    See Also
    --------
    .pen_wls_emission
    utils.g4gps_write_emission_spectrum
    """
    from pygeomoptics.pyg4utils import pyg4_sample_λ

    # sample the measured emission spectrum.
    λ_scint = pyg4_sample_λ(350 * u.nm, 650 * u.nm, 200)
    scint_em = InterpolatingGraph(*pen_wls_emission(), min_idx=350 * u.nm)(λ_scint)
    # make sure that the scintillation spectrum is zero at the boundaries.
    scint_em[0] = 0
    scint_em[-1] = 0

    g4gps_write_emission_spectrum(
        filename, output_macro, λ_scint, scint_em, "pen_emissions_spectrum"
    )
