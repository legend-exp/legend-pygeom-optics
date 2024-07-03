"""Liquid Argon (LAr).

.. [Heindl2010] T. Heindl et al. “The scintillation of liquid argon.” In: EPL 91.6 (Sept. 2010).
    https://doi.org/10.1209/0295-5075/91/62002
.. [Doke1976] Doke et al. “Estimation of Fano factors in liquid argon, krypton, xenon and
    xenon-doped liquid argon.” NIM 134 (1976)353, https://doi.org/10.1016/0029-554X(76)90292-5
.. [Bideau-Mehu1980] Bideau-Mehu et al. “Measurement of refractive indices of neon, argon, krypton
    and xenon in the 253.7-140.4 nm wavelength range. Dispersion relations and
    estimated oscillator strengths of the resonance lines.” In: Journal of Quantitative
    Spectroscopy and Radiative Transfer 25.5 (1981)).
    https://doi.org/10.1016/0022-4073(81)90057-1
.. [Seidel2002] G. M. Seidel at al. “Rayleigh scattering in rare-gas liquids.” In: Nuclear
    Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers
    Detectors and Associated Equipment 489.1-3 (Aug. 2002).
    https://doi.org/10.1016/s0168-9002(02)00890-2
.. [Babicz2020] M. Babicz et al. “A measurement of the group velocity of scintillation light in liquid
    argon.” In: Journal of Instrumentation 15.09 (Sept. 2020).
    https://doi.org/10.1088/1748-0221/15/09/P09009
.. [Doke2002] T. Doke et al. “Absolute Scintillation Yields in Liquid Argon and Xenon for Various Particles”
    Jpn. J. Appl. Phys. 41 1538, https://doi.org/10.1143/JJAP.41.1538
.. [Hitachi1983] A. Hitachi et al. “Effect of ionization density on the time dependence of luminescence
    from liquid argon and xenon.” In: Phys. Rev. B 27 (9 May 1983), pp. 5279-5285,
    https://doi.org/10.1103/PhysRevB.27.5279
.. [Schwarz2024] M. Schwarz “Tracing impurities and illuminating their impact: Surveying and characterizing
    liquid argon with LLAMA for LEGEND and beyond” (PhD thesis). https://mediatum.ub.tum.de/?id=1741523
"""

from __future__ import annotations

import logging
from typing import Literal, NamedTuple

import numpy as np
import pint
from pint import Quantity

from legendoptics.scintillate import ScintConfig, ScintParticle
from legendoptics.utils import (
    InterpolatingGraph,
    g4gps_write_emission_spectrum,
    readdatafile,
)

log = logging.getLogger(__name__)
u = pint.get_application_registry()

ArDielectricMethods = Literal["cern2020", "bideau-mehu"]
ArLifetimeMethods = Literal["legend200-llama"]


class ArScintLiftime(NamedTuple):
    singlet: Quantity
    triplet: Quantity

    def as_tuple(self) -> tuple[Quantity, Quantity]:
        return (self.singlet, self.triplet)


def lar_dielectric_constant_bideau_mehu(
    λ: Quantity,
) -> Quantity:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    From the Bideau-Sellmeier formula [Bideau-Mehu1980]_ in gaseous argon,
    density-corrected for liquid argon.

    .. optics-plot:: {'call_x': True}
    """
    if not λ.check("[length]"):
        msg = "input does not look like a wavelength"
        raise ValueError(msg)

    if np.any(λ < 110 * u.nm):
        msg = f"this parametrization is not meaningful below {110*u.nm}"
        raise ValueError(msg)

    # equation for n-1
    ϵ = 1.2055e-2 * (
        0.2075 * λ**2 / (91.012 * λ**2 - 1 * u.um**2)
        + 0.0415 * λ**2 / (87.892 * λ**2 - 1 * u.um**2)
        + 4.3330 * λ**2 / (214.02 * λ**2 - 1 * u.um**2)
    )
    ϵ *= 2 / 3  # Bideau-Sellmeier -> Clausius-Mossotti
    ϵ *= 1.396 / 1.66e-3  # density correction (Ar gas -> LAr liquid)

    # solve Clausius-Mossotti
    return (1 + 2 * ϵ) / (1 - ϵ)


def lar_dielectric_constant_cern2020(
    λ: Quantity,
) -> Quantity:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    From [Babicz2020]_ (measurements in LAr).

    .. optics-plot:: {'call_x': True}
    """
    if not λ.check("[length]"):
        msg = "input does not look like a wavelength"
        raise ValueError(msg)

    λ_uv = 106.6 * u.nm
    λ_ir = 908.3 * u.nm

    if np.any(λ < λ_uv + 1 * u.nm) or np.any(λ > λ_ir - 1 * u.nm):
        msg = f"this parametrization holds only between {λ_uv+1*u.nm} and {λ_ir-1*u.nm}"
        raise ValueError(msg)

    x = 0.334 + ((0.100 * λ**2) / (λ**2 - λ_uv**2) + (0.008 * λ**2) / (λ**2 - λ_ir**2))

    # solve Clausius-Mossotti
    return (3 + 2 * x) / (3 - x)


def lar_dielectric_constant(
    λ: Quantity, method: ArDielectricMethods = "cern2020"
) -> Quantity:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    See Also
    --------
    .lar_dielectric_constant_bideau_mehu .lar_dielectric_constant_cern2020
    """
    if method == "bideau-mehu":
        return lar_dielectric_constant_bideau_mehu(λ)
    if method == "cern2020":
        return lar_dielectric_constant_cern2020(λ)
    msg = f"Unknown LAr dielectric constant method {method}"
    raise ValueError(msg)


def lar_refractive_index(
    λ: Quantity, method: ArDielectricMethods = "cern2020"
) -> Quantity:
    """Calculate the refractive index of LAr for a given photon wavelength.

    See Also
    --------
    .lar_dielectric_constant


    .. optics-plot:: {'call_x': True}
    """
    return np.sqrt(lar_dielectric_constant(λ, method))


def lar_emission_spectrum(λ: Quantity) -> Quantity:
    """Return the LAr emission spectrum, adapted from [Heindl2010]_.

    .. optics-plot:: {'call_x': True, 'xlim': [116, 141]}
    """
    heindl = readdatafile("lar_emission_heindl2010.dat")

    # sample the measured emission spectrum and avoid the fluctuations below 115 nm.
    return InterpolatingGraph(
        *heindl,
        min_idx=115 * u.nm,
        max_idx=150 * u.nm,
    )(λ)


def lar_fano_factor() -> float:
    """Fano factor.

    Statistical yield fluctuation can be broadened or narrower
    (impurities, fano factor). Value from [Doke1976]_.

    .. optics-const::
    """
    return 0.11


def lar_rayleigh(
    λ: Quantity,
    temperature: Quantity = 90 * u.K,
    method: ArDielectricMethods = "cern2020",
) -> Quantity:
    """Calculate the Rayleigh scattering length using the equations given in [Seidel2002]_.

    This uses the dielectric constant created using the specified method. This calculation
    leads to about 90cm length at 128nm (using cern2020), but keep in mind that the value
    changes drastically out the scintillation peak.

    See Also
    --------
    .lar_dielectric_constant


    .. optics-plot:: {'call_x': True}
    """
    if not temperature.check("[temperature]"):
        msg = "input does not look like a temperature"
        raise ValueError(msg)

    dyne = 1.0e-5 * u.newton
    κ = 2.18e-10 * u.cm**2 / dyne  # LAr isothermal compressibility
    k = 1.380658e-23 * u.joule / u.kelvin  # the Boltzmann constant

    ϵ = lar_dielectric_constant(λ, method)
    assert not np.any(ϵ < 1.00000001)

    inv_l = ((ϵ - 1.0) * (ϵ + 2.0)) ** 2
    inv_l *= κ * temperature * k
    inv_l /= λ**4
    inv_l *= (2 / 3 * np.pi) ** 3

    assert not np.any(inv_l < 1 / (10.0 * u.km))
    assert not np.any(inv_l > 1 / (0.1 * u.nm))

    return (1 / inv_l).to("cm")  # simplify units


def lar_abs_length(λ: Quantity) -> Quantity:
    """Absorption length (not correctly scaled).

    We don't know how the attenuation length actually varies with the wavelength, so here
    we use a custom exponential function connecting the LAr Scintillation peak and the VIS
    range just to avoid a step-like function. Still is a guess.
    This function has to be re-scaled with the intended attenuation length at the VUV
    emission peak.

    .. optics-plot:: {'call_x': True}
    """
    λ = np.maximum(λ, 141 * u.nm)
    absl = 5.976e-12 * np.exp(0.223 * λ.to("nm").m) * u.cm
    return np.minimum(absl, 100000 * u.cm)  # avoid large numbers


def lar_peak_attenuation_length(
    attenuation_method: ArLifetimeMethods | Quantity = "legend200-llama",
) -> Quantity:
    """Attenuation length in the LEGEND-argon, as measured with LLAMA,
    see [Schwarz2024]_ (p. 124).
    """
    if isinstance(attenuation_method, str):
        if attenuation_method == "legend200-llama":
            return 33 * u.cm

        msg = f"unknown attenuation_method {attenuation_method}"
        raise ValueError(msg)

    assert attenuation_method.check("[length]")
    return attenuation_method


def lar_lifetimes(
    triplet_lifetime_method: float | ArLifetimeMethods = "legend200-llama",
) -> ArScintLiftime:
    """Singlet and triplet lifetimes of liquid argon.

    Singlet time from [Hitachi1983]_ and triplet time as measured by LLAMA in LEGEND-200,
    see [Schwarz2024]_ (p. 117).
    """
    triplet = 1 * u.us
    if isinstance(triplet_lifetime_method, str):
        if triplet_lifetime_method == "legend200-llama":
            triplet = 1.16 * u.us
        else:
            msg = f"unknown triplet_lifetime_method {triplet_lifetime_method}"
            raise ValueError(msg)
    else:
        triplet = triplet_lifetime_method * u.us

    return ArScintLiftime(singlet=5.95 * u.ns, triplet=triplet)


def lar_scintillation_params(
    flat_top_yield: Quantity = 31250 / u.MeV,
) -> ScintConfig:
    r"""Scintillation yield (approx. inverse of the mean energy to produce a UV photon).

    This depends on the nature of the impinging particles, the field configuration
    and the quencher impurities. We set here just a reference value that is lower than
    the value provided by [Doke2002]_, that probably does not represent experimental
    reality.

    For flat-top response particles the mean energy to produce a photon is 19.5 eV

    .. math:: Y = 1/(19.5 \mathrm{eV}) = 0.051 \mathrm{eV}^{-1}

    At zero electric field, for not-flat-top particles, the scintillation yield,
    relative to the one of flat top particles is:

    .. math::
        Y_\textrm{e} &= 0.8 Y

        Y_\textrm{alpha} &= 0.7 Y

        Y_\textrm{recoils} &= 0.2\textrm{--}0.4

    Excitation ratio:
    For example, for nuclear recoils it should be 0.75
    nominal value for electrons and gammas: 0.23 (WArP data)

    See Also
    --------
    .lar_fano_factor
    """
    return ScintConfig(
        flat_top=flat_top_yield,
        fano_factor=lar_fano_factor(),
        particles=[
            ScintParticle("electron", yield_factor=0.8, exc_ratio=0.23),
            ScintParticle("alpha", yield_factor=0.7, exc_ratio=1),
            ScintParticle("ion", yield_factor=0.3, exc_ratio=0.75),
        ],
    )


def pyg4_lar_attach_rindex(
    lar_mat, reg, lar_dielectric_method: ArDielectricMethods = "cern2020"
) -> None:
    """Attach the refractive index to the given LAr material instance.

    Parameters
    ----------
    lar_dielectric_method
        Choose which calculation method is used for calculation of the refractive
        index.

    See Also
    --------
    .lar_refractive_index
    .lar_dielectric_constant
    """
    from legendoptics.pyg4utils import pyg4_sample_λ

    λ_full = pyg4_sample_λ(112 * u.nm, 650 * u.nm)
    rindex = lar_refractive_index(λ_full, lar_dielectric_method)
    with u.context("sp"):
        lar_mat.addVecPropertyPint("RINDEX", λ_full.to("eV"), rindex)


def pyg4_lar_attach_attenuation(
    lar_mat,
    reg,
    lar_temperature: Quantity,
    lar_dielectric_method: ArDielectricMethods = "cern2020",
    attenuation_method_or_length: ArLifetimeMethods | Quantity = "legend200-llama",
    rayleigh_enabled_or_length: bool | Quantity = True,
    absorption_enabled_or_length: bool | Quantity = True,
) -> None:
    """Attach all attenuation-related optical properties to the given LAr material instance.

    Parameters
    ----------
    lar_temperature
        liquid phase temperature for rayleigh scattering length calculation.
    lar_dielectric_method
        Choose which calculation method is used for calculation of the dielectric
        function, which is used for deriving the rayleigh scattering length.
    attenuation_method_or_length
        Change the method/measurement used to define the LAr triplet state lifetime.
        If set to a length-Quantity, this value is used directly as attenuation length at
        the scintillation peak.
    rayleigh_enabled_or_length
        If set to a boolean value, it enables or disables the default rayleigh scattering.

        If set to a length-Quantity, the given value will be used as the scattering length at
        the scintillation peak.
    absorption_enabled_or_length
        If set to a boolean value, the default absorption length is used (i.e. it is derived
        from the scattering length and the total attenuation length).

        If set to a length-Quantity, the given value will be used as the absorption length at
        the scintillation peak.

    Notes
    -----
    If all three of rayleigh length, absorption length and attenuation length are set via the function
    parameters, the parameter on total attenuation length will be ignored!

    See Also
    --------
    .lar_rayleigh
    .lar_peak_attenuation_length
    .lar_abs_length
    """
    from legendoptics.pyg4utils import pyg4_sample_λ

    if (
        isinstance(absorption_enabled_or_length, Quantity)
        and isinstance(rayleigh_enabled_or_length, Quantity)
        and isinstance(attenuation_method_or_length, Quantity)
    ):
        log.warning(
            "All three of attenuation, absorption and rayleigh scattering length are constrained manually. The specified attenuation length will be ignored."
        )

    λ_full = pyg4_sample_λ(112 * u.nm, 650 * u.nm)

    # rayleigh scattern is a (theoretically) defined property.
    rayleigh = lar_rayleigh(λ_full, lar_temperature, lar_dielectric_method)
    peak_rayleigh_length = lar_rayleigh(126.8 * u.nm, lar_temperature)
    if isinstance(rayleigh_enabled_or_length, Quantity):
        assert rayleigh_enabled_or_length.check("[length]")
        peak_abs_length = rayleigh_enabled_or_length

    # absorption length and rayleigh add up inversely to the measured attenuation length.
    peak_att_length = lar_peak_attenuation_length(attenuation_method_or_length)

    peak_abs_length = 1 / (1 / peak_att_length - 1 / peak_rayleigh_length)
    if isinstance(absorption_enabled_or_length, Quantity):
        assert absorption_enabled_or_length.check("[length]")
        peak_abs_length = absorption_enabled_or_length

    absl_scale = peak_abs_length / lar_abs_length(126.8 * u.nm)
    abslength = lar_abs_length(λ_full) * absl_scale

    with u.context("sp"):
        if rayleigh_enabled_or_length is not False:
            lar_mat.addVecPropertyPint("RAYLEIGH", λ_full.to("eV"), rayleigh)
        if absorption_enabled_or_length is not False:
            lar_mat.addVecPropertyPint("ABSLENGTH", λ_full.to("eV"), abslength)


def pyg4_lar_attach_scintillation(
    lar_mat,
    reg,
    flat_top_yield: Quantity = 31250 / u.MeV,
    triplet_lifetime_method: float | ArLifetimeMethods = "legend200-llama",
) -> None:
    """Attach all properties for LAr scintillation response to the given LAr material instance.

    Parameters
    ----------
    flat_top_yield_per_mev
        Change the flat-top light yield of the scintillation response. Note that for
        different particle types, the value might be lower (see :func:`lar_scintillation_params`).
    triplet_lifetime_method
        Change the method/measurement used to define the LAr triplet state lifetime.
        If set to a number, this value is used directly as lifetime in µs.

    See Also
    --------
    .lar_scintillation_params
    .lar_emission_spectrum
    .lar_lifetimes
    """
    from legendoptics.pyg4utils import pyg4_def_scint_by_particle_type, pyg4_sample_λ

    λ_peak = pyg4_sample_λ(116 * u.nm, 141 * u.nm)

    # sample the measured emission spectrum.
    scint_em = lar_emission_spectrum(λ_peak)
    # make sure that the scintillation spectrum is zero at the boundaries.
    scint_em[0] = 0
    scint_em[-1] = 0

    with u.context("sp"):
        lar_scint = lar_mat.addVecPropertyPint(
            "SCINTILLATIONCOMPONENT1", λ_peak.to("eV"), scint_em
        )
        lar_mat.addProperty("SCINTILLATIONCOMPONENT2", lar_scint)

    lifetimes = lar_lifetimes(triplet_lifetime_method)
    lar_mat.addConstPropertyPint("SCINTILLATIONTIMECONSTANT1", lifetimes.singlet)
    lar_mat.addConstPropertyPint("SCINTILLATIONTIMECONSTANT2", lifetimes.triplet)

    scint_params = lar_scintillation_params(flat_top_yield)

    # the fano factor is the ratio between variance and mean. Geant4 calculates
    # σ = RESOLUTIONSCALE × √mean, so we have to take the root of the fano factor
    # here to keep it consistent.
    lar_mat.addConstPropertyPint("RESOLUTIONSCALE", np.sqrt(scint_params.fano_factor))

    pyg4_def_scint_by_particle_type(lar_mat, scint_params)


def g4gps_lar_emissions_spectrum(filename: str) -> None:
    """Write a LAr emission energy spectrum for G4GeneralParticleSource.

    See Also
    --------
    .lar_emission_spectrum
    utils.g4gps_write_emission_spectrum
    """
    from legendoptics.pyg4utils import pyg4_sample_λ

    λ_peak = pyg4_sample_λ(116 * u.nm, 141 * u.nm)

    # sample the measured emission spectrum.
    scint_em = lar_emission_spectrum(λ_peak)
    # make sure that the scintillation spectrum is zero at the boundaries.
    scint_em[0] = 0
    scint_em[-1] = 0

    g4gps_write_emission_spectrum(filename, λ_peak, scint_em)
