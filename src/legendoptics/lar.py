"""Liquid Argon (LAr).

.. [Heindl2010] T. Heindl et al. “The scintillation of liquid argon.” In: EPL 91.6 (Sept. 2010).
    https://doi.org/10.1209/0295-5075/91/62002
.. [Doke1976] Doke et al. “Estimation of Fano factors in liquid argon, krypton, xenon and
    xenon-doped liquid argon.” NIM 134 (1976)353, https://doi.org/10.1016/0029-554X(76)90292-5
.. [Bideau-Mehu1980] Bideau-Mehu et al. “Measurement of refractive indices of neon, argon, krypton
    and xenon in the 253.7–140.4 nm wavelength range. Dispersion relations and
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
    from liquid argon and xenon.” In: Phys. Rev. B 27 (9 May 1983), pp. 5279–5285,
    https://doi.org/10.1103/PhysRevB.27.5279
.. [Pertoldi2020] L. Pertoldi “Search for new physics with two-neutrino double-beta decay in
    GERDA data.” 2020, https://www.mpi-hd.mpg.de/gerda/public/2020/phd2020_LuigiPertoldi.pdf
"""
from __future__ import annotations

import logging
from typing import NamedTuple

import numpy as np
import pint
from pint import Quantity

from legendoptics.utils import (
    InterpolatingGraph,
    ScintConfig,
    ScintParticle,
    readdatafile,
)

log = logging.getLogger(__name__)
u = pint.get_application_registry()


class ArScintLiftime(NamedTuple):
    singlet: Quantity
    triplet: Quantity


def lar_dielectric_constant_bideau_mehu(
    λ: Quantity,
) -> Quantity:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    From the Bideau-Sellmeier formula [Bideau-Mehu1980]_ in gaseous argon,
    density-corrected for liquid argon.
    """
    if not λ.check("[length]"):
        raise ValueError("input does not look like a wavelength")

    if np.any(λ < 110 * u.nm):
        raise ValueError(f"this parametrization is not meaningful below {110*u.nm}")

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
    """
    if not λ.check("[length]"):
        raise ValueError("input does not look like a wavelength")

    λ_uv = 106.6 * u.nm
    λ_ir = 908.3 * u.nm

    if np.any(λ < λ_uv + 1 * u.nm) or np.any(λ > λ_ir - 1 * u.nm):
        raise ValueError(
            f"this parametrization holds only between {λ_uv+1*u.nm} and {λ_ir-1*u.nm}"
        )

    x = 0.334 + (
        (0.100 * λ**2) / (λ**2 - λ_uv**2)
        + (0.008 * λ**2) / (λ**2 - λ_ir**2)
    )

    # solve Clausius-Mossotti
    return (3 + 2 * x) / (3 - x)


def lar_dielectric_constant(λ: Quantity, method: str = "cern2020") -> Quantity:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    See Also
    --------
    .lar_dielectric_constant_bideau_mehu .lar_dielectric_constant_cern2020
    """
    if method == "bideau-mehu":
        return lar_dielectric_constant_bideau_mehu(λ)
    elif method == "cern2020":
        return lar_dielectric_constant_cern2020(λ)
    raise ValueError(f"Unknown LAr dielectric constant method {method}")


def lar_refractive_index(λ: Quantity, method: str = "cern2020") -> Quantity:
    """Calculate the refractive index of LAr for a given photon wavelength.

    See Also
    --------
    .lar_dielectric_constant
    """
    return np.sqrt(lar_dielectric_constant(λ, method))


def lar_emission_spectrum() -> tuple[Quantity, Quantity]:
    """Return the LAr emission spectrum, adapted from [Heindl2010]_."""
    return readdatafile("lar_emission_heindl2010.dat")


def lar_fano_factor() -> float:
    """Fano factor.

    Statistical yield fluctuation can be broadened or narrower
    (impurities, fano factor). Value 0.11 from [Doke1976]_.
    """
    return 0.11


def lar_rayleigh(
    λ: Quantity, temperature: Quantity, method: str = "cern2020"
) -> Quantity:
    """Calculate the Rayleigh scattering length using the equations given in [Seidel2002]_.

    This uses the dielectric constant created using the specified method. This calculation
    leads to about 90cm length at 128nm (using cern2020), but keep in mind that the value
    changes drastically out the scintillation peak.

    See Also
    --------
    .lar_dielectric_constant
    """
    if not temperature.check("[temperature]"):
        raise ValueError("input does not look like a temperature")

    dyne = 1.0e-5 * u.newton
    κ = 2.18e-10 * u.cm**2 / dyne  # LAr isothermal compressibility
    k = 1.380658e-23 * u.joule / u.kelvin  # the Boltzmann constant

    ϵ = lar_dielectric_constant(λ, method)
    assert not np.any(ϵ < 1.00000001)

    inv_l = ((ϵ - 1.0) * (ϵ + 2.0)) ** 2
    inv_l *= κ * temperature * k
    inv_l /= λ**4
    inv_l *= (2 / 3 * np.pi) ** 3

    assert not np.any(inv_l < 1 / (10.0 * u.km)) and not np.any(
        inv_l > 1 / (0.1 * u.nm)
    )

    return (1 / inv_l).to("cm")  # simplify units


def lar_abs_length(λ: Quantity) -> Quantity:
    """Absorption length (not correctly scaled).

    We don't know how the attenuation length actually varies with the wavelength, so here
    we use a custom exponential function connecting the LAr Scintillation peak and the VIS
    range just to avoid a step-like function. Still is a guess.
    This function has to be re-scaled with the intended attenuation length at the VUV
    emission peak.
    """
    λ = np.maximum(λ, 141 * u.nm)
    absl = 5.976e-12 * np.exp(0.223 * λ.to("nm").m) * u.cm
    return np.minimum(absl, 100000 * u.cm)  # avoid large numbers


def lar_peak_attenuation_length() -> Quantity:
    """Attenuation length in the LEGEND-argon, as measured with LLAMA."""
    return 30 * u.cm


def lar_lifetimes() -> ArScintLiftime:
    """Singlet and triplet lifetimes of liquid argon.

    Singlet time from [Hitachi1983]_ and triplet time as measured in GERDA phase II
    ([Pertoldi2020]_, fig. 2.10).
    """
    return ArScintLiftime(singlet=5.95 * u.ns, triplet=1 * u.us)


def lar_scintillation_params() -> ScintConfig:
    """Scintillation yield (approx. inverse of the mean energy to produce a UV photon).

    This depends on the nature of the impinging particles, the field configuration
    and the quencher impurities. We set here just a reference value from [Doke2002]_,
    that probably does not represent the reality of GERDA/LEGEND:

    for flat top response particles the mean energy to produce a photon is 19.5 eV
    => Y = 1/19.5 = 0.051

    At zero electric field, for not-flat-top particles, the scintillation yield,
    relative to the one of flat top particles is:
    Y_e = 0.8 Y
    Y_alpha = 0.7 Y
    Y_recoils = 0.2-0.4

    Excitation ratio:
    For example, for nuclear recoils it should be 0.75
    nominal value for electrons and gammas: 0.23 (WArP data)
    """
    return ScintConfig(
        flat_top=31250 / u.MeV,
        particles=[
            ScintParticle("electron", yield_factor=0.8, exc_ratio=0.23),
            ScintParticle("alpha", yield_factor=0.7, exc_ratio=1),
            ScintParticle("ion", yield_factor=0.3, exc_ratio=0.75),
        ],
    )


def pyg4_lar_define_opticalproperties(
    lar_mat, lar_temperature: Quantity, reg, lar_dielectric_method="cern2020"
) -> None:
    """Define all liquid argon optical properties on a geant4 material, as defined by this module."""
    from legendoptics.pyg4utils import pyg4_def_scint_by_particle_type, pyg4_sample_λ

    λ_full = pyg4_sample_λ(112 * u.nm, 650 * u.nm)
    λ_peak = pyg4_sample_λ(116 * u.nm, 141 * u.nm)

    # build arrays with properties
    rindex = lar_refractive_index(λ_full, lar_dielectric_method)
    rayleigh = lar_rayleigh(λ_full, lar_temperature, lar_dielectric_method)

    #
    peak_rayleigh_length = lar_rayleigh(126.8 * u.nm, lar_temperature)

    # absorption length and rayleigh add up inversely to the attenuation length.
    peak_abs_length = 1 / (1 / lar_peak_attenuation_length() - 1 / peak_rayleigh_length)
    absl_scale = peak_abs_length / lar_abs_length(126.8 * u.nm)
    abslength = lar_abs_length(λ_full) * absl_scale

    scint_em = InterpolatingGraph(
        *lar_emission_spectrum(),
        min_idx=115 * u.nm,
        max_idx=150 * u.nm,
    )(λ_peak)
    # make sure that the scintillation spectrum is zero at the boundaries.
    scint_em[0] = 0
    scint_em[-1] = 0

    lar_mat.addConstProperty("RESOLUTIONSCALE", lar_fano_factor())

    with u.context("sp"):
        lar_mat.addVecProperty("RINDEX", λ_full.to("eV"), rindex)
        lar_mat.addVecProperty("RAYLEIGH", λ_full.to("eV"), rayleigh)
        lar_mat.addVecProperty("ABSLENGTH", λ_full.to("eV"), abslength)

        lar_scint = lar_mat.addVecProperty(
            "SCINTILLATIONCOMPONENT1", λ_peak.to("eV"), scint_em
        )
        lar_mat.addProperty("SCINTILLATIONCOMPONENT2", lar_scint)

    lar_mat.addConstProperty("SCINTILLATIONTIMECONSTANT1", lar_lifetimes().singlet)
    lar_mat.addConstProperty("SCINTILLATIONTIMECONSTANT2", lar_lifetimes().triplet)

    pyg4_def_scint_by_particle_type(lar_mat, lar_scintillation_params())
