"""
Liquid Argon (LAr)

.. [Heindl2010] T. Heindl et al. “The scintillation of liquid argon.” In: EPL 91.6 (Sept. 2010).
    https://doi.org/10.1209/0295-5075/91/62002
.. [Doke1976] Doke et al. “Estimation of Fano factors in liquid argon, krypton, xenon and
    xenon-doped liquid argon. NIM 134 (1976)353, https://doi.org/10.1016/0029-554X(76)90292-5
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


def lar_dielectric_constant_bideau_mehu(
    λ: Quantity[float | NDArray],
) -> Quantity[float | NDArray]:
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
    λ: Quantity[float | NDArray],
) -> Quantity[float | NDArray]:
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


def lar_dielectric_constant(
    λ: Quantity[float | NDArray], method: str = "cern2020"
) -> Quantity[float | NDArray]:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    See Also
    --------
    .lar_dielectric_constant_bideau_mehu .lar_dielectric_constant_cern2020
    """
    if method == "bideau-mehu":
        return lar_dielectric_constant_bideau_mehu(λ)
    elif method == "cern2020":
        return lar_dielectric_constant_cern2020(λ)


def lar_refractive_index(
    λ: Quantity[float | NDArray], method: str = "cern2020"
) -> Quantity[float | NDArray]:
    """Calculate the refractive index of LAr for a given photon wavelength.

    See Also
    --------
    .lar_dielectric_constant
    """
    return np.sqrt(lar_dielectric_constant(λ, method))


def lar_emission_spectrum() -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    """Return the LAr emission spectrum, adapted from [Heindl2010]_."""
    return readdatafile("lar_emission_heindl2010.dat")


def lar_fano_factor() -> float:
    """LAr Fano factor

    statistical yield fluctuation can be broadened or narrower
    (impurities, fano factor). Value 0.11 from [Doke1976]_.
    """
    return 0.11


def lar_rayleigh(
    λ: Quantity[float | NDArray], temperature: Quantity[float], method: str = "cern2020"
) -> Quantity[float | NDArray]:
    """Calculate the Rayleigh scattering length using the equations given in
    [Seidel2002]_, but using the dielectric constant created using the specified method.

    This calculation leads to about 90cm length at 128nm (using cern2020), but keep in
    mind that the value changes drastically out the scintillation peak.

    See Also
    --------
    .lar_dielectric_constant
    """

    dyne = 1.0e-5 * u.newton
    κ_T = 2.18e-10 * u.cm ** 2 / dyne  # LAr isothermal compressibility
    k = 1.380658e-23 * u.joule / u.kelvin  # the Boltzmann constant

    ϵ = lar_dielectric_constant(λ, method)
    assert not np.any(ϵ < 1.00000001)

    invL = ((ϵ - 1.0) * (ϵ + 2.0)) ** 2
    invL *= κ_T * temperature * k
    invL /= λ**4
    invL *= (2 / 3 * np.pi) ** 3

    assert not np.any(invL < 1 / (10.0 * u.km)) and not np.any(invL > 1 / (0.1 * u.nm))

    return (1 / invL).to("cm")  # simplify units
