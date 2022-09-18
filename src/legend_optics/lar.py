from __future__ import annotations

import logging

import numpy as np
import pint
from numpy.typing import NDArray
from pint import Quantity

from legend_optics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def lar_dielectric_constant_bideau_mehu(λ: Quantity[float | NDArray]) -> float:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    From the Bideau-Sellmeier formula [Bideau-Mehu1980]_.

    .. [Bideau-Mehu1980] https://doi.org/10.1016/0022-4073(81)90057-1
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


def lar_dielectric_constant_cern2020(λ: Quantity[float | NDArray]) -> float:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    From [Babicz2020]_.

    .. [Babicz2020] https://doi.org/10.1088/1748-0221/15/09/P09009
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
) -> float:
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
) -> float:
    """Calculate the refractive index of LAr for a given photon wavelength.

    See Also
    --------
    .lar_dielectric_constant
    """
    return np.sqrt(lar_dielectric_constant(λ, method))


def lar_emission_spectrum() -> tuple[Quantity[NDArray], Quantity[NDArray]]:
    """Return the LAr emission spectrum.

    Adapted from [Heindl2010]_.

    .. [Heindl2010] https://doi.org/10.1209/0295-5075/91/62002
    """
    return readdatafile("lar_emission_heindl2010.dat")
