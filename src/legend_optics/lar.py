import logging

import pint
from pint import Quantity

log = logging.getLogger(__name__)

u = pint.get_application_registry()


def lar_dielectric_constant_bideau_mehu(λ: Quantity) -> float:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    From the Bideau-Sellmeier formula [Bideau-Mehu1980]_.

    .. [Bideau-Mehu1980] https://doi.org/10.1016/0022-4073(81)90057-1
    """
    if not λ.check("[length]"):
        raise ValueError("input does not look like a wavelength")

    # wavelength must be > 110.0 nm
    if λ < 110 * u.nm:
        return 1.0e4

    # equation for n-1
    ϵ = 1.2055e-2 * (
        0.2075 * λ**2 / (91.012 * λ**2 - 1 * u.um**2)
        + 0.0415 * λ**2 / (87.892 * λ**2 - 1 * u.um**2)
        + 4.3330 * λ**2 / (214.02 * λ**2 - 1 * u.um**2)
    )
    ϵ *= 2 / 3  # Bideau-Sellmeier -> Clausius-Mossotti
    ϵ *= 1.396 / 1.66e-3  # density correction (Ar gas -> LAr liquid)

    if ϵ < 0.0 or ϵ > 0.999999:
        return 4.0e6
    else:
        return (1 + 2 * ϵ) / (1 - ϵ)  # solve Clausius-Mossotti


def lar_dielectric_constant_cern2020(λ: Quantity) -> float:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    From [Babicz2020]_.

    .. [Babicz2020] https://doi.org/10.1088/1748-0221/15/09/P09009
    """
    if not λ.check("[length]"):
        raise ValueError("input does not look like a wavelength")

    x = 0.334 + (
        (0.100 * λ**2) / (λ**2 - (106.6 * u.nm) ** 2)
        + (0.008 * λ**2) / (λ**2 - (908.3 * u.nm) ** 2)
    )

    return (3 + 2 * x) / (3 - x)  # solve Clausius-Mossotti


def lar_dielectric_constant(λ: Quantity, method: str = "cern2020") -> float:
    """Calculate the dielectric constant of LAr for a given photon wavelength.

    See Also
    --------
    .lar_dielectric_constant_bideau_mehu .lar_dielectric_constant_cern2020
    """
    if method == "bideau-mehu":
        return lar_dielectric_constant_bideau_mehu(λ)
    elif method == "cern2020":
        return lar_dielectric_constant_cern2020(λ)
