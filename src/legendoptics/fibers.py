"""
Scinillating fibers BCF91-A from Saint Gobain.

.. [SaintGobainDataSheet] https://www.crystals.saint-gobain.com/sites/hps-mac3-cma-crystals/files/2021-11/Fiber-Product-Sheet.pdf
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity

from legendoptics.utils import InterpolatingGraph, readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def fiber_cladding2_refractive_index() -> float:
    """Refractive index of second fiber cladding material [SaintGobainDataSheet]_.

    .. optics-const::
    """
    return 1.42


def fiber_cladding1_refractive_index() -> float:
    """Refractive index of first fiber cladding material [SaintGobainDataSheet]_.

    .. optics-const::
    """
    return 1.49


def fiber_core_refractive_index() -> float:
    """Refractive index of fiber core material [SaintGobainDataSheet]_.

    .. optics-const::
    """
    return 1.6


def fiber_wls_absorption(
    abs_at_400nm: Quantity = 0.7 * u.mm,
) -> tuple[Quantity, Quantity]:
    r"""[SaintGobainDataSheet]_ reports the absorption spectrum for BCF-91A.

    Knowing that the fibers are 1mm thick one can
    extract the absorption length: starting from the trivial relation:

    :math:`1 - P(E) = \exp(-x/l(E))`

    where :math:`P(E)` is the probability (thus proportional to the absorption spectrum) for a photon
    travelling a distance :math:`x` to be absorbed in the material given the attenuation length
    :math:`l(E)`, one can extract :math:`l(E)` from :math:`P(E)`. By integrating over the thickness of
    the material :math:`L` one obtains:

    :math:`(1 - P(E)) \cdot L = l(E) \cdot (1 - \exp(-L/l(E)))`

    but the problem now is that :math:`l(E)` cannot be extracted analytically (inhomogeneus expression).
    Luigi wrote a Mathematica script that solves it numerically.
    Remember that the units are arbitrary because the original absorption
    spectrum has arbitrary units.

    Measured an absorption length of 0.7 mm at 400 nm, the spectrum has been rescaled by
    that.

    .. optics-plot:: {'yscale': 'log'}
    """
    wvl, absorp = readdatafile("psfibers_wlsabslength.dat")  # arbitrary unit
    assert str(absorp.dimensionality) == "dimensionless"
    # scale factor for absorption lengths (abslength is 0.7mm at 400nm, see above)
    absorp *= abs_at_400nm / InterpolatingGraph(wvl, absorp)(400 * u.nm)
    return wvl, absorp


def fiber_wls_emission() -> tuple[Quantity, Quantity]:
    """[SaintGobainDataSheet]_ reports the emission spectrum for BCF-91A.

    .. optics-plot::
    """
    return readdatafile("psfibers_wlscomponent.dat")


def fiber_wls_timeconstant() -> Quantity:
    """WLS time constant [SaintGobainDataSheet]_.

    .. optics-const::
    """
    return 12 * u.ns


def fiber_absorption_length() -> Quantity:
    """Absorption length of fiber [SaintGobainDataSheet]_. Note this is a macroscopical value for a 1 mm fiber.

    See Also
    --------
    .fiber_absorption_path_length


    .. optics-const::
    """
    return 3.5 * u.m


def fiber_absorption_path_length() -> Quantity:
    """Absorption length of fiber [SaintGobainDataSheet]_, corrected for the geometry of a 1 mm square fiber.

    Multiplied by an empirical factor to account for the prolonged path length inside a square fiber with 1mm side length.

    See Also
    --------
    .fiber_absorption_length


    .. optics-const::
    """
    return fiber_absorption_length() * 1.21


def pyg4_fiber_cladding1_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given fiber cladding 1 material instance.

    See Also
    --------
    .fiber_cladding1_refractive_index
    """
    λ = np.array([650.0, 115.0]) * u.nm
    r = [fiber_cladding1_refractive_index()] * 2
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_fiber_cladding2_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given fiber cladding 2 material instance.

    See Also
    --------
    .fiber_cladding2_refractive_index
    """
    λ = np.array([650.0, 115.0]) * u.nm
    r = [fiber_cladding2_refractive_index()] * 2
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_fiber_core_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given fiber core material instance.

    See Also
    --------
    .fiber_core_refractive_index
    """
    λ = np.array([650.0, 115.0]) * u.nm
    r = [fiber_core_refractive_index()] * 2
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_fiber_core_attach_wls(
    mat,
    reg,
    wls_abs_at_400nm: Quantity = 0.7 * u.mm,
) -> None:
    """Attach wavelength shifting properties to the given material instance.

    See Also
    --------
    .fiber_wls_absorption
    .fiber_wls_emission
    .fiber_wls_timeconstant
    """
    from legendoptics.pyg4utils import pyg4_sample_λ

    λ_full = pyg4_sample_λ(112 * u.nm, 650 * u.nm)
    absorption = InterpolatingGraph(*fiber_wls_absorption())(λ_full)
    emission = InterpolatingGraph(*fiber_wls_emission())(λ_full)
    # make sure that the scintillation spectrum is zero at the boundaries.
    emission[0] = 0
    emission[-1] = 0

    with u.context("sp"):
        mat.addVecPropertyPint("WLSABSLENGTH", λ_full.to("eV"), absorption)
        mat.addVecPropertyPint("WLSCOMPONENT", λ_full.to("eV"), emission)

    mat.addConstPropertyPint("WLSTIMECONSTANT", fiber_wls_timeconstant())


def pyg4_fiber_core_attach_absorption(
    mat, reg, use_geometrical_absorption: bool = True
) -> None:
    """Attach absorption to the given material instance.

    Parameters
    ----------
    use_geometrical_absorption
        switch between the absorption length as specified by the manufacturer and the length
        corrected for the geometry of a 1x1mm fiber.

    See Also
    --------
    .fiber_absorption_path_length
    """
    from legendoptics.pyg4utils import pyg4_sample_λ

    λ_full = pyg4_sample_λ(112 * u.nm, 650 * u.nm)
    length = (
        fiber_absorption_path_length()
        if use_geometrical_absorption
        else fiber_absorption_length()
    )
    absorption = np.array([length.m] * λ_full.shape[0]) * length.u
    with u.context("sp"):
        mat.addVecPropertyPint("ABSLENGTH", λ_full.to("eV"), absorption)
