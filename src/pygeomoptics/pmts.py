"""
PMT components incorporating data for different PMT models.

Common parameters applicable to both models are included where relevant, like refractive index of borosilicate glass.
Currently the photocathode efficiencies for the ETL9354KB [ETEL2010]_ and R7081 [HAMAMATSU2019]_ PMT models are included.

.. [ETEL2010] ET Enterprises Limited 2010 “200 mm (8") photomultiplier 9354KB series data sheet”, 2010, http://lampes-et-tubes.info/pm/9354KB.pdf
.. [HAMAMATSU2019] Hammamatsu Photonics 2019 "Large Area PMT data sheet"
    https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/etd/LARGE_AREA_PMT_TPMH1376E.pdf
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal

import numpy as np
import pint
from pint import Quantity

if TYPE_CHECKING:
    import pyg4ometry.geant4 as g4

from pygeomoptics import store
from pygeomoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def pmt_acryl_refractive_index() -> float:
    """Refractive index.

    .. optics-const::
    """
    return 1.489


@store.register_pluggable
def pmt_acryl_absorption_length() -> tuple[Quantity, Quantity]:
    """Absorption length.

    .. optics-plot::
    """
    energy = np.array([1.0, 6.0]) * u.eV
    absorp_length = np.array([2.5, 3.5]) * u.m  # estimation
    return energy, absorp_length


@store.register_pluggable
def pmt_air_refractive_index() -> float:
    """Refractive index.

    .. optics-const::
    """
    return 1.0


@store.register_pluggable
def pmt_air_absorption_length() -> Quantity:
    """Absorption length.

    .. optics-const::
    """
    return 100 * u.m


@store.register_pluggable
def pmt_borosilicate_refractive_index() -> float:
    """Refractive index.

    .. optics-const::
    """
    return 1.49


@store.register_pluggable
def pmt_borosilicate_absorption_length() -> tuple[Quantity, Quantity]:
    """Absorption length.

    .. optics-plot::
    """
    energy = np.array([1.0, 6.0]) * u.eV
    absorp_length = np.array([2.0, 3.0]) * u.m  # estimation
    return energy, absorp_length


@store.register_pluggable
def pmt_steel_reflectivity() -> tuple[Quantity, Quantity]:
    """Reflectivity. Modeled after [STEEL1982]_.

    .. optics-plot::
    """
    λ = np.array([200, 300, 400, 600, 800]) * u.nm
    refl = np.array([0.35, 0.45, 0.55, 0.58, 0.60])
    return λ, refl


@store.register_pluggable
def pmt_steel_efficiency() -> float:
    """Efficiency.

    .. deprecated:: 0.17

            steel should not have a detection efficiency.
    """
    return 0.0


@store.register_pluggable
def pmt_etl9354kb_photocathode_collection_efficiency() -> float:
    """Collection efficiency photocathode for ETL9354KB.

    .. optics-const::
    """
    return 0.85  # estimation


@store.register_pluggable
def pmt_r7081_photocathode_collection_efficiency() -> float:
    """Collection efficiency photocathode for Hamamatsu R7081.

    .. optics-const::
    """
    return 0.9  # estimation


@store.register_pluggable
def pmt_etl9354kb_photocathode_efficiency() -> tuple[Quantity, Quantity]:
    """Efficiency for ETL9354KB.

    .. optics-plot::
    """

    return readdatafile("pmt_etl9354kb_qe.dat")


@store.register_pluggable
def pmt_r7081_photocathode_efficiency() -> tuple[Quantity, Quantity]:
    """Efficiency for Hamamatsu R7081.

    .. optics-plot::
    """

    return readdatafile("pmt_r7081_qe.dat")


@store.register_pluggable
def pmt_photocathode_reflectivity() -> tuple[Quantity, Quantity]:
    """Efficiency.

    .. optics-plot::

    See Also
    --------
    .borosilicate_refractive_index
    """

    λ = np.array([270, 700]) * u.nm

    reflectivity_max = (
        (1 - pmt_borosilicate_refractive_index())
        / (1 + pmt_borosilicate_refractive_index())
    ) ** 2
    reflectivity = np.full_like(λ, reflectivity_max - 0.01)
    return λ, reflectivity


def pyg4_pmt_attach_acryl_rindex(mat: g4.Material, reg: g4.Registry) -> None:
    """Attach the refractive index to the given acryl material instance of the PMT cap.

    See Also
    --------
    .pmt_acryl_refractive_index
    """
    energy = np.array([1.0, 6.0]) * u.eV
    r = [pmt_acryl_refractive_index()] * 2

    mat.addVecPropertyPint("RINDEX", energy, r)


def pyg4_pmt_attach_acryl_absorption_length(mat: g4.Material, reg: g4.Registry) -> None:
    """Attach the absorption length to the given acryl material instance of the PMT cap.

    See Also
    --------
    .pmt_acryl_absorption_length
    """

    energy, absorpt = pmt_acryl_absorption_length()

    mat.addVecPropertyPint("ABSLENGTH", energy, absorpt)


def pyg4_pmt_attach_air_rindex(mat: g4.Material, reg: g4.Registry) -> None:
    """Attach the refractive index to the given air material instance of the PMT cap.

    See Also
    --------
    .pmt_air_refractive_index
    """
    energy = np.array([1.0, 6.0]) * u.eV
    r = [pmt_air_refractive_index()] * 2

    mat.addVecPropertyPint("RINDEX", energy, r)


def pyg4_pmt_attach_air_absorption_length(mat: g4.Material, reg: g4.Registry) -> None:
    """Attach the absorption length to the given air material instance of the PMT cap.

    See Also
    --------
    .pmt_air_absorption_length
    """

    energy = np.array([1.0, 6.0]) * u.eV
    absorpt = np.full_like(energy, pmt_air_absorption_length())

    mat.addVecPropertyPint("ABSLENGTH", energy, absorpt)


def pyg4_pmt_attach_borosilicate_rindex(mat: g4.Material, reg: g4.Registry) -> None:
    """Attach the refractive index to the given borosilicate material instance of the PMT cap.

    See Also
    --------
    .pmt_borosilicate_refractive_index
    """
    energy = np.array([1.0, 6.0]) * u.eV
    r = [pmt_borosilicate_refractive_index()] * 2

    mat.addVecPropertyPint("RINDEX", energy, r)


def pyg4_pmt_attach_borosilicate_absorption_length(
    mat: g4.Material, reg: g4.Registry
) -> None:
    """Attach the absorption length to the given borosilicate material instance of the PMT cap.

    See Also
    --------
    .pmt_borosilicate_absorption_length
    """

    energy, absorpt = pmt_borosilicate_absorption_length()

    mat.addVecPropertyPint("ABSLENGTH", energy, absorpt)


def pyg4_pmt_attach_steel_reflectivity(mat: g4.Material, reg: g4.Registry) -> None:
    """Attach the reflectivity to the given PMT steel material instance.

    See Also
    --------
    .pmt_steel_reflectivity
    """
    λ, refl = pmt_steel_reflectivity()

    with u.context("sp"):
        mat.addVecPropertyPint("REFLECTIVITY", λ.to("eV"), refl)


def pyg4_pmt_attach_steel_efficiency(mat: g4.Material, reg: g4.Registry) -> None:
    """Attach the efficiency to the given PMT steel material instance.

    .. deprecated:: 0.17

            steel should not have a detection efficiency.
    """

    energy = np.array([1.0, 6.0]) * u.eV
    eff = [pmt_steel_efficiency()] * 2

    mat.addVecPropertyPint("EFFICIENCY", energy, eff)


def pyg4_pmt_attach_photocathode_reflectivity(
    mat: g4.Material, reg: g4.Registry
) -> None:
    """Attach the reflectivity to the given PMT photocathode material instance.

    See Also
    --------
    .pmt_photocathode_reflectivity
    """
    λ, refl = pmt_photocathode_reflectivity()

    with u.context("sp"):
        mat.addVecPropertyPint("REFLECTIVITY", λ.to("eV"), refl)


def pyg4_pmt_attach_photocathode_efficiency(
    mat: g4.Material,
    reg: g4.Registry,
    name: Literal["etl9354", "gerda", "r7081", "l1000"] = "etl9354",
) -> None:
    """Attach the efficiency to the given PMT photocathode material instance.

    See Also
    --------
    .pmt_photocathode_efficiency
    .pmt_photocathode_collection_efficiency
    """

    if name in {"etl9354", "gerda"}:
        λ, pmt_qe = pmt_etl9354kb_photocathode_efficiency()
        pmt_efficiency = (
            pmt_qe / 100 * pmt_etl9354kb_photocathode_collection_efficiency()
        )
    elif name in {"r7081", "l1000"}:
        λ, pmt_qe = pmt_r7081_photocathode_efficiency()
        pmt_efficiency = pmt_qe / 100 * pmt_r7081_photocathode_collection_efficiency()
    else:
        msg = f"PMT name {name} not known. There exists only r7081 or etl9354 data."
        raise ValueError(msg)

    with u.context("sp"):
        mat.addVecPropertyPint("EFFICIENCY", λ.to("eV"), pmt_efficiency)
