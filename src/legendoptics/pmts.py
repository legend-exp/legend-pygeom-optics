"""
PMT components mainly based on [ETEL2010]_.

.. [ETEL2010] ET Enterprises Limited 2010 “200 mm (8") photomultiplier 9354KB series data sheet”, 2010, http://lampes-et-tubes.info/pm/9354KB.pdf

"""

from __future__ import annotations

import logging

import numpy as np
import pint
from pint import Quantity

from legendoptics import store
from legendoptics.utils import readdatafile

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
def pmt_steel_reflectivity() -> float:
    """Reflectivity.

    .. optics-const::
    """
    return 0.9


@store.register_pluggable
def pmt_steel_efficiency() -> float:
    """Efficiency.

    .. optics-const::
    """
    return 1.0


@store.register_pluggable
def pmt_etl9350kb_photocathode_collection_efficiency() -> float:
    """Collection efficiency photocathode -Electron Tubes Limited 9350KB.

    .. optics-const::
    """
    return 0.85  # estimation


@store.register_pluggable
def pmt_r7081_photocathode_collection_efficiency() -> float:
    """Collection efficiency photocathode -Hamamatsu R7081.

    .. optics-const::
    """
    return 0.9  # estimation


@store.register_pluggable
def pmt_etl9350kb_photocathode_efficiency() -> tuple[Quantity, Quantity]:
    """Efficiency.

       The ETL9350KB PMTs.

    .. optics-plot::
    """

    return readdatafile("pmt_etl9350kb_qe.dat")


@store.register_pluggable
def pmt_r7081_photocathode_efficiency() -> tuple[Quantity, Quantity]:
    """Efficiency.

       The R7081 Hamamatsu PMTs.
       https://www.hamamatsu.com/us/en/product/optical-sensors/pmt/pmt_tube-alone/head-on-type/R7081.html

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

    wvl, _ = readdatafile("pmt_qe.dat")

    reflectivity_max = (
        (1 - pmt_borosilicate_refractive_index())
        / (1 + pmt_borosilicate_refractive_index())
    ) ** 2
    reflectivity = np.full_like(wvl, reflectivity_max - 0.01)
    return wvl, reflectivity


def pyg4_pmt_attach_acryl_rindex(mat, reg) -> None:
    """Attach the refractive index to the given acryl material instance of the PMT cap.

    See Also
    --------
    .pmt_acryl_refractive_index
    """
    energy = np.array([1.0, 6.0]) * u.eV
    r = [pmt_acryl_refractive_index()] * 2

    mat.addVecPropertyPint("RINDEX", energy, r)


def pyg4_pmt_attach_acryl_absorption_length(mat, reg) -> None:
    """Attach the absorption length to the given acryl material instance of the PMT cap.

    See Also
    --------
    .pmt_acryl_absorption_length
    """

    energy, absorpt = pmt_acryl_absorption_length()

    mat.addVecPropertyPint("ABSLENGTH", energy, absorpt)


def pyg4_pmt_attach_air_rindex(mat, reg) -> None:
    """Attach the refractive index to the given air material instance of the PMT cap.

    See Also
    --------
    .pmt_air_refractive_index
    """
    energy = np.array([1.0, 6.0]) * u.eV
    r = [pmt_air_refractive_index()] * 2

    mat.addVecPropertyPint("RINDEX", energy, r)


def pyg4_pmt_attach_air_absorption_length(mat, reg) -> None:
    """Attach the absorption length to the given air material instance of the PMT cap.

    See Also
    --------
    .pmt_air_absorption_length
    """

    energy = np.array([1.0, 6.0]) * u.eV
    absorpt = np.full_like(energy, pmt_air_absorption_length())

    mat.addVecPropertyPint("ABSLENGTH", energy, absorpt)


def pyg4_pmt_attach_borosilicate_rindex(mat, reg) -> None:
    """Attach the refractive index to the given borosilicate material instance of the PMT cap.

    See Also
    --------
    .pmt_borosilicate_refractive_index
    """
    energy = np.array([1.0, 6.0]) * u.eV
    r = [pmt_borosilicate_refractive_index()] * 2

    mat.addVecPropertyPint("RINDEX", energy, r)


def pyg4_pmt_attach_borosilicate_absorption_length(mat, reg) -> None:
    """Attach the absorption length to the given borosilicate material instance of the PMT cap.

    See Also
    --------
    .pmt_borosilicate_absorption_length
    """

    energy, absorpt = pmt_borosilicate_absorption_length()

    mat.addVecPropertyPint("ABSLENGTH", energy, absorpt)


def pyg4_pmt_attach_steel_reflectivity(mat, reg) -> None:
    """Attach the reflectivity to the given PMT steel material instance.

    See Also
    --------
    .pmt_steel_reflectivity
    """
    energy = np.array([1.0, 6.0]) * u.eV
    refl = [pmt_steel_reflectivity()] * 2

    mat.addVecPropertyPint("REFLECTIVITY", energy, refl)


def pyg4_pmt_attach_steel_efficiency(mat, reg) -> None:
    """Attach the efficiency to the given PMT steel material instance.

    See Also
    --------
    .pmt_steel_efficiency
    """

    energy = np.array([1.0, 6.0]) * u.eV
    eff = [pmt_steel_efficiency()] * 2

    mat.addVecPropertyPint("EFFICIENCY", energy, eff)


def pyg4_pmt_attach_photocathode_reflectivity(mat, reg) -> None:
    """Attach the reflectivity to the given PMT photocathode material instance.

    See Also
    --------
    .pmt_photocathode_reflectivity
    """
    wvl, refl = pmt_photocathode_reflectivity()

    with u.context("sp"):
        mat.addVecPropertyPint("REFLECTIVITY", wvl.to("eV"), refl)


def pyg4_pmt_attach_photocathode_efficiency(mat, reg, name="etl9350") -> None:
    """Attach the efficiency to the given PMT photocathode material instance.

    See Also
    --------
    .pmt_photocathode_efficiency
    .pmt_photocathode_collection_efficiency
    """

    if "etl9350" in name.lower() or "gerda" in name.lower():
        wvl, pmt_qe = pmt_etl9350kb_photocathode_efficiency()
        pmt_efficiency = (
            pmt_qe / 100 * pmt_etl9350kb_photocathode_collection_efficiency()
        )
    elif "r7081" in name.lower() or "l1000" in name.lower():
        wvl, pmt_qe = pmt_r7081_photocathode_efficiency()
        pmt_efficiency = pmt_qe / 100 * pmt_r7081_photocathode_efficiency()
    else:
        msg = "PMT name not known. There exists only R7081 or ETL9350 data."
        raise ValueError(msg)

    with u.context("sp"):
        mat.addVecPropertyPint("EFFICIENCY", wvl.to("eV"), pmt_efficiency)
