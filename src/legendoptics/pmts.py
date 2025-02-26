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
    """Absorption length (estimation).

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
def pmt_photocathode_collection_efficiency() -> float:
    """Collection efficiency photocathode.

    .. optics-const::
    """
    return 0.85


@store.register_pluggable
def pmt_photocathode_efficiency() -> tuple[Quantity, Quantity]:
    """Efficiency.

    .. optics-plot::
    """

    return readdatafile("pmt_qe.dat")


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
    reflectivity = [reflectivity_max - 0.01] * len(wvl)
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


def pyg4_pmt_attach_photocathode_efficiency(mat, reg) -> None:
    """Attach the efficiency to the given PMT photocathode material instance.

    See Also
    --------
    .pmt_photocathode_efficiency
    .pmt_photocathode_collection_efficiency
    """

    wvl, pmt_qe = pmt_photocathode_efficiency()
    pmt_efficiency = pmt_qe / 100 * pmt_photocathode_collection_efficiency()

    with u.context("sp"):
        mat.addVecPropertyPint("EFFICIENCY", wvl.to("eV"), pmt_efficiency)
