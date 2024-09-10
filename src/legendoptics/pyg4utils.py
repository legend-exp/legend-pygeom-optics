from __future__ import annotations

import logging

import numpy as np
import pint
import pyg4ometry.gdml.Units as gdml_u
import pyg4ometry.geant4 as g4
from pint import Quantity

from .scintillate import ScintConfig

log = logging.getLogger(__name__)
ureg = pint.get_application_registry().get()


@ureg.with_context("sp")
def pyg4_sample_λ(
    start_lambda: Quantity, end_lambda: Quantity, sample_count: int = 200
) -> Quantity:
    """Sample equally-spaced energies between the two specified wavelengths."""
    assert start_lambda <= end_lambda

    samples = np.linspace(end_lambda.to("eV"), start_lambda.to("eV"), num=sample_count)
    return samples.to("nm")


def _get_scint_yield_vector(yield_per_mev: Quantity):
    """In Geant4 11.0+, ScintillationByParticleType takes some sort of integrated scintillation yield.

    To fulfill this we use a simple linear function.
    """
    ye = ureg.Quantity(np.array([1, 10e6]), ureg.eV)
    yv = [(e * yield_per_mev).to_reduced_units().m for e in ye]
    return ye, yv


def _def_scint_particle(
    mat, particle: str, y: Quantity, yield_factor: float, exc_ratio: float | None
) -> None:
    """Define a single particle type used by Geant4's ScintillationByParticleType."""
    if not y.check("1/[energy]"):
        msg = "Scintillation yield must have dimensionality 1/energy"
        raise ValueError(msg)

    mat.addVecPropertyPint(
        particle + "SCINTILLATIONYIELD", *_get_scint_yield_vector(y * yield_factor)
    )
    if exc_ratio is None:
        mat.addConstProperty(particle + "SCINTILLATIONYIELD1", 1)
    else:
        mat.addConstProperty(particle + "SCINTILLATIONYIELD1", exc_ratio)
        mat.addConstProperty(particle + "SCINTILLATIONYIELD2", 1 - exc_ratio)


def pyg4_def_scint_by_particle_type(mat, scint_cfg: ScintConfig) -> None:
    """Define a full set of particles for scintillation."""
    for particle in scint_cfg.particles:
        if not particle.valid_geant_particle():
            msg = (
                f"Unknown particle type {particle.name} for ScintillationByParticleType"
            )
            raise ValueError(msg)

        _def_scint_particle(
            mat,
            particle.name.upper(),
            scint_cfg.flat_top,
            particle.yield_factor,
            particle.exc_ratio,
        )


def _gdml_unit(u: str) -> str:
    # Only as of Geant4 11.1.0, `um` and `nm` are supported.
    u = u.replace("µ", "u")
    if u == "nm":
        u = "nanometer"
    if u == "um":
        u = "micrometer"
    return u


@pint.register_unit_format("gdml")
def _gdml_format(unit, registry, **options):
    proc = {_gdml_unit(u): e for u, e in unit.items()}
    return pint.formatting.formatter(
        proc.items(),
        as_ratio=True,
        single_denominator=False,
        product_fmt="*",
        division_fmt="/",
        power_fmt="{}{}",  # TODO: validate only GDML-valid power units (e.g. mm2) get through.
        parentheses_fmt="({})",
        **options,
    )


def _patch_g4_pint_unit_support() -> None:
    """:py:mod:`pyg4ometry` does currently not support code::`pint` unit that we use extensively here.

    This function adds some helper functions to :py:class:`pyg4ometry.geant4.Material` and
    :py:class:`pyg4ometry.geant4.solid.OpticalSurface` in order to make adding material properties nicer.
    The new functions also check for some common properties to be used with the correct units.
    """

    def _val_pint_to_gdml(v):
        if not isinstance(v, pint.Quantity):
            return "", v

        base_unit = v.units

        unit = f"{base_unit:~gdml}"
        # assert unit == f"{base_unit:~}".replace(" ", "").replace("µ", "u")
        assert "dimensionless" not in unit

        msg = f"Unit pint->gdml: {unit} - {base_unit}"
        log.debug(msg)

        v = v.m_as(base_unit)
        return unit, v

    # Only as of Geant4 11.1.0, `um` and `nm` are supported.
    length_u = ["km", "m", "cm", "mm", "nanometer", "micrometer"]
    dimless_props = [
        "RINDEX",
        "WLSCOMPONENT",
        "REFLECTIVITY",
        "REALRINDEX",
        "IMAGINARYRINDEX",
        "SPECULARLOBECONSTANT",
        "SPECULARSPIKECONSTANT",
        "BACKSCATTERCONSTANT",
        "EFFICIENCY",
    ]
    length_props = ["ABSLENGTH", "WLSABSLENGTH", "RAYLEIGH"]

    def addVecPropertyPint(self, name, e, v):
        vunit, v = _val_pint_to_gdml(v)
        eunit, e = _val_pint_to_gdml(e)
        v = np.array(v)
        e = np.array(e)
        # assert that we have only numeric data after this:
        assert e.dtype.kind in "uif"
        assert v.dtype.kind in "uif"

        if name in length_props and vunit not in length_u:
            log.warning("Wrong unit %s for property %s", vunit, name)
        if name in dimless_props and vunit != "":
            log.warning("Wrong unit %s for property %s", vunit, name)
        if eunit not in ["", "eV", "keV", "MeV", "GeV", "TeV", "PeV"]:
            log.warning("Wrong energy unit %s", eunit)

        # reorder the values to be in ascending energy order.
        sort = np.argsort(e)
        e = e[sort]
        v = v[sort]

        return g4.Material.addVecProperty(self, name, e, v, eunit, vunit)

    g4.Material.addVecPropertyPint = addVecPropertyPint
    g4.solid.OpticalSurface.addVecPropertyPint = addVecPropertyPint

    def addConstPropertyPint(self, name, value):
        vunit, value = _val_pint_to_gdml(value)

        if name in ["SCINTILLATIONYIELD"]:
            log.warning("%s cannot be used with scintillationByParticleType", name)

        return g4.Material.addConstProperty(self, name, value, vunit)

    g4.Material.addConstPropertyPint = addConstPropertyPint
    g4.solid.OpticalSurface.addConstPropertyPint = addConstPropertyPint

    # those units are supported by G4, but not by pyg4ometry...
    # the first check is ugly, but necessary: we have a mock import of gdml_u for the docs build.
    if isinstance(gdml_u.units, dict):
        if "nanometer" not in gdml_u.units:
            gdml_u.units["nanometer"] = gdml_u.units["nm"]
        if "micrometer" not in gdml_u.units:
            gdml_u.units["micrometer"] = gdml_u.units["um"]


_patch_g4_pint_unit_support()


__all__ = ["pyg4_sample_λ", "pyg4_def_scint_by_particle_type"]
