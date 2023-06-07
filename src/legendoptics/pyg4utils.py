from __future__ import annotations

import numpy as np
import pint
import logging
import scipy.interpolate
from importlib_resources import files
from numpy.typing import NDArray
from pint import Quantity

import pyg4ometry
import pyg4ometry.geant4 as g4

from .utils import ScintConfig

log = logging.getLogger(__name__)
ureg = pint.get_application_registry().get()


NUMENTRIES = 69
NUMENTRIES_2 = 200 #500


@ureg.with_context('sp')
def pyg4_sample_λ(
    start_lambda: Quantity,
    end_lambda: Quantity,
    sample_count: int = NUMENTRIES_2
) -> Quantity:
    """Sample equally-spaced energies between the two specified wavelengths, and again return the energies
    converted to wavelengths."""
    assert(start_lambda <= end_lambda)

    samples = np.linspace(end_lambda.to('eV'), start_lambda.to('eV'), num=sample_count)
    return samples.to('nm')


def _get_scint_yield_vector(yield_per_mev: Quantity):
    """In Geant4 11.0, ScintillationByParticleType takes some sort of integrated scintiallation yield:

    ScintillationYield = yieldVector->Value(PreStepKineticEnergy)
            - yieldVector->Value(PreStepKineticEnergy - StepEnergyDeposit);

    and ScintillationYield is directly used as MeanNumberOfPhotons. To fulfill this we use a simple
    linear function.
    """
    ye = ureg.Quantity(np.array([1, 10e6]), ureg.eV)
    yv = [f'{(10*ureg.eV*yield_per_mev).to_reduced_units():~}', f'{(10*ureg.MeV*yield_per_mev).to_reduced_units():~}']
    return ye, yv


def _def_scint_particle(mat, particle: str, y: Quantity, yield_factor: float, exc_ratio: float) -> None:
    """Defines a single particle type used by Geant4's ScintillationByParticleType."""
    mat.addVecProperty(particle + 'SCINTILLATIONYIELD', *_get_scint_yield_vector(y * yield_factor))
    mat.addConstProperty(particle + 'SCINTILLATIONYIELD1', exc_ratio)
    mat.addConstProperty(particle + 'SCINTILLATIONYIELD2', 1 - exc_ratio)


def pyg4_def_scint_by_particle_type(mat, scint_cfg: ScintConfig) -> None:
    """Defines a full set of particles for scintillation."""
    for particle in scint_cfg.particles:
        _def_scint_particle(mat, particle.name.upper(), scint_cfg.flat_top, particle.yield_factor, particle.exc_ratio)


@pint.register_unit_format("gdml")
def _gdml_format(unit, registry, **options):
    proc = { u.replace("µ", "u"): e for u,e in unit.items() }
    return pint.formatter(
        proc.items(),
        as_ratio=True,
        single_denominator=False,
        product_fmt="*",
        division_fmt="/",
        power_fmt="{}{}", #TODO: validate only validate power units (e.g. mm2) get through.
        parentheses_fmt="({})",
        **options,
    )


def _patch_g4_pint_unit_support() -> None:
    """code::`pyg4ometry` does currently not support code::`pint` unit that we use extensively here.

    This function overrides some helper functions to make adding material properties nicer.
    The overridden functions also check for some common properties to be used with the correct units.
    """
    def _val_pint_to_gdml(v):
        if not isinstance(v, pint.Quantity):
            return '', v

        base_unit = v.units

        unit = f"{base_unit:~gdml}"
        assert(unit==f"{base_unit:~}".replace(" ", "").replace("µ", "u"))
        log.debug(f"Unit pint->gdml: {unit} - {base_unit}")

        v = v.m_as(base_unit)
        return unit, v

    orig_addVecProperty = g4.WithPropertiesBase.addVecProperty
    def addVecProperty(self, name, e, v):
        vunit, v = _val_pint_to_gdml(v)
        eunit, e = _val_pint_to_gdml(e)

        if name in ["ABSLENGTH", "WLSABSLENGTH", "RAYLEIGH"] and vunit not in ["m", "cm", "mm", "um"]:
            log.warn("Wrong unit %s for property %s", vunit, name)
        if name in ["RINDEX", "WLSCOMPONENT", "REFLECTIVITY"] and vunit != "":
            log.warn("Wrong unit %s for property %s", vunit, name)
        if eunit not in ["", "eV", "keV", "MeV", "GeV", "TeV" "PeV"]:
            log.warn("Wrong energy unit %s", eunit)

        return orig_addVecProperty(self, name, e, v, eunit, vunit)
    g4.WithPropertiesBase.addVecProperty = addVecProperty

    orig_addConstProperty = g4.WithPropertiesBase.addConstProperty
    def addConstProperty(self, name, value):
        vunit, value = _val_pint_to_gdml(value)

        if name in ["SCINTILLATIONYIELD"]:
            log.warn("%s cannot be used with scintillationByParticleType", name)

        return orig_addConstProperty(self, name, value, vunit)
    g4.WithPropertiesBase.addConstProperty = addConstProperty

_patch_g4_pint_unit_support()


__all__ = ['pyg4_sample_λ', 'pyg4_def_scint_by_particle_type']
