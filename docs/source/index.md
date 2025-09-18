# Welcome to legend-geom-optics's documentation!

This package contains a collection of optical properties of materials used in
the [LEGEND](https://legend-exp.org/) experiment.

As a common interface, each optical property gets its own defining function in
the material's module. Those functions can be used directly to just retrieve the
value(s) of the property. Most property definitions contain unit information via
the {mod}`pint` package. For a full list of defined properties see the API reference.

To ease the use in [Geant4](https://geant4.web.cern.ch/)-based simulations,
every module defines one or more functions prefixed with `pyg4_attach_`. Those
functions are to be used with
{mod}`pyg4ometry` and will
attach the listed properties to a material or surface instance.

## Features

- High-level helpers to attach optical properties to Geant4/pyg4ometry materials
  and optical surfaces (RINDEX, ABSLENGTH, RAYLEIGH, WLS*, SCINTILLATION*).
- Pluggable property store to override or swap implementations at runtime
  without forking your code.
- Unit-safe data handling using pint, including wavelength/energy conversion
  contexts for correct conversions.
- Data-driven spectra with a simple file format and interpolation utilities.
- Ready-made properties for many materials (LAr, PEN, TPB, scintillating fibers,
  reflectors like Tyvek/Tetratex/VM2000, metals/semiconductors/glass, water,
  nylon, Ultem, PMTs).
- CLI utility to write G4GeneralParticleSource emission spectra from built-in
  WLS/scintillation spectra.
- Sphinx plotting helpers to render optical property plots in the docs.

See the [User Guide](user_guide) for a short walkthrough of these features.


## Usage example

This example demonstrates how to use the {mod}`legendoptics.lar` submodule to add
important optical properties to a {class}`pyg4ometry.geant4.Material` instance using
{func}`legendoptics.lar.pyg4_lar_attach_rindex`, {func}`legendoptics.lar.pyg4_lar_attach_attenuation`,
and {func}`legendoptics.lar.pyg4_lar_attach_scintillation`:

```python
import legendoptics.lar
import pint
import pyg4ometry.geant4 as g4

g4_registry = g4.Registry()

_liquidargon = g4.Material(
    name="LiquidArgon",
    state="liquid",
    # [...]
    registry=g4_registry,
)

u = pint.get_application_registry().get()
legendoptics.lar.pyg4_lar_attach_rindex(
    _liquidargon,
    g4_registry,
)
legendoptics.lar.pyg4_lar_attach_attenuation(
    _liquidargon,
    g4_registry,
)
legendoptics.lar.pyg4_lar_attach_scintillation(
    _liquidargon,
    g4_registry,
    triplet_lifetime_method="legend200-llama",
)
```

Each property can be used separately:

```python
import legendoptics.lar

lifetimes = legendoptics.lar.lar_lifetimes(triplet_lifetime_method="legend200-llama")
print(lifetimes.triplet)
```
