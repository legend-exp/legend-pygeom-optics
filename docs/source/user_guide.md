# User Guide

This package provides reusable optical material property definitions for fast
integration with Geant4 geometries via pyg4ometry. It includes:

- Ready-to-use functions to attach optical properties (wavelength-dependent and
  constant) to materials and optical surfaces.
- A pluggable store mechanism to override or extend material properties without
  forking the package.
- Utilities to read spectra from data files and interpolate them with physical
  units using pint.

See concrete geometry integrations and usage examples in:

- LEGEND-200:
  [legend-exp/legend-pygeom-l200](https://github.com/legend-exp/legend-pygeom-l200)
- LEGEND-1000:
  [legend-exp/legend-pygeom-l1000](https://github.com/legend-exp/legend-pygeom-l1000)

These repositories demonstrate how this package is used to populate material
properties in full detector geometries.

## High-level API: attaching optical properties

Most materials expose convenience functions with the prefix
`pyg4_<material>_attach_*`. These functions add properties to a Geant4 material
(or optical surface) using pint-aware helpers that ensure correct units and
sorting. Typical property names include:

- `RINDEX`: refractive index (dimensionless)
- `ABSLENGTH`: absorption length (length)
- `RAYLEIGH`: Rayleigh scattering length (length)
- `REFLECTIVITY`, `EFFICIENCY`: surface/optical parameters (dimensionless)
- `WLSABSLENGTH`, `WLSCOMPONENT`, `WLSTIMECONSTANT`: wavelength-shifting
  properties
- `SCINTILLATIONCOMPONENT*`, `SCINTILLATIONTIMECONSTANT*`, `RESOLUTIONSCALE`:
  scintillation models

Examples of high-level attachers:

- Liquid argon (LAr): {func}`legendoptics.lar.pyg4_lar_attach_rindex`,
  {func}`legendoptics.lar.pyg4_lar_attach_attenuation`,
  {func}`legendoptics.lar.pyg4_lar_attach_scintillation`
- PEN: {func}`legendoptics.pen.pyg4_pen_attach_rindex`,
  {func}`legendoptics.pen.pyg4_pen_attach_attenuation`,
  {func}`legendoptics.pen.pyg4_pen_attach_wls`,
  {func}`legendoptics.pen.pyg4_pen_attach_scintillation`
- TPB: {func}`legendoptics.tpb.pyg4_tpb_attach_rindex`,
  {func}`legendoptics.tpb.pyg4_tpb_attach_wls`
- Fibers: {func}`legendoptics.fibers.pyg4_fiber_core_attach_rindex`,
  {func}`legendoptics.fibers.pyg4_fiber_core_attach_absorption`,
  {func}`legendoptics.fibers.pyg4_fiber_core_attach_wls`,
  {func}`legendoptics.fibers.pyg4_fiber_core_attach_scintillation` (plus
  cladding 1/2 rindex attachers)
- Reflectors: {func}`legendoptics.tetratex.pyg4_tetratex_attach_reflectivity`,
  {func}`legendoptics.tyvek.pyg4_tyvek_attach_reflectivity`,
  {func}`legendoptics.vm2000.pyg4_vm2000_attach_wls` (see module for more:
  rindex, reflectivity, border params)
- Substrates/metals:
  {func}`legendoptics.copper.pyg4_copper_attach_reflectivity`,
  {func}`legendoptics.germanium.pyg4_germanium_attach_reflectivity`
- Semiconductors/glass:
  {func}`legendoptics.silicon.pyg4_silicon_attach_complex_rindex`,
  {func}`legendoptics.silica.pyg4_silica_attach_rindex`
- Others: {mod}`legendoptics.nylon`, {mod}`legendoptics.ultem`,
  {mod}`legendoptics.water`, {mod}`legendoptics.pmts`

Minimal usage outline (pseudo-code):

```python
# Minimal outline with pyg4ometry (exact construction details may vary):
import pint
import pyg4ometry.geant4 as g4

from legendoptics.lar import (
    pyg4_lar_attach_rindex,
    pyg4_lar_attach_attenuation,
    pyg4_lar_attach_scintillation,
)

u = pint.get_application_registry().get()
reg = g4.Registry()

# Create your material in pyg4ometry here (example only; adjust to your setup)
lar_mat = g4.Material(name="LAr", density=1.396, registry=reg)  # density in g/cm**3

# Attach optical properties
pyg4_lar_attach_rindex(lar_mat, reg)
pyg4_lar_attach_attenuation(lar_mat, reg, lar_temperature=88.8 * u.K)
pyg4_lar_attach_scintillation(lar_mat, reg)
```

Under the hood, the helpers use {mod}`legendoptics.pyg4utils` to patch
{mod}`pyg4ometry` with pint-aware methods:

- `Material.addVecPropertyPint(name, energy, values)`
- `Material.addConstPropertyPint(name, value)`

These ensure proper unit handling and ascending-energy ordering.

## Modifying properties without forking (pluggable store)

Many property functions are decorated with
{func}`legendoptics.store.register_pluggable`. This makes them dynamically
replaceable at runtime.

Common operations:

```python
from legendoptics import store
from legendoptics.pen import pen_refractive_index


# Replace a property implementation
def my_pen_rindex() -> float:
    return 1.52  # new constant refractive index


pen_refractive_index.replace_implementation(my_pen_rindex)

# Query which functions were replaced
assert "pen_refractive_index" in store.get_replaced()

# Reset just this function
pen_refractive_index.reset_implementation()

# Or reset everything back to defaults
store.reset_all_to_original()
```

Note: Any imported function that is decorated with `@store.register_pluggable`
can be replaced. The wrapper preserves the original function and exposes:

- `wrap.replace_implementation(new_impl)`
- `wrap.reset_implementation()`
- `wrap.is_original()`
- `wrap.original_impl()`

## Adding a new material

To integrate a new material in the same style:

1. Provide data-backed property functions

- Use {func}`legendoptics.utils.readdatafile` to read spectra with units from a
  data file included in a Python package (default `legendoptics.data`).
- Use {class}`legendoptics.utils.InterpolatingGraph` for interpolation with
  correct unit handling.
- Decorate property functions you want to be overridable with
  `@store.register_pluggable`.

Example:

```python
import numpy as np
import pint
from legendoptics import store
from legendoptics.utils import readdatafile, InterpolatingGraph

u = pint.get_application_registry()


@store.register_pluggable
def mymat_refractive_index() -> float:
    return 1.37


@store.register_pluggable
def mymat_absorption() -> tuple[u.Quantity, u.Quantity]:
    λ, L = readdatafile("mymat_absorption.dat")  # "# unit1 unit2" header expected
    return λ, L
```

2. Implement pyg4 attachers

- Sample wavelengths with {func}`legendoptics.pyg4utils.pyg4_sample_λ` where
  needed.
- Convert wavelength to energy inside a pint context and attach via
  `addVecPropertyPint` and `addConstPropertyPint`.

Example:

```python
import numpy as np
from legendoptics.pyg4utils import pyg4_sample_λ
import pint

u = pint.get_application_registry()


def pyg4_mymat_attach_rindex(mat, reg):
    λ = np.array([650.0, 115.0]) * u.nm
    r = [mymat_refractive_index()] * 2
    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_mymat_attach_absorption(mat, reg):
    λ, L = mymat_absorption()
    with u.context("sp"):
        mat.addVecPropertyPint("ABSLENGTH", λ.to("eV"), L)
```

3. Optional: WLS and scintillation

- Follow patterns in `pen.py`, `tpb.py`, `fibers.py`, and `lar.py`:
  - WLS: `WLSABSLENGTH`, `WLSCOMPONENT`, `WLSTIMECONSTANT`, optional
    `WLSMEANNUMBERPHOTONS`
  - Scintillation: define a {class}`legendoptics.scintillate.ScintConfig` and
    use {func}`legendoptics.pyg4utils.pyg4_def_scint_by_particle_type`.

4. Include data files

- Place spectral data in an importable package (e.g., the `data` directory in
  the pakage `mypkg` is importable as `mypkg.data`)
- Format: first line header with units (`# unit1 unit2`), then pairs of numbers;
  comments allowed after `#` (after the header line).

## CLI helper

A small CLI (defined in {mod}`legendoptics.cli`) can generate
G4GeneralParticleSource emission spectra:

```
legend-pygeom-optics g4gps lar_emission out.mac
legend-pygeom-optics g4gps pen_emission out.mac
legend-pygeom-optics g4gps fiber_emission out.mac
```

This uses the same emission spectra as the attachers.

## Practical examples

- pygeomtools (central material definitions):
  [legend-exp/legend-pygeom-tools](https://github.com/legend-exp/legend-pygeom-tools)
- LEGEND-200:
  [legend-exp/legend-pygeom-l200](https://github.com/legend-exp/legend-pygeom-l200)
- LEGEND-1000:
  [legend-exp/legend-pygeom-l1000](https://github.com/legend-exp/legend-pygeom-l1000)

The three repositories illustrate how materials are defined with their
attachers, and how to manage optical surfaces and properties consistently across
a large detector model.

## Tips and best practices

- Always use pint quantities and the provided attach helpers to avoid unit
  mistakes.
- When interpolating spectra, use `InterpolatingGraph` to handle extrapolation
  bounds predictably and similar to Geant4.
- For wavelength/energy conversions, use a pint context, e.g.
  `with u.context("sp")`.
- Prefer overriding pluggable functions instead of patching code. This keeps
  your runs reproducible and centralized.
- Keep emission spectra zeroed at sampling boundaries to avoid artifacts (see
  `pen.py`, `fibers.py`, `lar.py` patterns).
