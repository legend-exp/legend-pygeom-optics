# legend-pygeom-optics

[![PyPI](https://img.shields.io/pypi/v/legend-pygeom-optics?logo=pypi)](https://pypi.org/project/legend-pygeom-optics/)
![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/legend-exp/legend-pygeom-optics?logo=git)
[![GitHub Workflow Status](https://img.shields.io/github/checks-status/legend-exp/legend-pygeom-optics/main?label=main%20branch&logo=github)](https://github.com/legend-exp/legend-pygeom-optics/actions)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Codecov](https://img.shields.io/codecov/c/github/legend-exp/legend-pygeom-optics?logo=codecov)](https://app.codecov.io/gh/legend-exp/legend-pygeom-optics)
![GitHub issues](https://img.shields.io/github/issues/legend-exp/legend-pygeom-optics?logo=github)
![GitHub pull requests](https://img.shields.io/github/issues-pr/legend-exp/legend-pygeom-optics?logo=github)
![License](https://img.shields.io/github/license/legend-exp/legend-pygeom-optics)
[![Read the Docs](https://img.shields.io/readthedocs/legend-pygeom-optics?logo=readthedocs)](https://legend-pygeom-optics.readthedocs.io)

This package contains a collection of optical properties of materials used in
the LEGEND experiment.

As a common interface, each optical property gets its own defining function in
the material's module. Those functions can be used directly to just retrieve the
value(s) of the property. Most property definitions contain unit information via
the [`pint`](https://pint.readthedocs.io) package. For a full list of defined
properties refer to the
[package documentation](https://legend-pygeom-optics.readthedocs.io).

## Export of Geant4 material properties

To ease the use in Geant4-based simulations, every module also provides
functions to be used with [`pyg4ometry`](https://pyg4ometry.readthedocs.io).
These functions' names are always prefixed with `pyg4_`. All other functions of
this module can be used without having `pyg4ometry` installed.

> [!NOTE]
>
> The provided functions for interoperability with Geant4 generally **only
> support Geant4 version 11.0 and later.**.
