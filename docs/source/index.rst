Welcome to legend-geom-optics's documentation!
==============================================

This package contains a collection of optical properties of materials used in the `LEGEND`_ experiment.

As a common interface, each optical property gets its own defining function in the material's module. Those functions can be used directly to just retrieve the value(s) of the property. Most property definitions contain unit information via the :py:mod:`pint` package. For a full list of defined properties see the API reference.

To ease the use in `Geant4`_-based simulations, every module defines one or more functions prefixed with :code:`pyg4_attach_`. Those functions are to be used with `pyg4ometry`_ and will attach the listed properties to a material or surface instance.

.. _LEGEND: https://legend-exp.org/
.. _Geant4: https://geant4.web.cern.ch/
.. _pyg4ometry: https://pyg4ometry.readthedocs.io/en/stable/index.html

Table of Contents
-----------------

.. toctree::
   :maxdepth: 1

   Package API reference <api/modules>

Usage example
-------------

This example demonstrates how to use the :py:mod:`legendoptics.lar` submodule to add important optical properties to a pyg4ometry module instance:

.. code-block:: python

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


Each property can be used separately:

.. code-block:: python

    import legendoptics.lar

    lifetimes = legendoptics.lar.lar_lifetimes(triplet_lifetime_method="legend200-llama")
    print(lifetimes.triplet)
