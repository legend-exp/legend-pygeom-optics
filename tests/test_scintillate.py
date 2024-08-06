"""Test the attaching of all properties to Geant4 materials."""

from __future__ import annotations

import numpy as np
import pint

from legendoptics import lar, pen
from legendoptics import scintillate as sc

u = pint.get_application_registry()


def test_scintillate_lar():
    rng = np.random.default_rng()

    params = sc.precompute_scintillation_params(
        lar.lar_scintillation_params(),
        lar.lar_lifetimes().as_tuple(),
    )
    part_e = sc.particle_to_index("electron")
    part_ion = sc.particle_to_index("ion")

    sc.scintillate_local(params, part_e, 10, rng)
    sc.scintillate_local(params, part_ion, 10, rng)

    x0 = np.array([0, 0, 0], dtype=np.float64)
    x1 = np.array([0, 0, 1], dtype=np.float64)

    sc.scintillate(params, x0, x1, 0.1, 0.09, 1234.5, part_e, -1, 1000, rng)
    sc.scintillate(params, x0, x1, 0.1, 0.09, 1234.5, part_ion, -1, 1000, rng)


def test_scintillate_pen():
    rng = np.random.default_rng()

    params = sc.precompute_scintillation_params(
        pen.pen_scintillation_params(),
        (pen.pen_scint_timeconstant(),),
    )
    part_e = sc.particle_to_index("electron")
    part_ion = sc.particle_to_index("ion")

    sc.scintillate_local(params, part_e, 10, rng)
    sc.scintillate_local(params, part_ion, 10, rng)
