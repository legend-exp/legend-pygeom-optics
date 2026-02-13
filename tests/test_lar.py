from __future__ import annotations

import numpy as np
import pint
import pytest

from pygeomoptics import lar

u = pint.get_application_registry()


def test_dielectric_constant():
    assert lar.lar_dielectric_constant_bideau_mehu(128 * u.nm) == pytest.approx(
        1.993, rel=1e-3
    )
    assert lar.lar_dielectric_constant_cern2020(128 * u.nm) == pytest.approx(
        1.846, rel=1e-3
    )


def test_two_component_absorption_length():
    _, peak_absl, λ_full, _, abslength, _ = lar.lar_calculate_attenuation(
        88.8 * u.K,
        lar_dielectric_method="cern2020",
        absorption_enabled_or_length="legend200-llama-two-components",
    )

    assert peak_absl.m == pytest.approx(5.6, rel=1e-3)

    assert abslength[λ_full == np.max(λ_full)][0].m == pytest.approx(1000, rel=1e-3)
