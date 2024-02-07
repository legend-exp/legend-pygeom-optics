from __future__ import annotations

import pint
import pytest

from legendoptics import lar

u = pint.get_application_registry()


def test_dielectric_constant():
    assert lar.lar_dielectric_constant_bideau_mehu(128 * u.nm) == pytest.approx(
        1.993, rel=1e-3
    )
    assert lar.lar_dielectric_constant_cern2020(128 * u.nm) == pytest.approx(
        1.846, rel=1e-3
    )
