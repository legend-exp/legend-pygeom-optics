"""Test the exporting of G4GPS spectra."""

from __future__ import annotations

import numpy as np
import pint

u = pint.get_application_registry()


def test_g4gps(tmp_path) -> None:
    from legendoptics.utils import g4gps_write_emission_spectrum

    λ_peak = np.array([200, 300, 400]) * u.nm
    scint_em = np.array([10, 11, 12]) * u.dimensionless
    tmp_file = tmp_path / "test.csv"
    g4gps_write_emission_spectrum(tmp_file, False, λ_peak, scint_em, "test")

    saved = np.loadtxt(tmp_file)
    # ensure monotonicity (increasing _energy_ order).
    assert saved[0, 0] < saved[1, 0]
    assert saved[1, 0] < saved[2, 0]
    assert saved[0, 1] > saved[1, 1]
    assert saved[1, 1] > saved[2, 1]


def test_g4gps_lar(tmp_path) -> None:
    import legendoptics.lar

    legendoptics.lar.g4gps_lar_emissions_spectrum(tmp_path / "lar.csv", False)
    legendoptics.lar.g4gps_lar_emissions_spectrum(tmp_path / "lar.mac", True)


def test_g4gps_pen(tmp_path) -> None:
    import legendoptics.pen

    legendoptics.pen.g4gps_pen_emissions_spectrum(tmp_path / "pen.csv", False)
    legendoptics.pen.g4gps_pen_emissions_spectrum(tmp_path / "pen.mac", True)
