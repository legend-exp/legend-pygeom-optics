"""Test the exporting of G4GPS spectra."""

from __future__ import annotations

import pint

u = pint.get_application_registry()


def test_g4gps_lar(tmp_path) -> None:
    import legendoptics.lar

    legendoptics.lar.g4gps_lar_emissions_spectrum(tmp_path / "lar.csv")


def test_g4gps_pen(tmp_path) -> None:
    import legendoptics.pen

    legendoptics.pen.g4gps_pen_emissions_spectrum(tmp_path / "pen.csv")
