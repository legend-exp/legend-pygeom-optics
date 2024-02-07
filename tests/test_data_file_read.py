from __future__ import annotations

from legendoptics.utils import readdatafile


def test_read_data_file():
    readdatafile("lar_emission_heindl2010.dat")
