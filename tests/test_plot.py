from __future__ import annotations

import ast

import matplotlib.pyplot as plt
import numpy as np
import pint

from pygeomoptics.copper import copper_reflectivity
from pygeomoptics.lar import lar_refractive_index
from pygeomoptics.plot import plot_callable, plot_continuous_prop, plot_discrete_prop

u = pint.get_application_registry()


def test_plot_continuous_prop():
    _fig, ax = plt.subplots()
    plot_continuous_prop(ax, lar_refractive_index, np.arange(110, 300, 1) * u.nm)


def test_plot_discrete_prop():
    _fig, ax = plt.subplots()
    plot_discrete_prop(ax, copper_reflectivity)


def test_plot_for_docs(tmp_path):
    from pygeomoptics.store import _optical_property_store

    # this should mirror the code in the sphinx extension.
    plot_token = ".. optics-plot::"
    for p in _optical_property_store:
        lines = p.__doc__.split("\n")
        lines = [line.strip() for line in lines if line.strip().startswith(plot_token)]
        for line in lines:
            opt_string = line[len(plot_token) :].strip()
            opts = {}
            if opt_string != "":
                opts = ast.literal_eval(opt_string)
            assert isinstance(opts, dict)

            plot_callable(p, tmp_path / "test.png", opts)
