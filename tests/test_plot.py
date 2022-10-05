import matplotlib.pyplot as plt
import numpy as np
import pint

from legend_optics.lar import lar_emission_spectrum, lar_refractive_index
from legend_optics.plot import plot_continuous_prop, plot_discrete_prop

u = pint.get_application_registry()


def test_plot_continuous_prop():
    fig, ax = plt.subplots()
    plot_continuous_prop(ax, lar_refractive_index, np.arange(110, 300, 1) * u.nm)


def test_plot_discrete_prop():
    fig, ax = plt.subplots()
    plot_discrete_prop(ax, lar_emission_spectrum)
