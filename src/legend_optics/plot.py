from typing import Callable

from numpy.typing import NDArray
import matplotlib.pyplot as plt
import pint
from pint import Quantity

pint.get_application_registry().setup_matplotlib(True)


def plot_continuous_prop(ax: plt.Axes, prop: Callable, x: Quantity[NDArray], param_dict={}):
    """Plot continuous property.

    Plots the `prop` function on values `x` with matplotlib's :func:`plot`.

    Parameters
    ----------
    ax
        axes instance.
    prop
        function defining continuous property, such as for example
        :func:`.lar.lar_refractive_index`.
    x
        array of domain values for which `prop` is shown.
    param_dict
        dictionary defining custom matplotlib settings to be passed to
        :func:`plot`.
    """
    ax.grid(True)
    out = ax.plot(x, prop(x), **param_dict)
    plt.tight_layout()
    return out


def plot_discrete_prop(ax: plt.Axes, prop: Callable, param_dict={}):
    """Plot discrete property.

    Unpacks what returned by the `prop` function and feeds it to matplotlib's
    :func:`plot`.

    Parameters
    ----------
    ax
        axes instance.
    prop
        function defining discrete properties, such as for example
        :func:`.lar.lar_emission_spectrum`.
    param_dict
        dictionary defining custom matplotlib settings to be passed to
        :func:`plot`.
    """
    ax.grid(True)
    out = ax.plot(*prop(), marker="o", **param_dict)
    plt.tight_layout()
    return out
