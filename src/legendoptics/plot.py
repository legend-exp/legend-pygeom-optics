from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pint
from pint import Quantity

u = pint.get_application_registry()
u.setup_matplotlib(True)


def plot_continuous_prop(ax: plt.Axes, prop: Callable, x: Quantity, param_dict=None):
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
    if param_dict is None:
        param_dict = {}
    ax.grid(True)
    out = ax.plot(x, prop(x), **param_dict)
    plt.tight_layout()
    return out


def plot_discrete_prop(ax: plt.Axes, prop: Callable, param_dict=None):
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
    if param_dict is None:
        param_dict = {}
    ax.grid(True)
    out = ax.plot(*prop(), marker="o", **param_dict)
    plt.tight_layout()
    return out


def plot_callable(
    obj: Callable,
    plot_file: Path,
    options: dict[str, Any],
    *,
    ax=None,
    plotoptions: dict | None = None,
) -> list[str]:
    """Create a plot from the given optical property function.

    By default, it will call ``obj()`` and unpack the result into an x-vector, and multiple
    y-vectors. All y-vectors will pe plotted together into one output file.

    Other Parameters
    ----------------
    options
        Change the behaviour of the plot.

        xlim
            Set the plot's x axis limits, see :py:func:`matplotlib.pyplot.xlim`
        ylim
            Set the plot's y axis limits, see :py:func:`matplotlib.pyplot.ylim`
        xscale
            Set the plot's x axis scaling, see :py:func:`matplotlib.pyplot.xscale`
        yscale
            Set the plot's y axis scaling, see :py:func:`matplotlib.pyplot.yscale`
        labels
            Tuple of labels that will be applied if more than one y vector is returned.
        call_x
            Differing from the default behavior above, the function will be called with an
            x vector of wavelengths in the optical range (as :py:class:`pint.Quantity`).
            All return values are interpreted as y vectors.
        ret_offset
            use the argument numbered by this (default: 0) as the first argument that will
            be treated as an x or y vector.
    """
    # init plot
    has_ax = ax is not None
    if not has_ax:
        fig = plt.figure(figsize=(4, 2))
        ax = plt.gca()
    ax.grid(True)

    plotoptions = plotoptions if plotoptions is not None else {}

    ret_offset = options.get("ret_offset", 0)
    if "call_x" in options:
        # special case for LAr properties
        lim = [112 * u.nm, 650 * u.nm]
        if "xlim" in options:
            lim = [xl * u.nm for xl in options["xlim"]]
        x = np.linspace(*lim, num=200)
        ys = obj(x)
        ys = ys[ret_offset:]
        # wrap the result in a tuple, if needed
        ys = ys if isinstance(ys, tuple) else (ys,)
    else:
        data = obj()
        x = data[ret_offset]
        # ys holds one or more
        ys = data[ret_offset + 1 :]

    if isinstance(x, pint.Quantity):
        if "xunit" in options:
            x = x.to(options["xunit"])
        ax.set_xlabel(f"{x.u:~}")
        x = x.magnitude

    # plot all supplied data vectors
    for i, val in enumerate(ys):
        y = val
        if isinstance(y, pint.Quantity):
            if "yunit" in options:
                y = y.to(options["yunit"])
            ax.set_ylabel(f"{y.u:~}")
            y = y.magnitude
        elif not isinstance(y, np.ndarray | list):
            msg = f"unsupported y-vector type {type(y)} for plot {plot_file}"
            raise ValueError(msg)

        po = {"marker": ".", "markersize": 2, "linewidth": 0.5}
        if "labels" in options:
            po["label"] = options["labels"][i]

        ax.plot(x, y, **(po | plotoptions))

    if len(ys) > 1 and "labels" in options:
        ax.legend()

    # adjust plotting options
    if "ylim" in options:
        ax.set_ylim(options["ylim"])
    if "xlim" in options:
        ax.set_xlim(options["xlim"])
    if "yscale" in options:
        ax.set_yscale(options["yscale"])
    if "xscale" in options:
        ax.set_xscale(options["xscale"])
    if "xlabel" in options:
        ax.set_xlabel(
            options["xlabel"]
            + (f" [{ax.get_xlabel()}]" if ax.get_xlabel() != "" else "")
        )
    if "ylabel" in options:
        ax.set_ylabel(
            options["ylabel"]
            + (f" [{ax.get_ylabel()}]" if ax.get_ylabel() != "" else "")
        )
    if "yright" in options:
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()

    # export figure to the filesystem
    if not has_ax and plot_file is not None:
        fig.tight_layout(pad=0.3)
        fig.savefig(plot_file, dpi=300)

        plt.close()
