from __future__ import annotations

import ast
import inspect
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
import numpy as np
import pint
from sphinx.application import Sphinx

u = pint.get_application_registry()


def do_plot(obj: Callable, fn: str, options: dict[str, Any]):
    """Create a plot from the given optical property function.

    By default, it will call ``obj()`` and unpack the result into an x-vector, and multiple
    y-vectors. All y-vectors will pe plotted together into one output file.

    Other parameters
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
    """
    # init plot
    fig = plt.figure(figsize=(4, 2))
    ax = plt.gca()
    plt.grid()

    if "call_x" in options:
        # special case for LAr properties
        x = np.linspace(112 * u.nm, 650 * u.nm, num=200)
        ys = obj(x)
        # wrap the result in a tuple, if needed
        ys = ys if isinstance(ys, tuple) else (ys,)
    else:
        data = obj()
        x = data[0]
        # ys holds one or more
        ys = data[1:]

    if isinstance(x, pint.Quantity):
        ax.set_xlabel(x.u)
        x = x.magnitude

    # plot all supplied data vectors
    i = 0
    for y in ys:
        if isinstance(y, pint.Quantity):
            ax.set_ylabel(y.u)
            y = y.magnitude

        plotoptions = {}
        if "labels" in options:
            plotoptions["label"] = options["labels"][i]

        plt.plot(x, y, marker=".", markersize=2, linewidth=0.5, **plotoptions)
        i += 1

    if len(ys) > 1 and "labels" in options:
        plt.legend()

    # adjust plotting options
    if "ylim" in options:
        ax.set_ylim(options["ylim"])
    if "xlim" in options:
        ax.set_xlim(options["xlim"])
    if "yscale" in options:
        ax.set_yscale(options["yscale"])
    if "xscale" in options:
        ax.set_xscale(options["xscale"])

    # export figure to the filesystem
    fig.tight_layout(pad=0.3)
    fig.savefig(fn)


def process_docstring(
    app: Sphinx, what: str, name: str, obj: Any, options: Any, lines: list[str]
) -> None:
    if inspect.isclass(obj) or what != "function":
        return
    if not callable(obj):
        return

    plots_dir = Path(app.srcdir) / "api" / "plots"

    # this only appears to be a sphix directive, but the parsing here is very simple.
    plot_token = ".. optics-plot::"

    i = 0
    for line in lines:
        if line.startswith(plot_token):
            plots_dir.mkdir(exist_ok=True, parents=True)
            safe_name = "".join(c for c in name if c.isalnum() or c in (".", "-", "_"))
            fn = plots_dir / (safe_name + ".png")

            opt_string = line[len(plot_token) :]
            opts = []
            if opt_string != "":
                opts = ast.literal_eval(opt_string)

            do_plot(obj, fn, opts)

            # replace the custom 'directive' with an actually supported reST directive.
            lines[i] = f".. image:: plots/{safe_name}.png"
        i += 1


def setup(app: Sphinx) -> dict[str, bool]:
    """Register this sphinx extension."""
    app.connect("autodoc-process-docstring", process_docstring)
    return dict(parallel_read_safe=True)
