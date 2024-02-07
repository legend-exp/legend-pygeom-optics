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


def do_plot(
    obj: Callable, plots_dir: Path, safe_name: str, options: dict[str, Any]
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
    """
    # init plot
    fig = plt.figure(figsize=(4, 2))
    ax = plt.gca()
    plt.grid()

    if "call_x" in options:
        # special case for LAr properties
        lim = [112 * u.nm, 650 * u.nm]
        if "xlim" in options:
            lim = [xl * u.nm for xl in options["xlim"]]
        x = np.linspace(*lim, num=200)
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
    for i, val in enumerate(ys):
        if isinstance(val, pint.Quantity):
            ax.set_ylabel(val.u)
            y = val.magnitude

        plotoptions = {}
        if "labels" in options:
            plotoptions["label"] = options["labels"][i]

        plt.plot(x, y, marker=".", markersize=2, linewidth=0.5, **plotoptions)

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
    fig.savefig(plots_dir / (safe_name + ".png"), dpi=300)

    return [
        ":returns:",
        f"    .. image:: plots/{safe_name}.png",
        "        :width: 400px",
    ]


def do_const(obj: Callable) -> list[str]:
    """Output the constant return value of the function.

    By default, it will call ``obj()`` and display the numerical value of the return value.
    """
    const = obj()
    description = None

    if isinstance(const, pint.Quantity):
        description = f":math:`{const.m}\\ {const.u:L~}`"
    elif isinstance(const, float):
        description = f":math:`{const}`"
    else:
        msg = ""
        raise ValueError(msg)

    if description is None:
        return [""]

    return [f":returns: constant value {description}"]


def process_docstring(
    app: Sphinx,
    what: str,
    name: str,
    obj: Any,
    options: Any,
    lines: list[str],
) -> None:
    if inspect.isclass(obj) or what != "function":
        return
    if not callable(obj):
        return

    plots_dir = Path(app.srcdir) / "api" / "plots"

    # this only appears to be a sphix directive, but the parsing here is very simple.
    plot_token = ".. optics-plot::"
    const_token = ".. optics-const::"

    i = 0
    orig_lines = lines.copy()
    for line in orig_lines:
        if line.startswith(plot_token):
            plots_dir.mkdir(exist_ok=True, parents=True)
            safe_name = "".join(c for c in name if c.isalnum() or c in (".", "-", "_"))

            opt_string = line[len(plot_token) :]
            opts = []
            if opt_string != "":
                opts = ast.literal_eval(opt_string)

            # replace the custom 'directive' with an actually supported reST directive.
            replace = do_plot(obj, plots_dir, safe_name, opts)
            lines[i : i + 1] = replace
            i += len(replace)
        elif line.startswith(const_token):
            # replace the custom 'directive' with actually supported reST content.
            replace = do_const(obj)
            lines[i : i + 1] = replace
            i += len(replace)
        else:
            i += 1


def setup(app: Sphinx) -> dict[str, bool]:
    """Register this sphinx extension."""
    app.connect("autodoc-process-docstring", process_docstring)
    return {"parallel_read_safe": True}
