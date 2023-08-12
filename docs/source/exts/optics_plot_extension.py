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
    """Create a plot from the given optical property function."""
    # init plot
    fig = plt.figure(figsize=(4, 2))
    ax = plt.gca()
    plt.grid()

    if "call_x" in options:
        # special case for LAr properties.
        x = np.linspace(112 * u.nm, 650 * u.nm, num=200)
        y = obj(x)
    else:
        x, y = obj()

    if isinstance(x, pint.Quantity):
        ax.set_xlabel(x.u)
        x = x.magnitude
    if isinstance(y, pint.Quantity):
        ax.set_ylabel(y.u)
        y = y.magnitude
    plt.plot(x, y, marker=".", markersize=2, linewidth=0.5)

    # adjust plotting options
    if "ylim" in options:
        ax.set_ylim(options["ylim"])
    if "xlim" in options:
        ax.set_xlim(options["xlim"])
    if "yright" in options:
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
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

    plot_token = "..optics-plot::"

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
            lines[i] = f".. image:: plots/{safe_name}.png"
        i += 1


def setup(app: Sphinx) -> dict[str, bool]:
    """Register this sphinx extension."""
    app.connect("autodoc-process-docstring", process_docstring)
    return dict(parallel_read_safe=True)
