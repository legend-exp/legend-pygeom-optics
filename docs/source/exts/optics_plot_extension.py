from __future__ import annotations

import ast
import inspect
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pint
from sphinx.application import Sphinx

from pygeomoptics.plot import plot_callable

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
        Change the behaviour of the plot. see :meth:`plot_callable`

        standalone
            if True, do not output the return value preamble (i.e. to embed more than one plot)
    """
    file = plots_dir / (safe_name + ".png")

    plot_callable(obj, file, options)

    lines = []
    if not options.get("standalone", False):
        lines += [":returns:"]
    return [
        *lines,
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
    elif hasattr(const, "__to_optics_const__"):
        description = const.__to_optics_const__()
    else:
        msg = f"unsupported type {type(const)}"
        raise RuntimeError(msg)

    if description is None:
        return [""]

    if not isinstance(description, str):
        return [f":returns: constant value {description[0]}", *description[1:]]
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

    # this only appears to be a sphinx directive, but the parsing here is very simple.
    plot_token = ".. optics-plot::"
    const_token = ".. optics-const::"

    i = 0
    plot_idx = 0
    orig_lines = lines.copy()
    for line in orig_lines:
        if line.startswith(plot_token):
            plots_dir.mkdir(exist_ok=True, parents=True)
            safe_name = "".join(c for c in name if c.isalnum() or c in (".", "-", "_"))
            safe_name += f"_{plot_idx}"

            opt_string = line[len(plot_token) :]
            opts = {}
            if opt_string != "":
                opts = ast.literal_eval(opt_string)
            assert isinstance(opts, dict)

            # replace the custom 'directive' with an actually supported reST directive.
            replace = do_plot(obj, plots_dir, safe_name, opts)
            lines[i : i + 1] = replace
            i += len(replace)
            plot_idx += 1
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
