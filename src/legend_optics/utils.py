import pint
from importlib_resources import files
from numpy.typing import NDArray

u = pint.get_application_registry()


def readdatafile(filename: str) -> tuple[NDArray, NDArray]:
    """Read ``(x, y)`` data points from `filename` with units.

    Accepted file format ::

        # unit1 unit2
        0.23453 2.3456
        0.49678 3.6841
        ...

    Units in the header must be parseable as :mod:`pint` units.
    """
    x = []
    y = []
    lines = files("legend_optics.data").joinpath(filename).read_text().split("\n")

    # parse header
    header = lines[0].lstrip()
    if header[0] != "#":
        raise RuntimeError(
            "input data file does not seem to contain header with (pint) units"
        )

    units = header.lstrip("#").split()

    lineno = 0
    for line in lines[1:-1]:
        lineno += 1
        if not line:
            continue

        val = line.split()
        if len(val) < 2:
            raise RuntimeError(f"could not parse line {lineno}: '{line}'")

        x.append(float(val[0]))
        y.append(float(val[1]))

    return (x * u[units[0]], y * u[units[1]])
