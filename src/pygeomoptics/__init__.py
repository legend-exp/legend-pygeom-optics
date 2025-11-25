"""A collection of optical properties of materials used in the LEGEND experiment."""

from __future__ import annotations

from pygeomoptics import (
    copper,
    fibers,
    germanium,
    lar,
    nylon,
    pen,
    pmts,
    silicon,
    tetratex,
    tpb,
    tyvek,
    vm2000,
    water,
)
from pygeomoptics._version import version as __version__

# do not import pyg4utils here, to not enforce the dependency on pyg4ometry.

__all__ = [
    "__version__",
    "copper",
    "fibers",
    "germanium",
    "lar",
    "nylon",
    "pen",
    "pmts",
    "pyg4utils",  # lazy import!
    "silicon",
    "tetratex",
    "tpb",
    "tyvek",
    "vm2000",
    "water",
]


# inspired by PEP 562.
def __getattr__(name: str):
    if name in __all__:
        import importlib

        return importlib.import_module(f".{name}", __name__)
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
