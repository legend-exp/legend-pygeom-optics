from __future__ import annotations

import sys
import warnings

import pygeomoptics  # noqa: F401

sys.modules[__name__] = sys.modules["pygeomoptics"]

warnings.warn(
    "Please use `pygeomoptics` instead of `legendoptics`.", FutureWarning, stacklevel=2
)
