"""Test the attaching of all properties to Geant4 materials."""

from __future__ import annotations

import inspect
from contextlib import contextmanager

import numpy as np
import pint

from legendoptics import lar, pen
from legendoptics import scintillate as sc

u = pint.get_application_registry()


@contextmanager
def numba_unwrap(mod, fn: str, extra_fns: tuple[str, ...] = ()):
    """Unwrap a numba function `<mod>.<fn>`, execute both and compare results.

    Parameters
    ----------
    extra_fns
        while executing the unwrapped function, also unwrap other member variables of `mod`.
    """
    numba_f = getattr(mod, fn)
    orig_f = inspect.unwrap(numba_f)

    def _compare(*args, **kwargs):
        # note: inject the same rng state (i.e. equal seeds) into both functions.
        numba_ret = numba_f(*args, rng=np.random.default_rng(0), **kwargs)

        # only replace the extra functions while the unwrapped function executes.
        extra_orig = [getattr(mod, extra_fn) for extra_fn in extra_fns]
        for extra_numba_fn, extra_fn in zip(extra_orig, extra_fns):
            setattr(mod, extra_fn, inspect.unwrap(extra_numba_fn))

        orig_ret = orig_f(*args, rng=np.random.default_rng(0), **kwargs)

        # rollback, so that the next numba invocation works.
        for extra_numba_fn, extra_fn in zip(extra_orig, extra_fns):
            setattr(mod, extra_fn, extra_numba_fn)

        assert np.allclose(numba_ret, orig_ret)

        return numba_ret

    yield _compare


def test_scintillate_lar():
    params = sc.precompute_scintillation_params(
        lar.lar_scintillation_params(),
        lar.lar_lifetimes().as_tuple(),
    )
    part_e = sc.particle_to_index("electron")
    part_ion = sc.particle_to_index("ion")

    with numba_unwrap(sc, "scintillate_local") as scintillate_local:
        scintillate_local(params, part_e, 10)
        scintillate_local(params, part_ion, 10)

    x0 = np.array([0, 0, 0], dtype=np.float64)
    x1 = np.array([0, 0, 1], dtype=np.float64)

    with numba_unwrap(sc, "scintillate", ("scintillate_local",)) as scintillate:
        scintillate(params, x0, x1, 0.1, 0.09, 1234.5, part_e, -1, 1000)
        scintillate(params, x0, x1, 0.1, 0.09, 1234.5, part_ion, -1, 1000)

def test_scintillate_point_like():
    params = sc.precompute_scintillation_params(
        lar.lar_scintillation_params(),
        lar.lar_lifetimes().as_tuple(),
    )
    part_e = sc.particle_to_index("electron")
    part_ion = sc.particle_to_index("ion")

    with numba_unwrap(sc, "scintillate_local") as scintillate_local:
        scintillate_local(params, part_e, 10)
        scintillate_local(params, part_ion, 10)

    x0 = np.array([0, 0, 0], dtype=np.float64)
    #x1 = np.array([0, 0, 1], dtype=np.float64)

    with numba_unwrap(sc, "scintillate", ("scintillate_local",)) as scintillate:
        scintillate(params, x0, None, None, None, 1234.5, part_e, -1, 1000)
        scintillate(params, x0, None, None, None, 1234.5, part_ion, -1, 1000)

def test_scintillate_pen():
    params = sc.precompute_scintillation_params(
        pen.pen_scintillation_params(),
        (pen.pen_scint_timeconstant(),),
    )
    part_e = sc.particle_to_index("electron")
    part_ion = sc.particle_to_index("ion")

    with numba_unwrap(sc, "scintillate_local") as scintillate_local:
        scintillate_local(params, part_e, 10)
        scintillate_local(params, part_ion, 10)
