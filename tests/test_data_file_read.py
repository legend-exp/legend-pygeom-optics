from __future__ import annotations

import numpy as np
import pint
import pytest

from legendoptics.pyg4utils import pyg4_sample_λ
from legendoptics.utils import InterpolatingGraph, readdatafile

u = pint.get_application_registry()


def test_readdatafile():
    qdim, qval = readdatafile("lar_emission_heindl2010.dat")

    assert isinstance(qdim, u.Quantity)
    assert isinstance(qval, u.Quantity)
    assert qdim.u == u.nm
    assert qval.u == u.dimensionless

    assert qdim.shape == (1601,)
    assert qdim.shape == qval.shape


def test_interpolating_init():
    """test various common errors in initialization."""

    with pytest.raises(ValueError, match="pint Quantities"):
        InterpolatingGraph(np.array([0, 1]), np.array([0, 1]))
    with pytest.raises(ValueError, match="1-dimensional"):
        InterpolatingGraph(
            np.array([[0, 1], [1, 1]]) * u.m, np.array([[0, 1], [1, 1]]) * u.m
        )
    with pytest.raises(ValueError, match="same shape"):
        InterpolatingGraph(np.array([0, 1]) * u.m, np.array([0, 1, 1]) * u.m)
    with pytest.raises(ValueError, match="only contain positive"):
        InterpolatingGraph(np.array([0, 1]) * u.m, np.array([0, -1]) * u.m)


def test_interpolating_options():
    """test that options do not change the behaviour in intervals they should not affect."""

    data = readdatafile("pen_wlscomponent.dat")

    interp_default = InterpolatingGraph(*data)
    interp_zero = InterpolatingGraph(*data, zero_outside=True)
    interp_bound_min = InterpolatingGraph(*data, min_idx=380 * u.nm)
    interp_bound_max = InterpolatingGraph(*data, max_idx=500 * u.nm)

    interp_zero_bounds = InterpolatingGraph(
        *data, zero_outside=True, min_idx=380 * u.nm, max_idx=500 * u.nm
    )

    # the shared interval is (380 nm, 500 nm) - all graphs should be equal there.
    points = pyg4_sample_λ(381 * u.nm, 499 * u.nm)
    assert np.all(interp_default(points) == interp_zero(points))
    assert np.all(interp_default(points) == interp_bound_min(points))
    assert np.all(interp_default(points) == interp_bound_max(points))

    assert np.all(interp_zero_bounds(points) > 0)

    # test left boundary.
    points = pyg4_sample_λ(360 * u.nm, 499 * u.nm)
    assert np.all(interp_default(points) >= interp_zero(points))
    assert np.all(interp_default(points) <= interp_bound_min(points))
    assert np.all(interp_default(points) == interp_bound_max(points))

    assert np.all((interp_zero_bounds(points) > 0) == (points > 380.56 * u.nm))

    # test right boundary.
    points = pyg4_sample_λ(381 * u.nm, 650 * u.nm)
    assert np.all(interp_default(points) >= interp_zero(points))
    assert np.all(interp_default(points) == interp_bound_min(points))
    assert np.all(interp_default(points) <= interp_bound_max(points))

    assert np.all((interp_zero_bounds(points) > 0) == (points < 499.99 * u.nm))


def test_interpolating_linear():
    """compare with a simple linear interpolation."""

    data = readdatafile("pen_wlscomponent.dat")
    interp_default = InterpolatingGraph(*data)

    def _interp(idx, vals, point):
        # Find the index j such that idx[j] <= point <= idx[j+1]
        j = np.searchsorted(idx, point) - 1
        f = (j >= 0) & (j < len(idx) - 1)
        assert np.all(f), "simple linear interp can only handle in-bounds points"

        # Perform linear interpolation between points j and j+1
        ret = (point - idx[j]) * (vals[j + 1] - vals[j]) / (idx[j + 1] - idx[j])
        ret += vals[j]
        assert ret.shape == point.shape
        return ret

    points = pyg4_sample_λ(381 * u.nm, 630 * u.nm, sample_count=25000)
    result_simple = _interp(interp_default.idx, interp_default.vals, points)
    result_graph = interp_default(points)

    assert np.all(np.isclose(result_simple, result_graph, rtol=1e-10, atol=1e-10))
