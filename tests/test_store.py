from __future__ import annotations

from legendoptics import store
from legendoptics.fibers import (
    fiber_cladding1_refractive_index,
    fiber_core_refractive_index,
)


def test_store():
    assert fiber_core_refractive_index() == 1.6
    assert fiber_core_refractive_index.is_original()

    # test replacing the implementation.
    orig = fiber_core_refractive_index.original_impl()
    fiber_core_refractive_index.replace_implementation(lambda: 1234)
    assert fiber_core_refractive_index() == 1234
    assert not fiber_core_refractive_index.is_original()
    assert fiber_core_refractive_index.original_impl() is orig

    fiber_cladding1_refractive_index.replace_implementation(lambda: 1)

    # test single reset.
    fiber_core_refractive_index.reset_implementation()
    assert fiber_core_refractive_index() == 1.6
    assert fiber_core_refractive_index.is_original()
    # ... other properties are not reset:
    assert not store.is_all_original()
    assert fiber_cladding1_refractive_index() == 1

    # test global reset.
    fiber_core_refractive_index.replace_implementation(lambda: 5678)
    assert not fiber_core_refractive_index.is_original()

    assert store.get_replaced() == [
        "fiber_cladding1_refractive_index",
        "fiber_core_refractive_index",
    ]

    store.reset_all_to_original()
    assert fiber_core_refractive_index() == 1.6
    assert fiber_core_refractive_index.is_original()
    assert fiber_cladding1_refractive_index() == 1.49  # now also reset.
    assert store.get_replaced() == []
    assert store.is_all_original()
