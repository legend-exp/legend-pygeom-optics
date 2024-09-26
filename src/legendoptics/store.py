"""
A store that allows to manipulate and swap individual material properties with custom
implementations.

:meth:`register_pluggable` is a decorator to use above all functions defining material
properties.

Users can then replace the original implementation by their own implementations using
``replace_implementation(new_impl: Callable)``, and also switch back to the original
implementation using ``reset_implementation()`` on the decorated function object. The
original implementation is always available as ``original_impl()``.

Apart from the decorator, this store provides functions to get and reset the status of
all registered pluggable functions.
"""

from __future__ import annotations

from functools import wraps
from types import MethodType
from typing import Callable

_optical_property_store = []


def register_pluggable(fn: Callable) -> Callable:
    """Decorator that registers this function as a pluggable property function."""

    # create the new wrapper object.
    @wraps(fn)
    def wrap(*args, **kwargs):
        return wrap._impl(*args, **kwargs)

    wrap._impl = fn
    wrap._orig_impl = fn

    # add "instance methods" of the new wrapper.
    def reset_implementation(self) -> None:
        """Reset to the original function implementation."""
        self._impl = self._orig_impl

    def replace_implementation(self, new_impl: Callable) -> None:
        """Replace the underlying function implementation."""
        self._impl = new_impl

    def is_original(self) -> bool:
        """Is the underlying function the original implementation"""
        return self._impl == self._orig_impl

    def original_impl(self) -> Callable:
        """The original function implementation."""
        return self._orig_impl

    wrap.reset_implementation = MethodType(reset_implementation, wrap)
    wrap.replace_implementation = MethodType(replace_implementation, wrap)
    wrap.is_original = MethodType(is_original, wrap)
    wrap.original_impl = MethodType(original_impl, wrap)

    # store the wrapper to use it in the functions of this module.
    _optical_property_store.append(wrap)

    return wrap


def reset_all_to_original() -> None:
    """Reset all pluggable material property functions to their original implementations."""
    for entry in _optical_property_store:
        entry.reset_implementation()


def is_all_original() -> bool:
    """Get whether all pluggable material property functions use their original implementation."""
    return all(p.is_original() for p in _optical_property_store)


def get_replaced() -> list[str]:
    """Get the names of all replaced pluggable material property functions."""
    return [p.__name__ for p in _optical_property_store if not p.is_original()]
