"""Tantalum.

.. [Werner2009] W. S. M. Werner, K. Glantschnig and C. Ambrosch-Draxl “Optical Constants and
    Inelastic Electron-Scattering Data for 17 Elemental Metals”. In: J. Phys. Chem. Ref. Data
    38, 1013-1092 (2009). https://doi.org/10.1063/1.3243762
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pint
from pint import Quantity

if TYPE_CHECKING:
    import pyg4ometry.geant4 as g4

from pygeomoptics import store
from pygeomoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def tantalum_complex_rindex() -> tuple[Quantity, Quantity, Quantity]:
    """Real and imaginary parts of tantalum refractive index. Measurements from [Werner2009]_.

    Returns a tuple(wavelength, Re, Im). Data files as digitized by [Polyanskiy2024]_,
    restricted to a for us useful wavelength range.

    .. optics-plot:: {'labels':('Re n','Im n')}
    """
    real = readdatafile("ta_rindex_real.dat")
    imag = readdatafile("ta_rindex_imag.dat")
    assert (real[0] == imag[0]).all()
    return real[0], real[1], imag[1]


def pyg4_tantalum_attach_complex_rindex(mat: g4.Material, reg: g4.Registry) -> None:
    """Attach the complex refractive index to the given tantalum material instance.

    See Also
    --------
    .tantalum_complex_rindex
    """
    λ, re, im = tantalum_complex_rindex()
    with u.context("sp"):
        mat.addVecPropertyPint("REALRINDEX", λ.to("eV"), re)
        mat.addVecPropertyPint("IMAGINARYRINDEX", λ.to("eV"), im)
