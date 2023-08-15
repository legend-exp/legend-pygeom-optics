"""Silicon.

.. [Phillip1960] H. R. Phillip and E. A. Taft “Optical Constants of Silicon in the Region 1 to 10 eV”.
   In: Phys. Rev. 120 (1 Oct. 1960)
   https://doi.org/10.1103/PhysRev.120.37
"""

from __future__ import annotations

import logging

import pint
from pint import Quantity

from legendoptics.utils import readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def silicon_complex_rindex() -> tuple[Quantity, Quantity, Quantity]:
    """Real and imaginary parts as tuple(wavelength, Re, Im). Measurements from [Phillip1960]_.

    .. optics-plot:: {'labels':('Re n','Im n')}
    """
    real = readdatafile("si_rindex_real.dat")
    imag = readdatafile("si_rindex_imag.dat")
    assert (real[0] == imag[0]).all()
    return real[0], real[1], imag[1]


def pyg4_silicon_attach_complex_rindex(mat, reg) -> None:
    """Attach the complex refractive index to the given silicon material instance.

    See Also
    --------
    .silicon_complex_rindex
    """
    λ, re, im = silicon_complex_rindex()
    with u.context("sp"):
        mat.addVecPropertyPint("REALRINDEX", λ.to("eV"), re)
        mat.addVecPropertyPint("IMAGINARYRINDEX", λ.to("eV"), im)
