import pint
import numpy as np

u = pint.get_application_registry()


def test_import():
    import legend_optics  # noqa: F401
    import legend_optics.fibers  # noqa: F401
    import legend_optics.lar  # noqa: F401
    import legend_optics.tetratex  # noqa: F401
    import legend_optics.tpb  # noqa: F401

    legend_optics.lar.lar_fano_factor()
    legend_optics.lar.lar_emission_spectrum()
    wvl = np.arange(110, 400, 20) * u.nm
    legend_optics.lar.lar_dielectric_constant_bideau_mehu(wvl)
    legend_optics.lar.lar_dielectric_constant_cern2020(wvl)
    legend_optics.lar.lar_dielectric_constant(wvl)
    legend_optics.lar.lar_refractive_index(wvl)
    legend_optics.lar.lar_rayleigh(wvl, 87*u.K)

    legend_optics.tpb.tpb_quantum_efficiency()
    legend_optics.tpb.tpb_refractive_index()
    legend_optics.tpb.tpb_wls_timeconstant()
    legend_optics.tpb.tpb_wls_emission()
    legend_optics.tpb.tpb_wls_absorption()

    legend_optics.fibers.fiber_cladding2_refractive_index()
    legend_optics.fibers.fiber_cladding1_refractive_index()
    legend_optics.fibers.fiber_core_refractive_index()
    legend_optics.fibers.fiber_wls_absorption()
    legend_optics.fibers.fiber_wls_emission()
    legend_optics.fibers.fiber_wls_timeconstant()
    legend_optics.fibers.fiber_absorption_length()

    legend_optics.tetratex.tetratex_reflectivity()
