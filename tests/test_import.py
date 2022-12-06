import numpy as np
import pint

u = pint.get_application_registry()


def test_import():
    import legendoptics  # noqa: F401
    import legendoptics.fibers
    import legendoptics.lar
    import legendoptics.tetratex
    import legendoptics.tpb
    import legendoptics.germanium
    import legendoptics.copper
    import legendoptics.nylon

    legendoptics.lar.lar_fano_factor()
    legendoptics.lar.lar_emission_spectrum()
    wvl = np.arange(110, 400, 20) * u.nm
    legendoptics.lar.lar_dielectric_constant_bideau_mehu(wvl)
    legendoptics.lar.lar_dielectric_constant_cern2020(wvl)
    legendoptics.lar.lar_dielectric_constant(wvl)
    legendoptics.lar.lar_refractive_index(wvl)
    legendoptics.lar.lar_rayleigh(wvl, 87 * u.K)

    legendoptics.tpb.tpb_quantum_efficiency()
    legendoptics.tpb.tpb_refractive_index()
    legendoptics.tpb.tpb_wls_timeconstant()
    legendoptics.tpb.tpb_wls_emission()
    legendoptics.tpb.tpb_wls_absorption()

    legendoptics.fibers.fiber_cladding2_refractive_index()
    legendoptics.fibers.fiber_cladding1_refractive_index()
    legendoptics.fibers.fiber_core_refractive_index()
    legendoptics.fibers.fiber_wls_absorption()
    legendoptics.fibers.fiber_wls_emission()
    legendoptics.fibers.fiber_wls_timeconstant()
    legendoptics.fibers.fiber_absorption_length()

    legendoptics.tetratex.tetratex_reflectivity()

    legendoptics.germanium.germanium_reflectivity()

    legendoptics.copper.copper_reflectivity()

    legendoptics.nylon.nylon_refractive_index()
    legendoptics.nylon.nylon_absorption()
