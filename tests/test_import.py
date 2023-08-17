import numpy as np
import pint

u = pint.get_application_registry()


def test_import():
    import legendoptics  # noqa: F401
    import legendoptics.copper
    import legendoptics.fibers
    import legendoptics.germanium
    import legendoptics.lar
    import legendoptics.nylon
    import legendoptics.pen
    import legendoptics.silicon
    import legendoptics.tetratex
    import legendoptics.tpb

    legendoptics.lar.lar_fano_factor()
    legendoptics.lar.lar_emission_spectrum()
    wvl = np.arange(110, 400, 20) * u.nm
    legendoptics.lar.lar_dielectric_constant_bideau_mehu(wvl)
    legendoptics.lar.lar_dielectric_constant_cern2020(wvl)
    legendoptics.lar.lar_dielectric_constant(wvl)
    legendoptics.lar.lar_refractive_index(wvl)
    legendoptics.lar.lar_rayleigh(wvl, 87 * u.K)
    legendoptics.lar.lar_abs_length(wvl)
    legendoptics.lar.lar_peak_attenuation_length()
    legendoptics.lar.lar_lifetimes()
    legendoptics.lar.lar_scintillation_params()

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
    legendoptics.fibers.fiber_absorption_path_length()

    legendoptics.tetratex.tetratex_reflectivity()

    legendoptics.germanium.germanium_reflectivity()

    legendoptics.copper.copper_reflectivity()

    legendoptics.silicon.silicon_complex_rindex()

    legendoptics.nylon.nylon_refractive_index()
    legendoptics.nylon.nylon_absorption()

    legendoptics.pen.pen_refractive_index()
    legendoptics.pen.pen_quantum_efficiency()
    legendoptics.pen.pen_scint_timeconstant()
    legendoptics.pen.pen_scint_light_yield()
    legendoptics.pen.pen_wls_emission()
    legendoptics.pen.pen_absorption()
    legendoptics.pen.pen_wls_absorption()
    legendoptics.pen.pen_wls_absorption()
    legendoptics.pen.pen_wls_absorption()
