from __future__ import annotations

import numpy as np
import pint

u = pint.get_application_registry()


def test_import():
    import legendoptics
    import legendoptics.copper
    import legendoptics.fibers
    import legendoptics.germanium
    import legendoptics.lar
    import legendoptics.nylon
    import legendoptics.pen
    import legendoptics.pmts
    import legendoptics.silicon
    import legendoptics.tetratex
    import legendoptics.tpb
    import legendoptics.tyvek
    import legendoptics.vm2000
    import legendoptics.water

    legendoptics.lar.lar_fano_factor()
    λ_peak = np.arange(116, 141, 20) * u.nm
    legendoptics.lar.lar_emission_spectrum(λ_peak)
    λ = np.arange(110, 400, 20) * u.nm
    legendoptics.lar.lar_dielectric_constant_bideau_mehu(λ)
    legendoptics.lar.lar_dielectric_constant_cern2020(λ)
    legendoptics.lar.lar_dielectric_constant(λ)
    legendoptics.lar.lar_refractive_index(λ)
    legendoptics.lar.lar_rayleigh(λ, 87 * u.K)
    legendoptics.lar.lar_abs_length(λ)
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
    legendoptics.fibers.fiber_core_scint_light_yield()
    legendoptics.fibers.fiber_core_scintillation_params()

    legendoptics.tetratex.tetratex_reflectivity()
    legendoptics.tyvek.tyvek_reflectivity()

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
    legendoptics.pen.pen_scintillation_params()

    legendoptics.water.water_refractive_index()
    legendoptics.water.water_absorption()

    legendoptics.vm2000.vm2000_calculate_wls_mfp(0.075)
    legendoptics.vm2000.vm2000_refractive_index()
    legendoptics.vm2000.vm2000_absorption_length()
    legendoptics.vm2000.vm2000_parameters()
    legendoptics.vm2000.vm2000_scint_timeconstant()

    legendoptics.pmts.pmt_acryl_absorption_length()
    legendoptics.pmts.pmt_acryl_refractive_index()
    legendoptics.pmts.pmt_air_absorption_length()
    legendoptics.pmts.pmt_air_refractive_index()
    legendoptics.pmts.pmt_borosilicate_absorption_length()
    legendoptics.pmts.pmt_borosilicate_refractive_index()
    legendoptics.pmts.pmt_etl9354kb_photocathode_collection_efficiency()
    legendoptics.pmts.pmt_r7081_photocathode_collection_efficiency()
    legendoptics.pmts.pmt_etl9354kb_photocathode_efficiency()
    legendoptics.pmts.pmt_r7081_photocathode_efficiency()
    legendoptics.pmts.pmt_photocathode_reflectivity()
    legendoptics.pmts.pmt_steel_efficiency()
    legendoptics.pmts.pmt_steel_reflectivity()
