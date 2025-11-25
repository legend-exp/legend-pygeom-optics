from __future__ import annotations

import numpy as np
import pint
import pytest

u = pint.get_application_registry()


def test_import_legacy():
    with pytest.deprecated_call():
        import legendoptics.copper
        import legendoptics.fibers
        import legendoptics.germanium
        import legendoptics.lar
        import legendoptics.nylon
        import legendoptics.pen
        import legendoptics.pmts
        import legendoptics.silica
        import legendoptics.silicon
        import legendoptics.tetratex
        import legendoptics.tpb
        import legendoptics.tyvek
        import legendoptics.ultem
        import legendoptics.vm2000
        import legendoptics.water

        import legendoptics  # noqa: F401


def test_import():
    import pygeomoptics
    import pygeomoptics.copper
    import pygeomoptics.fibers
    import pygeomoptics.germanium
    import pygeomoptics.lar
    import pygeomoptics.nylon
    import pygeomoptics.pen
    import pygeomoptics.pmts
    import pygeomoptics.silica
    import pygeomoptics.silicon
    import pygeomoptics.tetratex
    import pygeomoptics.tpb
    import pygeomoptics.tyvek
    import pygeomoptics.ultem
    import pygeomoptics.vm2000
    import pygeomoptics.water

    pygeomoptics.lar.lar_fano_factor()
    λ_peak = np.arange(116, 141, 20) * u.nm
    pygeomoptics.lar.lar_emission_spectrum(λ_peak)
    λ = np.arange(110, 400, 20) * u.nm
    pygeomoptics.lar.lar_dielectric_constant_bideau_mehu(λ)
    pygeomoptics.lar.lar_dielectric_constant_cern2020(λ)
    pygeomoptics.lar.lar_dielectric_constant(λ)
    pygeomoptics.lar.lar_refractive_index(λ)
    pygeomoptics.lar.lar_rayleigh(λ, 87 * u.K)
    pygeomoptics.lar.lar_abs_length(λ)
    pygeomoptics.lar.lar_peak_attenuation_length()
    pygeomoptics.lar.lar_lifetimes()
    pygeomoptics.lar.lar_scintillation_params()

    pygeomoptics.tpb.tpb_quantum_efficiency()
    pygeomoptics.tpb.tpb_refractive_index()
    pygeomoptics.tpb.tpb_wls_timeconstant()
    pygeomoptics.tpb.tpb_wls_emission()
    pygeomoptics.tpb.tpb_wls_absorption()

    pygeomoptics.fibers.fiber_cladding2_refractive_index()
    pygeomoptics.fibers.fiber_cladding1_refractive_index()
    pygeomoptics.fibers.fiber_core_refractive_index()
    pygeomoptics.fibers.fiber_wls_absorption()
    pygeomoptics.fibers.fiber_wls_emission()
    pygeomoptics.fibers.fiber_wls_timeconstant()
    pygeomoptics.fibers.fiber_absorption_length()
    pygeomoptics.fibers.fiber_absorption_path_length()
    pygeomoptics.fibers.fiber_core_scint_light_yield()
    pygeomoptics.fibers.fiber_core_scintillation_params()

    pygeomoptics.tetratex.tetratex_reflectivity()
    pygeomoptics.tyvek.tyvek_reflectivity()

    pygeomoptics.germanium.germanium_reflectivity()

    pygeomoptics.copper.copper_reflectivity()

    pygeomoptics.silicon.silicon_complex_rindex()

    pygeomoptics.nylon.nylon_refractive_index()
    pygeomoptics.nylon.nylon_absorption()

    pygeomoptics.pen.pen_refractive_index()
    pygeomoptics.pen.pen_quantum_efficiency()
    pygeomoptics.pen.pen_scint_timeconstant()
    pygeomoptics.pen.pen_scint_light_yield()
    pygeomoptics.pen.pen_wls_emission()
    pygeomoptics.pen.pen_absorption()
    pygeomoptics.pen.pen_wls_absorption()
    pygeomoptics.pen.pen_wls_absorption()
    pygeomoptics.pen.pen_wls_absorption()
    pygeomoptics.pen.pen_scintillation_params()

    pygeomoptics.ultem.ultem_refractive_index()
    pygeomoptics.ultem.ultem_absorption()

    λ_silica = np.arange(200, 400, 20) * u.nm
    pygeomoptics.silica.silica_refractive_index(λ_silica)

    pygeomoptics.water.water_refractive_index()
    pygeomoptics.water.water_absorption()

    pygeomoptics.vm2000.vm2000_calculate_wls_mfp(0.075)
    pygeomoptics.vm2000.vm2000_refractive_index()
    pygeomoptics.vm2000.vm2000_absorption_length()
    pygeomoptics.vm2000.vm2000_parameters()
    pygeomoptics.vm2000.vm2000_scint_timeconstant()

    pygeomoptics.pmts.pmt_acryl_absorption_length()
    pygeomoptics.pmts.pmt_acryl_refractive_index()
    pygeomoptics.pmts.pmt_air_absorption_length()
    pygeomoptics.pmts.pmt_air_refractive_index()
    pygeomoptics.pmts.pmt_borosilicate_absorption_length()
    pygeomoptics.pmts.pmt_borosilicate_refractive_index()
    pygeomoptics.pmts.pmt_etl9354kb_photocathode_collection_efficiency()
    pygeomoptics.pmts.pmt_r7081_photocathode_collection_efficiency()
    pygeomoptics.pmts.pmt_etl9354kb_photocathode_efficiency()
    pygeomoptics.pmts.pmt_r7081_photocathode_efficiency()
    pygeomoptics.pmts.pmt_photocathode_reflectivity()
    pygeomoptics.pmts.pmt_steel_efficiency()
    pygeomoptics.pmts.pmt_steel_reflectivity()
