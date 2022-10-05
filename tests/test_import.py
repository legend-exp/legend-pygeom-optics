def test_import():
    import legend_optics  # noqa: F401
    import legend_optics.fibers  # noqa: F401
    import legend_optics.lar  # noqa: F401
    import legend_optics.tetratex  # noqa: F401
    import legend_optics.tpb  # noqa: F401

    legend_optics.lar.lar_fano_factor()

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
