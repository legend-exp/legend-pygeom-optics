import pint
import pyg4ometry.geant4 as g4

u = pint.get_application_registry()


def _create_dummy_mat():
    reg = g4.Registry()
    mat = g4.MaterialSingleElement("dummy", 1, 1, 1, reg)
    return reg, mat


def test_pyg4_attach():
    """Test the attaching of all properties to Geant4 materials."""
    import legendoptics  # noqa: F401
    import legendoptics.copper
    import legendoptics.fibers
    import legendoptics.germanium
    import legendoptics.lar
    import legendoptics.nylon
    import legendoptics.silicon
    import legendoptics.tetratex
    import legendoptics.tpb

    u = pint.get_application_registry()

    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_rindex(mat, reg)
    legendoptics.lar.pyg4_lar_attach_attenuation(mat, reg, 90 * u.K)
    legendoptics.lar.pyg4_lar_attach_scintillation(mat, reg)

    reg, mat = _create_dummy_mat()
    legendoptics.tpb.pyg4_tpb_attach_rindex(mat, reg)
    legendoptics.tpb.pyg4_tpb_attach_wls(mat, reg)

    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_cladding1_attach_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_cladding1_attach_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_core_attach_rindex(mat, reg)
    legendoptics.fibers.pyg4_fiber_core_attach_wls(mat, reg)
    legendoptics.fibers.pyg4_fiber_core_attach_absorption(mat, reg)

    reg, mat = _create_dummy_mat()
    legendoptics.tetratex.pyg4_tetratex_attach_reflectivity(mat, reg)

    reg, mat = _create_dummy_mat()
    legendoptics.germanium.pyg4_germanium_attach_reflectivity(mat, reg)

    reg, mat = _create_dummy_mat()
    legendoptics.copper.pyg4_copper_attach_reflectivity(mat, reg)

    reg, mat = _create_dummy_mat()
    legendoptics.silicon.pyg4_silicon_attach_complex_rindex(mat, reg)

    reg, mat = _create_dummy_mat()
    legendoptics.nylon.pyg4_nylon_attach_rindex(mat, reg)
    legendoptics.nylon.pyg4_nylon_attach_absorption(mat, reg)
