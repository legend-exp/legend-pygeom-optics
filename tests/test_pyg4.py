"""Test the attaching of all properties to Geant4 materials."""

from __future__ import annotations

import pint
import pyg4ometry.geant4 as g4

u = pint.get_application_registry()


def _create_dummy_mat() -> tuple[g4.Registry, g4.Material]:
    reg = g4.Registry()
    mat = g4.MaterialSingleElement("dummy", 1, 1, 1, reg)
    return reg, mat


def test_pyg4_attach_lar() -> None:
    import legendoptics.lar

    u = pint.get_application_registry()
    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_rindex(mat, reg)
    legendoptics.lar.pyg4_lar_attach_attenuation(mat, reg, 90 * u.K)
    legendoptics.lar.pyg4_lar_attach_scintillation(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_attenuation(
        mat, reg, 90 * u.K, attenuation_method_or_length=10 * u.cm
    )
    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_attenuation(
        mat,
        reg,
        90 * u.K,
        attenuation_method_or_length=10 * u.cm,
        rayleigh_enabled_or_length=90 * u.cm,
    )
    legendoptics.lar.pyg4_lar_attach_scintillation(
        mat, reg, triplet_lifetime_method=1.1
    )
    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_attenuation(
        mat, reg, 90 * u.K, absorption_enabled_or_length=False
    )
    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_attenuation(
        mat, reg, 90 * u.K, absorption_enabled_or_length=30 * u.cm
    )
    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_attenuation(
        mat,
        reg,
        90 * u.K,
        absorption_enabled_or_length=10 * u.cm,
        rayleigh_enabled_or_length=90 * u.cm,
    )
    assert "ABSORPTION" not in mat.properties
    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_attenuation(
        mat, reg, 90 * u.K, rayleigh_enabled_or_length=False
    )
    assert "RAYLEIGH" not in mat.properties


def test_pyg4_attach_tpb() -> None:
    import legendoptics.tpb

    reg, mat = _create_dummy_mat()
    legendoptics.tpb.pyg4_tpb_attach_rindex(mat, reg)
    legendoptics.tpb.pyg4_tpb_attach_wls(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.tpb.pyg4_tpb_attach_wls(mat, reg, 0.1)
    reg, mat = _create_dummy_mat()
    legendoptics.tpb.pyg4_tpb_attach_wls(mat, reg, True, "polystyrene_matrix")


def test_pyg4_attach_fibers() -> None:
    import legendoptics.fibers

    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_cladding1_attach_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_cladding2_attach_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_core_attach_rindex(mat, reg)
    legendoptics.fibers.pyg4_fiber_core_attach_wls(mat, reg)
    legendoptics.fibers.pyg4_fiber_core_attach_absorption(mat, reg)
    legendoptics.fibers.pyg4_fiber_core_attach_scintillation(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_core_attach_absorption(
        mat, reg, use_geometrical_absorption=False
    )


def test_pyg4_attach_tetratex() -> None:
    import legendoptics.tetratex

    reg, mat = _create_dummy_mat()
    legendoptics.tetratex.pyg4_tetratex_attach_reflectivity(mat, reg)


def test_pyg4_attach_tyvek() -> None:
    import legendoptics.tyvek

    reg, mat = _create_dummy_mat()
    legendoptics.tyvek.pyg4_tyvek_attach_reflectivity(mat, reg)


def test_pyg4_attach_germanium() -> None:
    import legendoptics.germanium

    reg, mat = _create_dummy_mat()
    legendoptics.germanium.pyg4_germanium_attach_reflectivity(mat, reg)


def test_pyg4_attach_copper() -> None:
    import legendoptics.copper

    reg, mat = _create_dummy_mat()
    legendoptics.copper.pyg4_copper_attach_reflectivity(mat, reg)


def test_pyg4_attach_silicon() -> None:
    import legendoptics.silicon

    reg, mat = _create_dummy_mat()
    legendoptics.silicon.pyg4_silicon_attach_complex_rindex(mat, reg)


def test_pyg4_attach_nylon() -> None:
    import legendoptics.nylon

    reg, mat = _create_dummy_mat()
    legendoptics.nylon.pyg4_nylon_attach_rindex(mat, reg)
    legendoptics.nylon.pyg4_nylon_attach_absorption(mat, reg)


def test_pyg4_attach_pen() -> None:
    import legendoptics.pen

    reg, mat = _create_dummy_mat()
    legendoptics.pen.pyg4_pen_attach_rindex(mat, reg)
    legendoptics.pen.pyg4_pen_attach_attenuation(mat, reg)
    legendoptics.pen.pyg4_pen_attach_wls(mat, reg)
    legendoptics.pen.pyg4_pen_attach_scintillation(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pen.pyg4_pen_attach_wls(mat, reg, 0.1)


def test_pyg4_attach_water() -> None:
    import legendoptics.water

    reg, mat = _create_dummy_mat()
    legendoptics.water.pyg4_water_attach_rindex(mat, reg)
    legendoptics.water.pyg4_water_attach_absorption(mat, reg)


def test_pyg4_attach_vm2000() -> None:
    import legendoptics.vm2000

    reg, mat = _create_dummy_mat()
    legendoptics.vm2000.pyg4_vm2000_attach_reflectivity(mat, reg)
    legendoptics.vm2000.pyg4_vm2000_attach_absorption_length(mat, reg)
    legendoptics.vm2000.pyg4_vm2000_attach_particle_scintillationyields(mat, reg)
    legendoptics.vm2000.pyg4_vm2000_attach_wls(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.vm2000.pyg4_vm2000_attach_border_params(mat, reg)


def test_pyg4_attach_pmts() -> None:
    import legendoptics.pmts

    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_acryl_absorption_length(mat, reg)
    legendoptics.pmts.pyg4_pmt_attach_acryl_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_air_absorption_length(mat, reg)
    legendoptics.pmts.pyg4_pmt_attach_air_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_borosilicate_absorption_length(mat, reg)
    legendoptics.pmts.pyg4_pmt_attach_borosilicate_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_photocathode_efficiency(mat, reg)
    legendoptics.pmts.pyg4_pmt_attach_photocathode_reflectivity(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_photocathode_efficiency(mat, reg, name="etl9354")
    legendoptics.pmts.pyg4_pmt_attach_photocathode_reflectivity(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_photocathode_efficiency(mat, reg, name="gerda")
    legendoptics.pmts.pyg4_pmt_attach_photocathode_reflectivity(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_photocathode_efficiency(mat, reg, name="r7081")
    legendoptics.pmts.pyg4_pmt_attach_photocathode_reflectivity(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_photocathode_efficiency(mat, reg, name="l1000")
    legendoptics.pmts.pyg4_pmt_attach_photocathode_reflectivity(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.pmts.pyg4_pmt_attach_steel_efficiency(mat, reg)
    legendoptics.pmts.pyg4_pmt_attach_steel_reflectivity(mat, reg)
