#include "cfcs.h"
#include "collision_induced_absorption.h"
#include "gas_optics.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "spectral_bin.h"
#include "test_harness.h"
#include "tips2017.h"
#include "water_vapor_continuum.h"


#define NUM_TESTS 10


int test_create_gas_optics()
{
    GasOptics_t gas_optics;
    int num_levels = 7;
    SpectralGrid_t grid;
    Device_t device;
    char * hitran_path = "";
    char * h2o_ctm_dir = "";
    char * o3_ctm_file = "";
    double wcutoff = 25.;
    int optical_depth_method = line_sample;
    rc_check(create_gas_optics(&gas_optics, num_levels, &grid, &device, hitran_path,
                               h2o_ctm_dir, o3_ctm_file, &wcutoff, &optical_depth_method));
    return GRTCODE_SUCCESS;
}


int test_destroy_gas_optics()
{
    GasOptics_t gas_optics;
    rc_check(destroy_gas_optics(&gas_optics));
    return GRTCODE_SUCCESS;
}


int test_add_molecule()
{
    GasOptics_t gas_optics;
    int molecule_id = H2O;
    double min_line_center = 1.;
    double max_line_center = 1000.;
    rc_check(add_molecule(&gas_optics, molecule_id, &min_line_center, &max_line_center));
    return GRTCODE_SUCCESS;
}


int test_set_molecule_ppmv()
{
    GasOptics_t gas_optics;
    int molecule_id = CO2;
    fp_t ppmv[1] = {400.};
    rc_check(set_molecule_ppmv(&gas_optics, molecule_id, ppmv));
    return GRTCODE_SUCCESS;
}


int test_add_cfc()
{
    GasOptics_t gas_optics;
    int cfc_id = CFC11;
    char * filepath = "";
    rc_check(add_cfc(&gas_optics, cfc_id, filepath));
    return GRTCODE_SUCCESS;
}


int test_set_cfc_ppmv()
{
    GasOptics_t gas_optics;
    int cfc_id = CFC11;
    fp_t ppmv[1] = {10.};
    rc_check(set_cfc_ppmv(&gas_optics, cfc_id, ppmv));
    return GRTCODE_SUCCESS;
}


int test_add_cia()
{
    GasOptics_t gas_optics;
    int species1 = N2;
    int species2 = N2;
    char * filepath = "";
    rc_check(add_cia(&gas_optics, species1, species2, filepath));
    return GRTCODE_SUCCESS;
}


int test_set_cia_ppmv()
{
    GasOptics_t gas_optics;
    int cia_id = N2;
    fp_t ppmv[1] = {10.};
    rc_check(set_cia_ppmv(&gas_optics, cia_id, ppmv));
    return GRTCODE_SUCCESS;
}


int test_calculate_optical_depth()
{
    GasOptics_t gas_optics;
    fp_t pressure[1] = {700.};
    fp_t temperature[1] = {275.};
    Optics_t optics;
    rc_check(calculate_optical_depth(&gas_optics, pressure, temperature, &optics));
    return GRTCODE_SUCCESS;
}


int test_get_num_molecules()
{
    GasOptics_t gas_optics;
    int n;
    rc_check(get_num_molecules(&gas_optics, &n));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_create_gas_optics,
            "test_create_gas_optics",
            GRTCODE_SUCCESS
        },
        {
            test_destroy_gas_optics,
            "test_destroy_gas_optics",
            GRTCODE_SUCCESS
        },
        {
            test_add_molecule,
            "test_add_molecule",
            GRTCODE_SUCCESS
        },
        {
            test_set_molecule_ppmv,
            "test_set_molecule_ppmv",
            GRTCODE_SUCCESS
        },
        {
            test_add_cfc,
            "test_add_cfc",
            GRTCODE_SUCCESS
        },
        {
            test_set_cfc_ppmv,
            "test_set_cfc_ppmv",
            GRTCODE_SUCCESS
        },
        {
            test_add_cia,
            "test_add_cia",
            GRTCODE_SUCCESS
        },
        {
            test_set_cia_ppmv,
            "test_set_cia_ppmv",
            GRTCODE_SUCCESS
        },
        {
            test_calculate_optical_depth,
            "test_calculate_optical_depth",
            GRTCODE_SUCCESS
        },
        {
            test_get_num_molecules,
            "test_get_num_molecules",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_gas_optics", NUM_TESTS, tests);
}
