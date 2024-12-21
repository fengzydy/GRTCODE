#include <stdint.h>
#include "grtcode_utilities.h"
#include "kernels.h"
#include "kernel_utils.h"
#include "line_shape.h"
#include "molecules.h"
#include "RFM_voigt.h"
#include "spectral_bin.h"
#include "spectral_bin-internal.h"
#include "test_harness.h"
#include "tips2017.h"


#define NUM_TESTS 15


int test_calc_line_centers()
{
    uint64_t num_lines;
    int num_layers;
    fp_t v0[5];
    fp_t delta[5];
    fp_t pressure[1];
    fp_t vnn[5];
    rc_check(calc_line_centers(num_lines, num_layers, v0, delta, pressure, vnn));
    return GRTCODE_SUCCESS;
}


int test_calc_partition_functions()
{
    int num_layers;
    int mol_id = H2O;
    int num_iso = 1;
    fp_t temperature[1];
    fp_t q[1];
    rc_check(calc_partition_functions(num_layers, mol_id, num_iso, temperature, q));
    return GRTCODE_SUCCESS;
}


int test_calc_line_strengths()
{
    uint64_t num_lines;
    int num_layers;
    int num_iso;
    int iso[10];
    fp_t s0[10];
    fp_t vnn[10];
    fp_t en[10];
    fp_t temperature[1];
    fp_t q[1];
    fp_t snn[10];
    rc_check(calc_line_strengths(num_lines, num_layers, num_iso, iso, s0, vnn, en,
                                 temperature, q, snn));
    return GRTCODE_SUCCESS;
}


int test_calc_lorentz_hw()
{
    uint64_t num_lines;
    int num_layers;
    fp_t n[10];
    fp_t yair[10];
    fp_t yself[10];
    fp_t temperature[1];
    fp_t pressure[1];
    fp_t ps[1];
    fp_t gamma[10];
    rc_check(calc_lorentz_hw(num_lines, num_layers, n, yair, yself, temperature,
                             pressure, ps, gamma));
    return GRTCODE_SUCCESS;
}


int test_calc_doppler_hw()
{
    uint64_t num_lines;
    int num_layers;
    fp_t mass;
    fp_t vnn[10];
    fp_t temperature[1];
    fp_t alpha[10];
    rc_check(calc_doppler_hw(num_lines, num_layers, mass, vnn, temperature, alpha));
    return GRTCODE_SUCCESS;
}


int test_sort_lines()
{
    uint64_t num_lines;
    int num_layers;
    fp_t vnn[10];
    fp_t snn[10];
    fp_t gamma[10];
    fp_t alpha[10];
    rc_check(sort_lines(num_lines, num_layers, vnn, snn, gamma, alpha));
    return GRTCODE_SUCCESS;
}


int test_calc_optical_depth_bin_sweep()
{
    uint64_t num_lines;
    int num_layers;
    fp_t vnn[10];
    fp_t snn[10];
    fp_t gamma[10];
    fp_t alpha[10];
    fp_t n[10];
    SpectralBins_t bins;
    fp_t tau[10];
    rc_check(calc_optical_depth_bin_sweep(num_lines, num_layers, vnn, snn, gamma,
                                          alpha, n, bins, tau));
    return GRTCODE_SUCCESS;
}


int test_calc_optical_depth_line_sweep()
{
    uint64_t num_lines;
    int num_layers;
    fp_t vnn[10];
    fp_t snn[10];
    fp_t gamma[10];
    fp_t alpha[10];
    fp_t n[10];
    SpectralBins_t bins;
    fp_t tau[10];
    rc_check(calc_optical_depth_line_sweep(num_lines, num_layers, vnn, snn, gamma,
                                           alpha, n, bins, tau));
    return GRTCODE_SUCCESS;
}


int test_calc_optical_depth_line_sample()
{
    uint64_t num_lines;
    int num_layers;
    fp_t vnn[10];
    fp_t snn[10];
    fp_t gamma[10];
    fp_t alpha[10];
    fp_t n[10];
    SpectralBins_t bins;
    fp_t tau[10];
    rc_check(calc_optical_depth_line_sample(num_lines, num_layers, vnn, snn, gamma,
                                            alpha, n, bins, tau, NULL, NULL));
    return GRTCODE_SUCCESS;
}


int test_calc_water_vapor_ctm_optical_depth()
{
    uint64_t num_wpoints;
    int num_layers;
    fp_t tau[60];
    fp_t CS[50];
    fp_t temperature[5];
    fp_t ps[5];
    fp_t n[5];
    fp_t T0[60];
    fp_t CF[60];
    fp_t pressure[5];
    fp_t T0F[50];
    rc_check(calc_water_vapor_ctm_optical_depth(num_wpoints, num_layers,
                                                tau, CS, temperature, ps,
                                                n, T0, CF, pressure, T0F));
    return GRTCODE_SUCCESS;
}


int test_calc_ozone_ctm_optical_depth()
{
    uint64_t num_wpoints;
    int num_layers;
    fp_t cross_section[50];
    fp_t n[5];
    fp_t tau[100];
    rc_check(calc_ozone_ctm_optical_depth(num_wpoints, num_layers, cross_section,
                                          n, tau));
    return GRTCODE_SUCCESS;
}


int test_interpolate()
{
    SpectralBins_t bins;
    fp_t tau[50];
    rc_check(interpolate(bins, tau));
    return GRTCODE_SUCCESS;
}


int test_interpolate_last_bin()
{
    SpectralBins_t bins;
    fp_t tau[1000];
    rc_check(interpolate_last_bin(bins, tau));
    return GRTCODE_SUCCESS;
}


int test_calc_cfc_optical_depth()
{
    uint64_t num_wpoints = 1000;
    int num_layers = 5;
    fp_t n[5];
    fp_t x[5];
    fp_t cross_section[1000];
    fp_t tau[1000];
    rc_check(calc_cfc_optical_depth(num_wpoints, num_layers, n, x,
                                    cross_section, tau));
    return GRTCODE_SUCCESS;
}


int test_calc_cia_optical_depth()
{
    uint64_t num_wpoints = 1000;
    int num_layers = 5;
    fp_t pressure[5];
    fp_t temperature[5];
    fp_t x1[5];
    fp_t x2[5];
    fp_t cross_section[1000];
    fp_t tau[1000];
    rc_check(calc_cia_optical_depth(num_wpoints, num_layers, pressure,
                                    temperature, x1, x2, cross_section, tau));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_calc_line_centers,
            "test_calc_line_centers",
            GRTCODE_SUCCESS
        },
        {
            test_calc_partition_functions,
            "test_calc_partition_functions",
            GRTCODE_SUCCESS
        },
        {
            test_calc_line_strengths,
            "test_calc_line_strengths",
            GRTCODE_SUCCESS
        },
        {
            test_calc_lorentz_hw,
            "test_calc_lorentz_hw",
            GRTCODE_SUCCESS
        },
        {
            test_calc_doppler_hw,
            "test_calc_doppler_hw",
            GRTCODE_SUCCESS
        },
        {
            test_sort_lines,
            "test_sort_lines",
            GRTCODE_SUCCESS
        },
        {
            test_calc_optical_depth_bin_sweep,
            "test_calc_optical_depth_bin_sweep",
            GRTCODE_SUCCESS
        },
        {
            test_calc_optical_depth_line_sweep,
            "test_calc_optical_depth_line_sweep",
            GRTCODE_SUCCESS
        },
        {
            test_calc_optical_depth_line_sample,
            "test_calc_optical_depth_line_sample",
            GRTCODE_SUCCESS
        },
        {
            test_calc_water_vapor_ctm_optical_depth,
            "test_calc_water_vapor_ctm_optical_depth",
            GRTCODE_SUCCESS
        },
        {
            test_calc_ozone_ctm_optical_depth,
            "test_calc_ozone_ctm_optical_depth",
            GRTCODE_SUCCESS
        },
        {
            test_interpolate,
            "test_interpolate",
            GRTCODE_SUCCESS
        },
        {
            test_interpolate_last_bin,
            "test_interpolate_last_bin",
            GRTCODE_SUCCESS
        },
        {
            test_calc_cfc_optical_depth,
            "test_calc_cfc_optical_depth",
            GRTCODE_SUCCESS
        },
        {
            test_calc_cia_optical_depth,
            "test_calc_cia_optical_depth",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_kernels", NUM_TESTS, tests);
}
