#include "circ1.h"
#include "curtis_godson.h"
#include "floating_point_type.h"
#include "test_harness.h"
#include "return_codes.h"


#define NUM_TESTS 3


/*Calculate integrated number densities.*/
int test_calc_number_densities()
{
    fp_t p[NUM_LEVELS];
    int i;
    for (i=0; i<num_levels; ++i)
    {
        fp_t const mbtoatm = 0.000986923;
        p[i] = level_pressure[i]*mbtoatm;
    }
    fp_t n[NUM_LAYERS];
    rc_check(calc_number_densities(num_layers, p, n));

    fp_t n_ref[NUM_LAYERS] = {
    };
    rc_check(check_array(n_ref, n, num_layers, 1.e-7, 1));
    return GRTCODE_SUCCESS;
}


/*Calculate layer pressures and temperatures.*/
int test_calc_pressures_and_temperatures()
{
    fp_t pavg[NUM_LAYERS];
    fp_t tavg[NUM_LAYERS];
    rc_check(calc_pressures_and_temperatures(num_layers, level_pressure, level_temperature,
                                             pavg, tavg));

    fp_t pavg_ref[NUM_LAYERS] = {
    };
    fp_t tavg_ref[NUM_LAYERS] = {
    };
    rc_check(check_array(pavg_ref, pavg, num_layers, 1.e-7, 1));
    rc_check(check_array(tavg_ref, tavg, num_layers, 1.e-7, 1));
    return GRTCODE_SUCCESS;
}


/*Calculate partial pressures and number densities.*/
int test_calc_partial_pressures_and_number_densities()
{
    fp_t abundance[NUM_LEVELS];
    fp_t number_density[NUM_LAYERS];
    fp_t ps[NUM_LAYERS];
    fp_t ns[NUM_LAYERS];
    rc_check(calc_partial_pressures_and_number_densities(num_layers, level_pressure,
                                                         abundance, number_density,
                                                         ps, ns));

    fp_t ps_ref[NUM_LAYERS];
    fp_t ns_ref[NUM_LAYERS];
    rc_check(check_array(ps_ref, ps, num_layers, 1.e-7, 1));
    rc_check(check_array(ns_ref, ns, num_layers, 1.e-7, 1));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_calc_number_densities,
            "test_calc_number_densities",
            GRTCODE_SUCCESS
        },
        {
            test_calc_pressures_and_temperatures,
            "test_calc_pressures_and_temperatures",
            GRTCODE_SUCCESS
        },
        {
            test_calc_partial_pressures_and_number_densities,
            "test_calc_partial_pressures_and_number_densities",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_curtis_godson", NUM_TESTS, tests);
}
