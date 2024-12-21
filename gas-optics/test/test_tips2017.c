#include "grtcode_utilities.h"
#include "molecules.h"
#include "test_harness.h"
#include "tips2017.h"


#define NUM_TESTS 6


int test_inittips_d()
{
    rc_check(inittips_d());
    return GRTCODE_SUCCESS;
}


int test_Q_H2O()
{
    fp_t temperature = 275;
    int iso = 1;
    fp_t qt = Q(H2O, temperature, iso);
    return GRTCODE_SUCCESS;
}


int test_Q_CO2()
{
    fp_t temperature = 275;
    int iso = 1;
    fp_t qt = Q(CO2, temperature, iso);
    return GRTCODE_SUCCESS;
}


int test_Q_CH4()
{
    fp_t temperature = 275;
    int iso = 1;
    fp_t qt = Q(CH4, temperature, iso);
    return GRTCODE_SUCCESS;
}


int test_Q_N2O()
{
    fp_t temperature = 275;
    int iso = 1;
    fp_t qt = Q(N2O, temperature, iso);
    return GRTCODE_SUCCESS;
}


int test_Q_O3()
{
    fp_t temperature = 275;
    int iso = 1;
    fp_t qt = Q(O3, temperature, iso);
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_inittips_d,
            "test_inittips_d",
            GRTCODE_SUCCESS
        },
        {
            test_Q_H2O,
            "test_Q_H2O",
            GRTCODE_SUCCESS
        },
        {
            test_Q_CO2,
            "test_Q_CO2",
            GRTCODE_SUCCESS
        },
        {
            test_Q_CH4,
            "test_Q_CH4",
            GRTCODE_SUCCESS
        },
        {
            test_Q_N2O,
            "test_Q_N2O",
            GRTCODE_SUCCESS
        },
        {
            test_Q_O3,
            "test_Q_O3",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_tips2017", NUM_TESTS, tests);
}
