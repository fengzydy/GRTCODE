#include <stdlib.h>
#include "debug.h"
#include "device.h"
#include "floating_point_type.h"
#include "return_codes.h"
#include "spectral_grid.h"
#include "test_harness.h"
#include "utilities.h"


#define NUM_TESTS 5


int test_compare_spectral_grids()
{
    double w0 = 1.;
    double wn = 10000.;
    double dw = 0.1;
    SpectralGrid_t one;
    rc_check(create_spectral_grid(&one, w0, wn, dw));
    SpectralGrid_t two;
    rc_check(create_spectral_grid(&two, w0, wn, dw));
    int result;
    rc_check(compare_spectral_grids(&one, &two, &result));
    if (!result)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_create_spectral_grid()
{
    double w0 = 1.;
    double wn = 10000.;
    double dw = 0.1;
    SpectralGrid_t grid;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    if (grid.w0 != w0 || grid.wn != wn || grid.dw != dw || grid.n != 100000)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_grid_point_index()
{
    double w0 = 1.;
    double wn = 10000.;
    double dw = 0.1;
    SpectralGrid_t grid;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    double w = 5555.;
    uint64_t index;
    rc_check(grid_point_index(grid, w, &index));
    if (index != (uint64_t)(((int)(1./dw))*(w - w0)))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_grid_points()
{
    double w0 = 1.;
    double wn = 10000.;
    double dw = 0.1;
    SpectralGrid_t grid;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    Device_t device;
    rc_check(create_device(&device, NULL));
    fp_t * buffer;
    rc_check(grid_points(grid, &buffer, device));
    int returncode;
    if (buffer[0] == w0 && buffer[grid.n - 1] == wn)
    {
        returncode = GRTCODE_SUCCESS;
    }
    else
    {
        returncode = GRTCODE_VALUE_ERR;
    }
    gfree(buffer, device);
    return returncode;
}


int test_interpolate_to_grid()
{
    double w0 = 1.;
    double wn = 10000.;
    double dw = 0.1;
    SpectralGrid_t grid;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));

    size_t n = 5;
    fp_t x[5];
    fp_t y[5];
    int i;
    for (i=0; i<n; ++i)
    {
        x[i] = w0 + (fp_t)i + 2.;
        y[i] = (fp_t)i;
    }
    fp_t * newy = (fp_t *)malloc(sizeof(*newy)*grid.n);
    Sample1d_t interp = linear_sample;
    Sample1d_t extrap = constant_extrapolation;
    rc_check(interpolate_to_grid(grid, x, y, n, newy, interp, extrap));
    free(newy);
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_compare_spectral_grids,
            "test_compare_spectral_grids",
            GRTCODE_SUCCESS
        },
        {
            test_create_spectral_grid,
            "test_create_spectral_grid",
            GRTCODE_SUCCESS
        },
        {
            test_grid_point_index,
            "test_grid_point_index",
            GRTCODE_SUCCESS
        },
        {
            test_grid_points,
            "test_grid_points",
            GRTCODE_SUCCESS
        },
        {
            test_interpolate_to_grid,
            "test_interpolate_to_grid",
            GRTCODE_SUCCESS
        },
    };
    return test_harness("test_spectral_grid", NUM_TESTS, tests);
}
