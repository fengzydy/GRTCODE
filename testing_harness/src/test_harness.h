#ifndef TEST_HARNESS_H_
#define TEST_HARNESS_H_


#include <stdlib.h>

#include "floating_point_type.h"


#define rc_check(function) { \
    int rc = function; \
    if (rc != GRTCODE_SUCCESS) \
    {  \
        return GRTCODE_VALUE_ERR; \
    } \
}


/*Check values in an array.*/
int check_array(
    fp_t const * const actual, /*Array of actual values.*/
    fp_t const * const expected, /*Array of expected values.*/
    size_t const size, /*Size of the input arrays.*/
    fp_t const tolerance, /*Tolerance.*/
    int const absolute /*Use absolute differences or percent differences?*/
);


/*Helper structure so tests can be looped through.*/
typedef struct Test
{
    int (* test_function_pointer)(void);
    char * name;
    int returncode;
} Test_t;


/*Testing harness.*/
int test_harness(
    char const * const title, /*Name of the test suite.*/
    int const num_tests, /*Size of the input tests array.*/
    Test_t const * const tests /*Array of tests.*/
);


#endif
