#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "floating_point_type.h"
#include "test_harness.h"


/*Check values in an array.*/
int check_array(fp_t const * const actual, fp_t const * const expected, size_t const size,
                fp_t const tolerance, int const absolute)
{
    int failures = 0;
    size_t i;
    for (i=0; i<size; ++i)
    {
        fp_t error;
        if (absolute)
        {
            error = actual[i] - expected[i];
        }
        else
        {
            error = 100.*(actual[i] - expected[i])/expected[i];
        }
        if (fabs(error) > tolerance)
        {
            failures++;
        }
    }
    return failures;
}


/*Testing harness.*/
int test_harness(char const * const title, int const num_tests, Test_t const * const tests)
{
    printf("Running %s tests.\n", title);
    int passed = 0;
    int i;
    for (i=0; i<num_tests; ++i)
    {
        printf("Running test %s:", tests[i].name);
        if (tests[i].test_function_pointer() == tests[i].returncode)
        {
            printf(" passed.\n");
            passed++;
        }
        else
        {
            printf(" failed.\n");
        }
    }
    printf("%d/%d tests passed.\n", passed, num_tests);
    if (passed == num_tests)
    {
        return EXIT_SUCCESS;
    }
    else
    {
        return EXIT_FAILURE;
    }
}
