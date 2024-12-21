#include "return_codes.h"
#include "parse_csv.h"
#include "test_harness.h"


#define NUM_TESTS 1


/*Parse a csv file, assuming that the first line in the file contains
  headers for each of the columns.*/
int test_parse_csv()
{
    char * path = "";
    int num_lines;
    int num_columns;
    int ignore_headers = 0;
    char ** data;
    rc_check(parse_csv(path, &num_lines, &num_columns, ignore_headers, &data));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_parse_csv,
            "test_parse_csv",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_curtis_godson", NUM_TESTS, tests);
}
