#include "device.h"
#include "return_codes.h"
#include "test_harness.h"


#define NUM_TESTS 2


int test_get_num_gpus()
{
    int num_devices;
    rc_check(get_num_gpus(&num_devices, 0));
    if (num_devices != 0)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_create_device()
{
    Device_t device;
    rc_check(create_device(&device, NULL));
    if (device != HOST_ONLY)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_get_num_gpus,
            "test_get_num_gpus",
            GRTCODE_SUCCESS
        },
        {
            test_create_device,
            "test_create_device",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_device", NUM_TESTS, tests);
}
