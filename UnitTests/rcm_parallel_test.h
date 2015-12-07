#include "../CommonFiles/protos_parallel.h"
#include <assert.h>

extern void test_prefix_sum();
extern void test_prefix_sum_parallel_serial(int* counts, int max_level);