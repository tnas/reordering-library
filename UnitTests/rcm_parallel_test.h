#include "../CommonFiles/protos_parallel.h"
#include <assert.h>

#define TEST_EXEC_TIMES 5

typedef enum EXECUTION { ALL_TESTS, ONE_INSTANCE} EXECUTION;
typedef struct {
	double time;
	long int bandwidth;
} test_result;

extern void test_prefix_sum();
extern void test_prefix_sum_parallel_serial(int* counts, int max_level);
extern int  get_node_peripheral(const char* path_matrix_file);
test_result run_test_serial_rcm(const char* path_matrix_file, int root);
test_result run_test_unordered_rcm(const char* path_matrix_file, const int num_threads, const float bfs_chunk_size, int root);
test_result run_test_leveled_rcm(const char* path_matrix_file, int root);
extern void run_test_serial_parallel_rcm(const char* path_matrix_file, const int num_threads, const float bfs_chunk_size);
extern void run_all_tests();