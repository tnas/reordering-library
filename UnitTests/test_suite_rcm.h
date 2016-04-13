#include "../CommonFiles/protos_parallel.h"
#include <assert.h>
#include <limits.h>
#include <float.h>

#define TEST_EXEC_TIMES 5

typedef enum {
	serial,
	parallel
} reordering_strategy;

typedef enum {
	serial_rcm,
	hsl_rcm,
	hsl_spectral,
	serial_sloan,
	unordered_rcm,
	leveled_rcm,
	bucket_rcm,
	parallel_sloan
} reorder_algorithm;


typedef struct {
	double time;
	long int bandwidth;
	const char* path_matrix_file;
	int root;
	int start_node;
	int end_node;
	char* algorithm_name;
	int num_threads;
	reorder_algorithm algorithm;
	reordering_strategy strategy;
	float percent_chunk;
} test_def;


typedef enum EXECUTION { ALL_TESTS, ONE_INSTANCE } EXECUTION;



extern void test_prefix_sum();
extern void test_prefix_sum_parallel_serial(int* counts, int max_level);
int get_node_peripheral(const char* path_matrix_file);

test_def test_serial_rcm(const char* path_matrix_file, int root);
test_def test_hsl_spectral(const char* path_matrix_file);
test_def test_hsl_rcm(const char* path_matrix_file);
test_def test_serial_sloan(const char* path_matrix_file);

test_def test_unordered_rcm(const char* path_matrix_file, const int num_threads, const float bfs_chunk_percent, int root);
test_def test_leveled_rcm(const char* path_matrix_file, const int num_threads, int root);
test_def test_bucket_rcm(const char* path_matrix_file, const int num_threads, int root);
test_def test_parallel_sloan(const char* path_matrix_file, const int num_threads);

extern void run_all_tests();