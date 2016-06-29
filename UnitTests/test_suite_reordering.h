#include "../CommonFiles/protos_parallel.h"
#include <assert.h>
#include <limits.h>
#include <float.h>

#define TEST_EXEC_TIMES 5

typedef enum {
	serial_rcm,
	serial_sloan,
	hsl_rcm,
	hsl_spectral,
	hsl_sloan,
	unordered_rcm,
	leveled_rcm,
	bucket_rcm,
	parallel_sloan
} reorder_algorithm;


typedef struct {
	double time_peripheral;
	double time_reordering;
	double time_permutation;
	long int original_band;
	long int reorder_band;
	long int wavefront;
	const char* path_matrix_file;
	int root;
	int start_node;
	int end_node;
	char* algorithm_name;
	int num_threads;
	reorder_algorithm algorithm;
	float percent_chunk;
} test;


typedef enum { START, END } PERIPHERAL_NODES;


int  get_node_peripheral     (const char* path_matrix_file);
int* get_node_peripheral_hsl (const char* path_matrix_file);
test test_reorder_algorithm  (test defs);
void run_all_tests           ();