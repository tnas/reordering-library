#include <unistd.h>
#include <ctype.h>
#include "./UnitTests/rcm_parallel_test.h"


int main (int argc, char* argv[]){
  
	int opt, num_threads;
	float bfs_chunk_size;
	char* matrix_name;
	EXECUTION exec_type;
// 	int algorithm;

	while ((opt = getopt(argc, argv, "m:p:s:t:a")) != -1)
	{
		switch (opt)
		{
			case 'm' :
				matrix_name = optarg;
				exec_type = ONE_INSTANCE;
				break;
/*			case 't' :
				algorithm = optarg;
				break;*/				
			case 'p':
				num_threads = atoi(optarg);
				break;
				
			case 's':
				bfs_chunk_size = atof(optarg);
				break;
				
			case 'a' :
				exec_type = ALL_TESTS;
				break;
		}
	}
	
	
	int root = get_node_peripheral(matrix_name);
// 	int root = 1;
	run_test_leveled_rcm(matrix_name, num_threads, root);
	run_test_serial_rcm(matrix_name, root);
	
// 	if (exec_type == ALL_TESTS)
// 	{
// 		run_all_tests();
// 	}
// 	else 
// 	{
// 		run_test_serial_parallel_rcm(matrix_name, num_threads, bfs_chunk_size);
// 	}
	
	return 0;
}