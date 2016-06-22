#include <unistd.h>
#include <ctype.h>
#include "./UnitTests/test_suite_reordering.h"


int main (int argc, char* argv[]){
  
	int opt, num_threads, algorithm;
	float bfs_chunk_size;
	char* matrix_name;
	EXECUTION exec_type;
	test defs;

	while ((opt = getopt(argc, argv, "m:p:b:t:a")) != -1)
	{
		switch (opt)
		{
			case 'm' :
				matrix_name = optarg;
				exec_type = ONE_INSTANCE;
				break;
			case 't' :
				algorithm = atoi(optarg);
				break;		
				
			case 'p':
				num_threads = atoi(optarg);
				break;
				
			case 'b':
				bfs_chunk_size = atof(optarg);
				break;
				
			case 'a' :
				exec_type = ALL_TESTS;
				break;
		}
	}
	
	if (exec_type == ALL_TESTS)
	{
		run_all_tests();
	}
	else 
	{
		defs.path_matrix_file = matrix_name;
		defs.algorithm        = algorithm;
		defs.percent_chunk    = bfs_chunk_size;
		defs.num_threads      = num_threads;
		
		test_reorder_algorithm(defs);
	}
	
	return 0;
}