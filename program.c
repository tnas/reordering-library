#include <unistd.h>
#include <ctype.h>
#include "./UnitTests/test_suite_rcm.h"


int main (int argc, char* argv[]){
  
	int opt, num_threads, root;
	float bfs_chunk_size;
	char* matrix_name;
	EXECUTION exec_type;
	int algorithm;
	int* peripheral_nodes;

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
		peripheral_nodes = get_node_peripheral_hsl(matrix_name);
		root = peripheral_nodes[0];
// 		root = get_node_peripheral(matrix_name);
		
		switch (algorithm)
		{
			case serial_rcm : 
				// t = 0
				test_serial_rcm(matrix_name, root);
				break;
				
			case serial_sloan :
				// t = 1
				test_serial_sloan(matrix_name);
				break;
				
			case hsl_rcm : 
				// t = 2
				test_hsl_rcm(matrix_name);
				break;
				
			case hsl_spectral :
				// t = 3
				test_hsl_spectral(matrix_name);
				break;
				
			case hsl_sloan :
				// t = 4
				test_hsl_sloan(matrix_name, peripheral_nodes);
				break;
				
			case unordered_rcm :
				// t = 5
				test_unordered_rcm(matrix_name, num_threads, bfs_chunk_size, root);
				break;
				
			case leveled_rcm :
				// t = 6
				test_leveled_rcm(matrix_name, num_threads, root);
				break;
				
			case bucket_rcm :
				// t = 7
				test_bucket_rcm(matrix_name, num_threads, root);
				break;
				
			case parallel_sloan :
				// t = 8
				test_parallel_sloan(matrix_name, num_threads, peripheral_nodes);
				break;
				
			default :
				printf("*** [Error] Algorithm must be between 0 and 8 ***\n");
				exit(1);
		}
		
		free(peripheral_nodes);
	}
	
	return 0;
}