#include <omp.h>
#include <iostream>

int main () {
	#pragma omp parallel
	{
		#pragma omp master
		{
			std::cout << "Processors available: " << omp_get_num_procs() << std::endl;
			std::cout << "Threads being used: " << omp_get_num_threads() << std::endl;
			std::cout << "Max threads available: " << omp_get_max_threads() << std::endl;
			std::cout << "Parallel region? " << omp_in_parallel() << std::endl;
			std::cout << "Dynamic threads enabled? " << omp_get_dynamic() << std::endl;
			std::cout << "Nested parallelism supported? " << omp_get_nested() << std::endl;
		}
	}
	return 0;
}
