#include <iostream>			// cout
#include <algorithm>		// qsort, swap
#include <random>			// random_device, mt19937, shuffle
#include <string>			// string
#include <vector>			// Vector
#include <mpi.h>			// MPI everything

#define MCW MPI_COMM_WORLD
#define VERBOSE true
#define MASTER 0
/*
Reference: https://www.youtube.com/watch?v=s6FhG--P7z0
*/

/* Main function */
int main(int argc, char *argv[]) {
	// MPI Setup
	int mpi_rank, mpi_size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MCW, &mpi_rank);
	MPI_Comm_size(MCW, &mpi_size);

	// Set and sum size
	int* set;
	int job_size = 10;
	int set_size = mpi_size * job_size;
	int set_sum = set_size * job_size;

	// Allocate new array
	set = new int[set_size];

	// Construct data set to be broadcast
	if (mpi_rank == MASTER) {
		// Random initializer
		std::random_device rd;
		std::mt19937 g(rd());

		// Allocate numbers between [1, set_sum)
		for (int i = 0; i < set_size; ++i) {
			set[i] = g() % (set_sum / 2 - 1) + 1;
		}
	}

	// MPI shares info about the data
	MPI_Bcast(set, set_size, MPI_INT, MASTER, MCW);

	if (VERBOSE && mpi_rank == MASTER) {
		std::string output = "{";
		for (int i = 0; i < set_size; ++i) {
			output += std::to_string(set[i]) + ",";
		}
		std::cout << output << "}" << std::endl;
		std::cout << "Attempting to fill a subset with size " << set_sum << std::endl;
	}

	// Allocate for subset
	bool** subset_grid = new bool*[set_size];
	for (int i = 0; i < set_size; ++i) {
		subset_grid[i] = new bool[set_sum + 1];
		// setup before dynamic programming
		subset_grid[i][0] = true;
	}

	// Start dynamic programming process
	for (int i = 0; i < set_size; ++i) {
		// determine previous solution set

		// 
	}

	// Deallocate variable pointers
	for (int i = 0; i < set_size; ++i) {
		delete subset_grid[i];
	}
	delete subset_grid;
	delete set;

	// Finalize and quit
	MPI_Finalize();
	return 0;
}
