#include <iostream>			// cout
#include <algorithm>		// qsort, swap
#include <random>			// random_device, mt19937, shuffle
#include <string>			// striog
#include <vector>			// Vector
#include <mpi.h>			// MPI everything

#define MCW MPI_COMM_WORLD
#define VERBOSE true
#define MASTER 0


int main(int argc, char *argv[]) {
	// MPI Setup
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MCW, &rank);
	MPI_Comm_size(MCW, &size);

	// Set and sum size
	int* set;
	int set_size = 100, set_sum = 100;
	set = new int[set_size];

	// Construct data set to be broadcast
	if (rank == MASTER) {
		// Random initializer
		std::random_device rd;
		std::mt19937 g(rd());


		// Allocate numbers between [1, set_sum)
		for (int i = 0; i < set_size; ++i) {
			set[i] = g() % (set_sum - 1) + 1;
		}
	}

	// MPI shares info about the data
	MPI_Bcast(set, set_size, MPI_INT, MASTER, MCW);

	if (VERBOSE) {
		std::string output = "[" + std::to_string(rank) + "] = {";
		for (int i = 0; i < set_sum; ++i) {
			output += std::to_string(set[i]) + ",";
		}
		std::cout << output << "}" << std::endl;
	}

	// Allocate for subset
	bool** subset_grid = new bool*[set_size];
	for (int i = 0; i < set_size; ++i) {
		subset_grid[i] = new bool[set_sum + 1];
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
