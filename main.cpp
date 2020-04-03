#include <iostream>			// cout
#include <algorithm>		// qsort, swap
#include <random>			// random_device, mt19937, shuffle
#include <time.h>
#include <string>			// string
#include <vector>			// Vector
#include <mpi.h>			// MPI everything

#define MCW MPI_COMM_WORLD
#define VERBOSE false
#define MASTER 0
/*
Reference: https://www.youtube.com/watch?v=s6FhG--P7z0
*/

std::vector<int> traceback(int** cache, int* set, int row, int col) {
    std::vector<int> subset;
    while(row > 0 && col > 0) {
        int num = set[row];
        if(cache[row][col] == 1) {
            row--;
            if(cache[row][col] == 0){
                subset.insert(subset.begin(), num);
                col -= num;
            }
        } else 
            throw -1;
    }
    if(col > 0)
        subset.insert(subset.begin(), set[row]);
    return subset;
}

/* Main function */
int main(int argc, char *argv[]) {
	// MPI Setup
	int mpi_rank, mpi_size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MCW, &mpi_rank);
	MPI_Comm_size(MCW, &mpi_size);

    // Random initializer
	srand(time(0));

	// Set and sum size
	int* set;
	int job_size = pow(2, 2);
	int set_size = mpi_size * job_size; // N
	int target;
    int row_size;

	// Allocate new array
	set = new int[set_size];

	// Construct data set to be broadcast
	if (mpi_rank == MASTER) {
        // Allocate numbers between [1, 2N]
		for (int i = 0; i < set_size; ++i) {
			set[i] = rand() % (2*set_size) + 1;
		}
        // Initialize Target
        int tmp = pow(set_size, 2) / 2;
        target = ((((rand() % tmp) + tmp)/mpi_size)*mpi_size)-1;
	}

	// MPI shares info about the data
	MPI_Bcast(set, set_size, MPI_INT, MASTER, MCW);
    MPI_Bcast(&target, 1, MPI_INT, MASTER, MCW);
    row_size = target + 1;

	int total = 0;
	if (mpi_rank == MASTER) {
		std::string output = "{";
		for (int i = 0; i < set_size; ++i) {
			output += std::to_string(set[i]) + ",";
			total += set[i];
		}
		std::cout << output << "}" << std::endl;
		std::cout << "Attempting to find a subset with target sum " << target << std::endl;
		std::cout << "Total Set Sum: " << total << std::endl;
	}

	// Allocate for subset Initialize Cache
	int** subset_grid = new int*[set_size]; // Cache
	for (int i = 0; i < set_size; ++i) {
		subset_grid[i] = new int[row_size];
		// setup before dynamic programming
        for(int j = 0; j < row_size; j++)
            subset_grid[i][j] = false;
		subset_grid[i][0] = true;
	}
    subset_grid[0][set[0]] = true;

	// Start dynamic programming process
    int portion_size = (row_size)/mpi_size;    
    int* row_portion = new int[portion_size];
	for (int i = 1; i < set_size; ++i) {
        if(VERBOSE && mpi_rank == MASTER)
            std::cout << "Processing Row " << i << " out of " << set_size << std::endl;
		int new_num = set[i];
        int offset = mpi_rank * portion_size;
        for(int j = 0; j < portion_size; j++){
            if(VERBOSE && mpi_rank == MASTER)
                std::cout << "Processing Col " << j << " out of " << portion_size << std::endl;
            if(new_num > j + offset)
                row_portion[j] = subset_grid[i-1][j+offset];
            else 
                row_portion[j] = subset_grid[i-1][j+offset] || subset_grid[i-1][j+offset-new_num];
        }

        // Share data between Processes
        if(VERBOSE && mpi_rank == MASTER)
            std::cout << "\nSending Data\n";
        int* new_row = new int[row_size];
        for(int dest = 0; dest < mpi_size; dest++)
            MPI_Send(row_portion, portion_size, MPI_INT, dest, i, MCW);

        MPI_Barrier(MCW);

        if(VERBOSE && mpi_rank == MASTER)
            std::cout << "\nReceiving Data\n";
        for(int source = 0; source < mpi_size; source++) {
            int* received = new int[portion_size];
            MPI_Recv(received, portion_size, MPI_INT, source, i, MCW, MPI_STATUS_IGNORE);
            for(int j = 0; j < portion_size; j++) {
                subset_grid[i][j+(source*portion_size)] = received[j];
            }
        }
	}

    // // Print Cache
    // if(!mpi_rank) {   
    //     std::cout << std::endl; 
    //     for(int i = 0; i < set_size; i++) {
    //         for(int j = 0; j < row_size; j++){
    //             std::cout << subset_grid[i][j] << ' ';
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    if(mpi_rank == MASTER){
    	std::vector<int> subset;
        if(subset_grid[set_size-1][target]) {
            subset = traceback(subset_grid, set, set_size-1, target);
        } else {
            // THERE WAS NO COMBINATION TO EQUAL THE TARGET SUM
            // EXECUTE "CLOSEST SUM" EXTENSION
            // FEEL FREE TO USE TRACEBACK TO SHOW THE CLOSEST SUM SUBSET
            for (int x = target; x > 0; x--) {
            	if (subset_grid[set_size-1][x]) {
            		subset = traceback(subset_grid, set, set_size-1, x);
            		break;
            	}
            }
        }
    	// Print things
        int total = 0;
        std::cout << "\nThere is a subset which sums to " << target << ".\n";
        std::string output = "{";
        for (int i = 0; i < subset.size(); i++) {
            output += std::to_string(subset[i]) + ", ";
            total += subset[i];
        }
        std::cout << output << "}" << std::endl;
        if (total == target) {
	        std::cout << "Verified total is " << total << std::endl << std::endl;
	    }
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
