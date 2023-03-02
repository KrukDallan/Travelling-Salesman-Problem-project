#ifndef TSP_H_

#define TSP_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define EPSILON 1e-5

typedef struct
{

	// Basic components
	char input_file[2000];
	char output_file[2000];
	char solution_file[2000];
	char heuristic[100]; // greedy ; extram ; tabu ; VNS ; genetic
	char variant[100]; //2opt ; grasp ; 
	//For debugging pourpose
	int verbose;
	// Time management
	clock_t tstart;
	double timeLimit;
	//"dimension"
	int nnodes;
	//these pointers will point to an (future) allocated memory that will contain x and y coordinates
	double* xcoord;
	double* ycoord;
	double* dist;
	double best_cost;
	double current_cost;
	double current_value;
	int* best_sol;
	int* best_succ;
	int starting_node;

	// Flags useful for algos
	int twoOptFlag;
	int graspFlag;
	int geneticFlag;

	// Greedy
	int* greedySucc;
	double greedyCost;

	// Extra-mileage
	int* emSucc;
	double emCost;

	// Tabu 
	int* tabuSucc;
	int numIter;
	double tabuCost;

	// VNS 
	int* vnsSucc;
	double vnsCost;

	// Genetic
	int* geneticSucc;
	double geneticCost;

	// cplex
	int* cplexSucc;
	double cplexCost;
	int ncols;
	int patchFlag;

} Instance;

typedef struct
{
	int succ_node;
	int index_uncovered_node;
} emstruct;


typedef struct
{
	int* chromosome;
	int startingNode;
	double cost;
	int childSize;
} Population;


#endif