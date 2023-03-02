#include <time.h>
#include <errno.h>
#include "tsp.h"

#include "utilities.h"
#include "constructing_heuristics.h"
#include "meta-heuristics.h"

int main(int argc, char** argv)
{
	// variable used to store "parse_command_line" return value
	int pcl = 0;
	Instance instance;
	instance.verbose = 0;
	// Start taking track of "time"
	instance.tstart = clock();
	// randomizer
	srand((unsigned)time(NULL));

	// if pcl == 0 then generate random coordinates, else the file must have specified points coordinates
	pcl = parse_command_line(argc, argv, &instance);
	if (pcl)
	{
		read_input(&instance);
	}
	else
	{
		generate_random_points(&instance, instance.nnodes);
		strcat(strcat(strcat(strcat(instance.solution_file, "_"), instance.heuristic), instance.variant), ".dat");
	}

	for (int i = 0; i < instance.nnodes; i++)
	{
		for (int j = 0; j < instance.nnodes; j++)
		{
			distance(i, j, &instance);
		}
	}
	//----------------------------
	// Initial values for some variables of "instance"
	instance.best_sol = NULL;
	instance.best_succ = NULL;
	instance.greedySucc = NULL;
	instance.emSucc = NULL;
	instance.tabuSucc = NULL;
	instance.vnsSucc = NULL;
	instance.geneticSucc = NULL;
	instance.cplexSucc = NULL;

	instance.greedyCost = DBL_MAX;
	instance.emCost = DBL_MAX;
	instance.geneticCost = DBL_MAX;
	instance.tabuCost = DBL_MAX;
	instance.vnsCost = DBL_MAX;
	instance.cplexCost = DBL_MAX;
	//------------------------------
	if (strcmp(instance.heuristic, "greedy") == 0)
	{
		instance.greedySucc = (int*)calloc(instance.nnodes, sizeof(int));
		greed_search(&instance, 0, 1, NULL, 0, NULL);
		clock_t end = clock();
		printf("\nGreedy solution cost: %f\n", instance.greedyCost);
		printf("\nElapsed time: %.3f\n", (double)((end - instance.tstart)) / CLOCKS_PER_SEC);
	}
	if (strcmp(instance.heuristic, "extram") == 0)
	{
		instance.emSucc = (int*)calloc(instance.nnodes, sizeof(int));
		extra_mileage(&instance, NULL, NULL, NULL);
		clock_t end = clock();
		printf("\nExtra-mileage solution cost: %f\n", instance.emCost);
		printf("\nElapsed time: %.3f\n", (double)((end - instance.tstart)) / CLOCKS_PER_SEC);
	}
	if (strcmp(instance.heuristic, "tabu") == 0)
	{
		instance.tabuSucc = (int*)calloc(instance.nnodes, sizeof(int));
		tabuSearch(&instance);
		clock_t end = clock();
		printf("\nTabu solution cost: %f\n", instance.tabuCost);
		printf("\nElapsed time: %.3f\n", (double)(end - instance.tstart) / CLOCKS_PER_SEC);
	}
	if (strcmp(instance.heuristic, "vns") == 0)
	{
		instance.vnsSucc = (int*)calloc(instance.nnodes, sizeof(int));
		VNS(&instance);
		clock_t end = clock();
		printf("\nVNS solution cost: %f\n", instance.vnsCost);
		printf("\nElapsed time: %.3f\n", (double)((end - instance.tstart)) / CLOCKS_PER_SEC);
	}
	if (strcmp(instance.heuristic, "genetic") == 0)
	{
		instance.geneticSucc = (int*)calloc(instance.nnodes, sizeof(int));
		genetic(&instance);
		clock_t end = clock();
		printf("\nGenetic solution cost: %f\n", instance.geneticCost);
		printf("\nElapsed time: %.3f\n", (double)((end - instance.tstart)) / CLOCKS_PER_SEC);
	}
	if (strcmp(instance.heuristic, "cplex") == 0)
	{
		instance.cplexSucc = (int*)calloc(instance.nnodes, sizeof(int));
		instance.patchFlag = 0;
		//TSPopt(&instance);
		//TSPadHoc(&instance);
		//TSPLocalBranching(&instance);
		TSPHardFixing(&instance);
		//TSPCallback(&instance);
		clock_t end = clock();
		if (instance.verbose >= 10)
		{
			printSolution(&instance, &instance.cplexSucc);
		}
		printf("\nCplex solution cost: %f\n", instance.cplexCost);
		printf("\nElapsed time: %.3f\n", (double)((end - instance.tstart)) / CLOCKS_PER_SEC);
	}
	// free allocated memory in "instance"
	free_memory(&instance);
}