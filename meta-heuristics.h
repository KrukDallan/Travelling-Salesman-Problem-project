#ifndef META_H
#define META_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#include "tsp.h"
#include "utilities.h"
#include "constructing_heuristics.h"

#define MAX_VAL DBL_MAX;

//Tabù search
int tabu_opt(Instance* instance, int* succ, double* cost, int* tabuIter);
void tabuSearch(Instance* instance);

//Variable Neighborhood Search
void VNS(Instance* instance);

//Genetic algorithm
void genetic(Instance* instance);
void generatePopulation(Population* population, int nnodes, int pSize, Instance* instance);
void crossover(Instance* instance, Population* population, int parent1, int parent2, int* chromosome, int* childSize, int startingNode, double* cost);
void repair(Population* child, Instance* instance, int cSize, int nnodes);
void calc_fitness(Instance* instance, Population* population, int childIndex);
void kill_weak(Instance* instance, Population* population, const int pSize, Population* child, const int cSize);
void gen_stats(Population* population, int pSize, double* best_cost, double* average, int* champion_idx);

#endif
