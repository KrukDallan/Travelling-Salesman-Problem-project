#ifndef CONSTRUCTIVE_H
#define CONSTRUCTIVE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#include "tsp.h"
#include "utilities.h"

#define MAX_VAL DBL_MAX;

//Greedy heuristic
void greed_search(Instance* instance, int grasp_flag, int two_opt_flag, int* successors, int startingNode, double* cost);
int min_dist_node(int node, int uncovered_nodes[], int size, Instance* instance, double* current_cost, int grasp_flag, int second_choice_flag);

//extra-mileage
void extra_mileage(Instance* instance, int* successors, int* size, double* cost);
emstruct min_mileage(Instance* instance, int succ[], int uncovered_nodes[], double* current_cost);
int max_dist_node(int node, int uncovered_nodes[], Instance* instance, double* current_cost);


#endif
