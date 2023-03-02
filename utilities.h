#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsp.h"
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#define MAX_VAL DBL_MAX;

int parse_command_line(int argc, char** argv, Instance* instance);
void read_input(Instance* instance);
void generate_random_points(Instance* instance, int size);
void free_memory(Instance* instance);
void distance(int i, int j, Instance* instance);
void printSolution(Instance* inst, int* succ);

void succ2sol(Instance* instance, int* succ, int* sol);

int two_opt(Instance* instance, int* successors, double* cost, int* tabuIter, int fiveOptFlag);

#endif