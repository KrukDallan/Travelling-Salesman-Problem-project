#include "utilities.h"
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

void free_memory(Instance* instance)
{
	free(instance->xcoord);
	free(instance->ycoord);
	free(instance->dist);

	if (instance->best_sol != NULL)
		free(instance->best_sol);
	if (instance->best_succ != NULL)
		free(instance->best_succ);

	if (instance->greedySucc != NULL)
		free(instance->greedySucc);
	if (instance->emSucc != NULL)
		free(instance->emSucc);
	if (instance->tabuSucc != NULL)
		free(instance->tabuSucc);
	if (instance->vnsSucc != NULL)
		free(instance->vnsSucc);
	if (instance->geneticSucc != NULL)
		free(instance->geneticSucc);
	if (instance->cplexSucc != NULL)
		free(instance->cplexSucc);
}

int parse_command_line(int argc, char** argv, Instance* instance)
{
	int return_value = 0;
	// default values
	strcpy(instance->input_file, "NULL");
	strcpy(instance->variant, "");
	// not enough arguments
	if (argc < 2)
	{
		return 0;
	}

	for (int i = 1; i < argc; i++)
	{
		if (strcmp("-f", argv[i]) == 0)
		{
			strcpy(instance->input_file, argv[++i]);
			return_value += 1;
			continue;
		}
		if (strcmp("-verbose", argv[i]) == 0)
		{
			char* tmp = argv[++i];
			instance->verbose = atoi(tmp);
			printf("Verbose: %d\n", instance->verbose);
			tmp = NULL;
			continue;
		}
		if (strcmp("-h", argv[i]) == 0) // heuristic
		{
			strcpy(instance->heuristic, argv[++i]);
			if (strncmp(instance->heuristic, "genetic", 7) == 0)
			{
				instance->geneticFlag = 1;
			}
			else
			{
				instance->geneticFlag = 0;
			}
			continue;
		}
		if (strcmp("-v", argv[i]) == 0) // variant
		{
			strcpy(instance->variant, argv[++i]);
			if (strcmp(instance->variant, "2opt") == 0)
			{
				instance->twoOptFlag = 1;
				instance->graspFlag = 0;
			}
			if (strcmp(instance->variant, "grasp") == 0)
			{
				instance->twoOptFlag = 0;
				instance->graspFlag = 1;
			}
			continue;
		}
		if (strcmp("-c", argv[i]) == 0) // Cplex
		{
			strcpy(instance->heuristic, argv[++i]);
			strcpy(instance->variant, "");
			instance->geneticFlag = 0;
			instance->twoOptFlag = 0;
			instance->graspFlag = 0;
			continue;
		}
		if (strcmp("-t", argv[i]) == 0) // time limit
		{
			char* tmp = argv[++i];
			instance->timeLimit = atof(tmp); // da minuti a secondi
			instance->timeLimit *= 60.0;
			tmp = NULL;
			continue;
		}
		if (strcmp("-r", argv[i]) == 0)
		{
			char* tmp = argv[++i];
			instance->nnodes = atoi(tmp);
			strcat(strcpy(instance->solution_file, "RND_"), tmp);
			instance->dist = (double*)calloc(instance->nnodes * instance->nnodes, sizeof(double));
			tmp = NULL;
			continue;
		}
	}
	return return_value;
}

void read_input(Instance* instance)
{
	FILE* finput = fopen(instance->input_file, "r");
	FILE* fout = fopen("Output.txt", "w");
	// buffer for fgets()
	char line[200];
	// pointer for strtok()
	char* token;
	// flag to start storing coordinates
	int coordinates = 0;

	while (fgets(line, 200, finput) != NULL)
	{
		token = strtok(line, " :");

		if (strncmp(token, "NAME", 4) == 0)
		{
			char* name = strtok(NULL, " :");
			printf("NAME: %s", name);
			coordinates = 0;
			int len = strlen(name);
			if (len > 0 && name[len - 1] == '\n')
				name[len - 1] = '\0';
			strcpy(instance->solution_file, name);
			strcat(strcat(strcat(instance->solution_file, instance->heuristic), instance->variant), ".dat");
			printf("%s\n", instance->solution_file);
			continue;
		}

		if (strncmp(token, "DIMENSION", 9) == 0)
		{
			char* tmp = strtok(NULL, " :");
			instance->nnodes = atoi(tmp);
			int nn = instance->nnodes;
			instance->xcoord = (double*)calloc(nn, sizeof(double));
			instance->ycoord = (double*)calloc(nn, sizeof(double));
			instance->dist = (double*)calloc(nn * nn, sizeof(double));
			coordinates = 0;
		}

		if (strncmp(token, "NODE_COORD_SECTION", 18) == 0)
		{
			coordinates = 1;
			continue;
		}

		if (coordinates == 1)
		{
			char* x;
			char* y;
			char out[50];
			int i = atoi(token) - 1;
			if (i < 0 || i >= instance->nnodes)
			{
				coordinates = 0;
				continue;
			}
			x = strtok(NULL, " ");
			y = strtok(NULL, " ");
			instance->xcoord[i] = atof(x);
			instance->ycoord[i] = atof(y);
			snprintf(out, sizeof(out), "%s%s%s%s%s%s", token, ": ", x, " ", y, "\n");
			fputs(out, fout);
			continue;
		}
	}
	fclose(finput);
	fclose(fout);
}

void generate_random_points(Instance* instance, int size)
{
	// random x coordinates
	double* xpoints = calloc(size, sizeof(double));
	for (int i = 0; i < size; i++)
	{
		xpoints[i] = rand() % RAND_MAX;
	}
	// random y coordinates
	double* ypoints = calloc(size, sizeof(double));
	for (int i = 0; i < size; i++)
	{
		ypoints[i] = rand() % RAND_MAX;
	}

	instance->xcoord = (double*)calloc(size, sizeof(double));
	instance->ycoord = (double*)calloc(size, sizeof(double));
	for (int i = 0; i < size; i++)
	{
		instance->xcoord[i] = xpoints[i];
		instance->ycoord[i] = ypoints[i];
		// printf("Node %d at coordinates %f , %f\n", i, instance->xcoord[i], instance->ycoord[i]);
	}
	free(xpoints);
	free(ypoints);
}

void distance(int i, int j, Instance* instance)
{
	if (i == j)
	{
		instance->dist[i * instance->nnodes + j] = 0;
	}
	else
	{
		double xdist = instance->xcoord[i] - instance->xcoord[j];
		double ydist = instance->ycoord[i] - instance->ycoord[j];
		int distance = sqrt(xdist * xdist + ydist * ydist) + 0.49999999;
		instance->dist[i * instance->nnodes + j] = distance + 0.0;
	}
}

void printSolution(Instance* inst, int* succ)
{
	FILE* fout = fopen(inst->solution_file, "w");
	// array in which we write the coordinates; may be useless because we can just use best_sol to get the coordinates.
	double solution[2000];

	// create and fill the solution array
	int* sol = (int*)calloc(inst->nnodes, sizeof(int));
	sol[0] = 0;
	int tmp = succ[sol[0]];
	for (int i = 1; i < inst->nnodes; i++)
	{
		sol[i] = tmp;
		tmp = succ[tmp];
	}
	
	// this for and the two lines following it are used to write the minimum cost solution into the output file
	for (int i = 0; i < inst->nnodes; i++)
	{
		snprintf(solution, sizeof(solution), "%f%s%f%s", inst->xcoord[sol[i]], " ", inst->ycoord[sol[i]], "\n");
		fputs(solution, fout);
	}

	snprintf(solution, sizeof(solution), "%f%s%f%s", inst->xcoord[sol[0]], " ", inst->ycoord[sol[0]], "\n");
	fputs(solution, fout);
	fclose(fout);

	FILE* pipe = _popen("C:/Programmi/gnuplot/bin/gnuplot.exe -persist", "w");
	if (pipe != NULL)
	{
		fprintf(pipe, "plot '%s' with linespoints linestyle 6 \n", inst->solution_file);
		fflush(pipe);
	}
	else
		puts("Could not open the file\n");
	_pclose(pipe);
	free(sol);
}

int two_opt(Instance* instance, int* successors, double* cost, int* tabuIter, int fiveOptFlag)
{
	double delta = MAX_VAL;
	int node_a = 0;
	int node_b = 0;
	int flag = 0;
	// tabu-search
	if (strcmp(instance->heuristic, "tabu") == 0)
	{

		while (1)
		{
			flag = 0;
			delta = MAX_VAL;
			node_a = 0;
			node_b = 0;
			for (int i = 0; i < instance->nnodes; i++)
			{
				for (int j = 0; j < instance->nnodes; j++)
				{
					if (j != i)
					{
						double cost_a = instance->dist[i * instance->nnodes + successors[i]];
						double cost_b = instance->dist[j * instance->nnodes + successors[j]];
						double old_cost = cost_a + cost_b;

						double ncost_a = instance->dist[i * instance->nnodes + j];
						double ncost_b = instance->dist[successors[i] * instance->nnodes + successors[j]];
						double new_cost = ncost_a + ncost_b;

						double tmp = new_cost - old_cost;
						if (tmp < 0 && tmp < delta)
						{
							delta = tmp;
							node_a = i;
							node_b = j;
							flag = 1;
						}
					}
				}
			}
			int succ_a_old = successors[node_a];
			int succ_b_old = successors[node_b];
			if (flag != 0)
			{
				*cost += delta;

				// store in tmp the successor of node A
				int tmp = successors[node_a];
				// size of the stack in which the successors will be saved
				int size = 0;
				int* stack = calloc(instance->nnodes, sizeof(int));
				for (int i = 0; i < instance->nnodes; i++)
				{
					if (tmp == node_b)
					{
						break;
					}
					else
					{
						stack[i] = tmp;
						tmp = successors[tmp];
						size += 1;
					}
				}
				while (size > 0)
				{
					successors[tmp] = stack[size - 1];
					tmp = stack[size - 1];
					size -= 1;
				}
				successors[node_a] = node_b;
				successors[succ_a_old] = succ_b_old;
				free(stack);
			}
			else
			{
				if (* cost < instance->tabuCost)
				{
					for (int i = 0; i < instance->nnodes; i++)
					{
						instance->tabuSucc[i] = successors[i];
					}
					instance->tabuCost = *cost;
				}
				tabuIter[node_b] = instance->numIter;
				return -2;
			}
		}
	}

	// vns
	else if (strcmp(instance->heuristic, "vns") == 0)
	{
		if (fiveOptFlag == 0)
		{
			while (1)
			{
				flag = 0;
				double delta = MAX_VAL;
				int node_a = 0;
				int node_b = 0;
				for (int i = 0; i < instance->nnodes; i++)
				{
					for (int j = 0; j < instance->nnodes; j++)
					{
						if (j != i)
						{
							double cost_a = instance->dist[i * instance->nnodes + successors[i]];
							double cost_b = instance->dist[j * instance->nnodes + successors[j]];
							double old_cost = cost_a + cost_b;

							double ncost_a = instance->dist[i * instance->nnodes + j];
							double ncost_b = instance->dist[successors[i] * instance->nnodes + successors[j]];
							double new_cost = ncost_a + ncost_b;

							double tmp = new_cost - old_cost;
							if (tmp < 0 && tmp < delta)
							{
								delta = tmp;
								node_a = i;
								node_b = j;
								flag = 1;
							}
						}
					}
				}
				int succ_a_old = successors[node_a];
				int succ_b_old = successors[node_b];
				if (flag != 0)
				{
					*cost += delta;

					// store in tmp the successor of node A
					int tmp = successors[node_a];
					// size of the stack in which the successors will be saved
					int size = 0;
					int* stack = calloc(instance->nnodes, sizeof(int));
					for (int i = 0; i < instance->nnodes; i++)
					{
						if (tmp == node_b)
						{
							break;
						}
						else
						{
							stack[i] = tmp;
							tmp = successors[tmp];
							size += 1;
						}
					}
					while (size > 0)
					{
						successors[tmp] = stack[size - 1];
						tmp = stack[size - 1];
						size -= 1;
					}
					successors[node_a] = node_b;
					successors[succ_a_old] = succ_b_old;
					free(stack);
				}
				else
				{
					if (*cost < instance->vnsCost)
					{
						for (int i = 0; i < instance->nnodes; i++)
						{
							instance->vnsSucc[i] = successors[i];
						}
						instance->vnsCost = *cost;
					}
					return 1;
				}
			}
		}

		else
		{
			int size = 10;
			int* pool = (int*)calloc(size, sizeof(int));
			int* isUsed = (int*)calloc(instance->nnodes, sizeof(int));

			for (int i = 0; i < instance->nnodes; i++)
			{
				isUsed[i] = -1;
			}

			for (int i = 0; i < size; i++)
			{
				while (1)
				{
					int rv = rand() % instance->nnodes;
					if (isUsed[rv] == -1)
					{
						pool[i] = rv;
						isUsed[rv] = 1;
						break;
					}
				}
			}
			free(isUsed);

			for (int i = 0; i < size - 1; i += 2)
			{
				node_a = pool[i];
				node_b = pool[i + 1];

				double cost_a = instance->dist[node_a * instance->nnodes + successors[node_a]];
				double cost_b = instance->dist[node_b * instance->nnodes + successors[node_b]];
				double old_cost = cost_a + cost_b;

				double ncost_a = instance->dist[node_a * instance->nnodes + node_b];
				double ncost_b = instance->dist[successors[node_a] * instance->nnodes + successors[node_b]];
				double new_cost = ncost_a + ncost_b;

				double diff = new_cost - old_cost;
				// instance->current_value += diff;
				*cost += diff;
				int succ_a_old = successors[node_a];
				int succ_b_old = successors[node_b];
				int tmp = successors[node_a];
				int ssize = 0;
				int* stack = calloc(instance->nnodes, sizeof(int));
				// this loop will not go through every node in every case, it will just fill the stack array until we reach node_b
				for (int i = 0; i < instance->nnodes; i++)
				{
					if (tmp == node_b)
					{
						break;
					}
					else
					{
						stack[i] = tmp;
						tmp = successors[tmp];
						ssize += 1;
					}
				}
				// at this point we reverse the order of successors
				while (ssize > 0)
				{
					successors[tmp] = stack[ssize - 1];
					tmp = stack[ssize - 1];
					ssize -= 1;
				}
				successors[node_a] = node_b;
				successors[succ_a_old] = succ_b_old;
				free(stack);
			}
			//printf("Cost5 opt: %f\n", *cost);
			free(pool);
			for (int i = 0; i < instance->nnodes; i++)
			{
				// instance->best_succ[i] = instance->vnsSucc[i];
			}
			// printSolution(instance);
			return -1; // we return -1 so in the next call 2-opt will be applied
		}
	}

	// greedy or genetic or cplex
	else
	{
		for (int i = 0; i < instance->nnodes; i++)
		{
			for (int j = 0; j < instance->nnodes; j++)
			{
				if (j != i)
				{
					double cost_a = instance->dist[i * instance->nnodes + successors[i]];
					double cost_b = instance->dist[j * instance->nnodes + successors[j]];
					double old_cost = cost_a + cost_b;

					double ncost_a = instance->dist[i * instance->nnodes + j];
					double ncost_b = instance->dist[successors[i] * instance->nnodes + successors[j]];
					double new_cost = ncost_a + ncost_b;

					double tmp = new_cost - old_cost;
					if (tmp < 0 && tmp < delta)
					{
						delta = tmp;
						node_a = i;
						node_b = j;
						flag = 1;
					}
				}
			}
		}
		int succ_a_old = successors[node_a];
		int succ_b_old = successors[node_b];
		// if a shorter connection has been found:
		if (flag != 0)
		{
			// update cost
			*cost += delta;
			// store in tmp the successor of node A
			int tmp = successors[node_a];
			// size of the stack in which the successors will be saved
			int size = 0;
			int* stack = calloc(instance->nnodes, sizeof(int));
			for (int i = 0; i < instance->nnodes; i++)
			{
				if (tmp == node_b)
				{
					break;
				}
				else
				{
					stack[i] = tmp;
					tmp = successors[tmp];
					size += 1;
				}
			}
			// update connection
			while (size > 0)
			{
				successors[tmp] = stack[size - 1];
				tmp = stack[size - 1];
				size -= 1;
			}
			successors[node_a] = node_b;
			successors[succ_a_old] = succ_b_old;
			free(stack);
			return flag;
		}
		// else flag is 0, hence the 2 opt method is stopped
		else
		{
			return flag;
		}
	}
}

void succ2sol(Instance* instance, int* succ, int* sol)
{
	sol[0] = succ[0];
	int tmp = succ[sol[0]];
	for (int i = 1; i < instance->nnodes; i++)
	{
		sol[i] = tmp;
		tmp = succ[sol[i]];
	}
}
