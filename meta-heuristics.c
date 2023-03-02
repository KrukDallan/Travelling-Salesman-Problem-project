#include "meta-heuristics.h"
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )


int tabu_opt(Instance* instance, int* succ, double* cost, int* tabuIter)
{
	double delta = MAX_VAL;
	int node_a = 0;
	int node_b = 0;
	int flag = 0;
	int size = 0;
	if (instance->nnodes <= 100)
	{
		size = (int)ceil((instance->nnodes) / 2);
	}
	else
	{
		size = (int)ceil(sqrt(instance->nnodes));
	}
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

	for (int i = 0; i < size; i++)
	{
		if (tabuIter[i] != INT_MIN)
		{
			continue;
		}
		for (int j = 0; j < size; j++)
		{
			if (j == i || tabuIter[j] != INT_MIN)
			{
				continue;
			}
			else
			{
				double cost_a = instance->dist[pool[i] * instance->nnodes + succ[pool[i]]];
				double cost_b = instance->dist[pool[j] * instance->nnodes + succ[pool[j]]];
				double old_cost = cost_a + cost_b;

				double ncost_a = instance->dist[pool[i] * instance->nnodes + pool[j]];
				double ncost_b = instance->dist[succ[pool[i]] * instance->nnodes + succ[pool[j]]];
				double new_cost = ncost_a + ncost_b;

				double tmp = new_cost - old_cost;
				//printf("newc_cost - old_cost = %f\n", tmp);
				if (tmp < delta && tmp > 0)
				{
					delta = tmp;
					//printf("delta is: %f\n", delta);
					node_a = pool[i];
					node_b = pool[j];
					flag = 0;
				}
				if (tmp < 0)
				{
					delta = tmp;
					//printf("delta is: %f\n", delta);
					node_a = pool[i];
					node_b = pool[j];
					flag = -1;
				}
			}
		}
	}
	//printf("flag in tabu is %d and node_b is %d\n", flag, node_b);
	int succ_a_old = succ[node_a];
	int succ_b_old = succ[node_b];

	int tmp = succ[node_a];
	int ssize = 0;
	int* stack = calloc(instance->nnodes, sizeof(int));
	//this loop will not go through every node in every case, it will just fill the stack array until we reach node_b
	for (int i = 0; i < instance->nnodes; i++)
	{
		if (tmp == node_b)
		{
			break;
		}
		else
		{
			stack[i] = tmp;
			tmp = succ[tmp];
			ssize += 1;
		}
	}
	//at this point we reverse the order of successors
	while (ssize > 0)
	{
		succ[tmp] = stack[ssize - 1];
		tmp = stack[ssize - 1];
		ssize -= 1;
	}
	succ[node_a] = node_b;
	succ[succ_a_old] = succ_b_old;

	*cost += delta;

	/*for (int i = 0; i < instance->nnodes; i++)
	{
		instance->best_succ[i] = instance->tabuSucc[i];
	}*/
	//instance->best_cost = instance->current_value;
	//printf("Tabu cost: %f\n", instance->current_cost);

	free(stack);
	free(pool);
	if (flag == -1)
	{
		return -1;
	}
	return node_b;
}

void tabuSearch(Instance* instance)
{
	//We start by looking for a generic solution, which will be the starting solution
	//We can use the a greedy heuristic since it is pretty fast
	//initialization of some arrays

	int* uncovered_nodes = (int*)calloc(instance->nnodes, sizeof(int));
	//will contain the current solution
	int* succ = (int*)calloc(instance->nnodes, sizeof(int));

	//will hold the index returned by min_dist_node
	int next_node_index = 0;

	for (int i = 0; i < instance->nnodes; i++)
	{
		//filling the uncovered_nodes array
		uncovered_nodes[i] = i;
		succ[i] = -1;
	}

	int starting_node = rand() % instance->nnodes;
	if (instance->verbose > 50)
	{
		printf("\nTabu first node: %d\n", starting_node);
	}
	int sn_backup = starting_node;
	//cost value of the j-th cycle
	double current_cost = 0;

	//setting the first element of the current solution as the starting node
	succ[starting_node] = starting_node;

	int size = instance->nnodes - 1;
	//if j==size then j==47, so uncovered_nodes[j]=uncovered_nodes[size] is useless
	uncovered_nodes[starting_node] = -1;

	//used to fill node_list array
	int second_flag = 0;

	while (size > 0)
	{
		next_node_index = min_dist_node(starting_node, uncovered_nodes, size, instance, &current_cost, 0, second_flag);
		succ[starting_node] = uncovered_nodes[next_node_index];
		
		//re-assigning starting node so in the next call of min_dist_node the starting_node value passed will be the current "next_node"
		starting_node = uncovered_nodes[next_node_index];
		if (instance->verbose > 50)
		{
			printf("Next node: %d\n", starting_node);
		}
		uncovered_nodes[next_node_index] = -1;

		size = size - 1;
		if (second_flag == 9)
		{
			second_flag = 0;
		}
		else
		{
			second_flag += 1;
		}

	}
	for (int i = 0; i < instance->nnodes; i++)
	{
		if (succ[i] == -1)
		{
			succ[i] = sn_backup;
			break;
		}
	}
	if (current_cost < instance->tabuCost)
	{
		//printf("Cost at iteration %d is: %f\n", j, current_cost);
		instance->tabuCost = current_cost;
		if (instance->verbose >= 10)
		{
			printf("New best cost: %f\n", instance->tabuCost);
		}
		instance->starting_node = sn_backup;
		for (int i = 0; i < instance->nnodes; i++)
		{
			instance->tabuSucc[i] = succ[i];
		}
	}
	free(uncovered_nodes);
	//end of greedy search: we now have our starting solution

	//we can now apply 2-opt on a subset of the neighbors of the starting solution,
	//that is, we choose some random nodes on which to apply 2-opt and then select the minimum solution among those
	//It seems reasonable to choose as the number of nodes the square root of the total number of nodes, so that the time consumed is not too large,
	//even if the number of neighbors is not that big
	int numIter = 0;
	//instance->numIter = numIter;
	int* tabuIter = (int*)calloc(instance->nnodes, sizeof(int));
	int tenure = 0;
	int maxTenure = 0;
	int minTenure = 0;
	int increase = 1;
	if (instance->nnodes <= 50)
	{
		tenure = (int)ceil(instance->nnodes / 10);
		maxTenure = (int)ceil(tenure + tenure / 2);
		minTenure = tenure;
	}
	else
	{
		tenure = MIN((int)ceil(instance->nnodes / 10), 100);
		maxTenure = (int)ceil(tenure + tenure / 2);
		minTenure = tenure;
	}

	for (int i = 0; i < instance->nnodes; i++)
	{
		tabuIter[i] = INT_MIN;
	}
	// initial call: 2-opt
	double cost = instance->tabuCost;
	//int flag = two_opt(instance, succ, &cost, tabuIter, 0); // flag = -2
	int flag = 0;
	int tabuNode = 0;
	while (1)
	{
		//update tabu iter
		for (int i = 0; i < instance->nnodes; i++)
		{
			if (numIter - tabuIter[i] >= tenure)
			{
				tabuIter[i] = INT_MIN;
			}
		}
		if (flag == 0)
		{
			flag = two_opt(instance, succ, &cost, tabuIter, 0);
		}

		if (flag == -2)
		{
			
			tabuNode = tabu_opt(instance, succ, &cost, tabuIter);
			
			if (tabuNode != -1)
			{
				tabuIter[tabuNode] = numIter;
			}
			else
			{
				flag = 0;
			}
		}
		if (numIter % 50 == 0)
		{
			flag = 0;
		}

		//dynamic tenure
		tenure = (increase == 1) ? tenure + 1 : tenure - 1;
		if (tenure == maxTenure)
		{
			increase = 0;
		}
		else
		{
			if (tenure == minTenure)
			{
				increase = 1;
			}
		}
		numIter += 1;
		//instance->numIter = numIter;
		clock_t end = clock();
		double diff = (double)((end - instance->tstart)) / CLOCKS_PER_SEC;
		if (diff > instance->timeLimit)
		{
			break;
		}

	}
	if (instance->verbose >= 10)
	{
		printSolution(instance, instance->tabuSucc);
	}


	free(succ);
	free(tabuIter);
}


void VNS(Instance* instance)
{
	//We start by looking for a generic solution, which will be the starting solution
	//We can use the greedy heuristic since it is pretty fast
	//initialization of best_sol and of uncovered nodes

	//will contain the current solution
	int* succ = (int*)calloc(instance->nnodes, sizeof(int));

	//will hold the index returned by min_dist_node
	int next_node_index = 0;

	int* uncovered_nodes = (int*)calloc(instance->nnodes, sizeof(int));
	for (int i = 0; i < instance->nnodes; i++)
	{
		//filling the uncovered_nodes array
		uncovered_nodes[i] = i;
		succ[i] = -1;
	}

	int starting_node = rand() % instance->nnodes;
	if (instance->verbose > 50)
	{
		printf("\nVNS first node: %d\n", starting_node);
	}
	int sn_backup = starting_node;
	//cost value of the j-th cycle
	double current_cost = 0;

	//setting the first element of the current solution as the starting node
	succ[starting_node] = starting_node;

	int size = instance->nnodes - 1;
	//if j==size then j==47, so uncovered_nodes[j]=uncovered_nodes[size] is useless
	uncovered_nodes[starting_node] = -1;

	//used to fill node_list array
	int second_flag = 0;

	while (size > 0)
	{
		next_node_index = min_dist_node(starting_node, uncovered_nodes, size, instance, &current_cost, 0, second_flag);
		succ[starting_node] = uncovered_nodes[next_node_index];

		//re-assigning starting node so in the next call of min_dist_node the starting_node value passed will be the current "next_node"
		starting_node = uncovered_nodes[next_node_index];
		if (instance->verbose > 50)
		{
			printf("\nNext node: %d\n", starting_node);
		}
		uncovered_nodes[next_node_index] = -1;

		size = size - 1;
		if (second_flag == 9)
		{
			second_flag = 0;
		}
		else
		{
			second_flag += 1;
		}

	}
	for (int i = 0; i < instance->nnodes; i++)
	{
		if (succ[i] == -1)
		{
			succ[i] = sn_backup;
			break;
		}
	}
	if (current_cost < instance->vnsCost)
	{
		//printf("Cost at iteration %d is: %f\n", j, current_cost);
		instance->vnsCost = current_cost;
		if (instance->verbose >= 10)
		{
			printf("New best cost: %f\n", instance->vnsCost);
		}
		for (int i = 0; i < instance->nnodes; i++)
		{
			instance->vnsSucc[i] = succ[i];
		}
	}
	free(uncovered_nodes);
	//end of greedy search: we have our starting solution

	//we can now apply 2-opt on a subset of the neighbors of the starting solution,
	//that is, we choose some random nodes on which apply 2-opt and then select the minimum solution among those
	//It seems reasonable to choose as number nodes the square root of the total number of nodes, so that the time consumed is not so large,
	//even if the number of neighbors is not so big
	int numIter = 0;
	int fiveOptFlag = 0;

	//printSolution(instance);
	while (1)
	{
		int flag = two_opt(instance, succ, &current_cost, NULL, fiveOptFlag);
		if (flag == 1)
		{
			fiveOptFlag = 1; //apply 5-opt
		}
		else if (flag == -1)
		{
			fiveOptFlag = 0; //apply 2-opt
		}
		else
		{
			break;
		}
		numIter += 1;
		if (numIter == 100)
		{
			two_opt(instance, succ, &current_cost, NULL, 0);
		}
		clock_t end = clock();
		double diff = (double)((end - instance->tstart)) / CLOCKS_PER_SEC;
		if (diff > instance->timeLimit)
		{
			break;
		}

	}
	free(succ);
	if (instance->verbose >= 10)
	{
		printSolution(instance, instance->vnsSucc);
	}
}


void generatePopulation(Population* population, int nnodes, int pSize, Instance* instance)
{
	for (int i = 0; i < nnodes; i++)
	{
		// Initializations
		population[i].chromosome = (int*)calloc(nnodes, sizeof(int));
		int starting_node = i;
		population[i].startingNode = starting_node;
		population[i].cost = 0.0;

		greed_search(instance, 0, 0, population[i].chromosome, population[i].startingNode, &population[i].cost);
	}

	int index = 0;
	if (nnodes < pSize) // We now use greedy again but this time with GRASP
	{
		index = nnodes;
		for (int i = index; i < pSize; i++)
		{
			// Initializations
			population[i].chromosome = (int*)calloc(nnodes, sizeof(int));
			int starting_node = rand() % nnodes;
			population[i].startingNode = starting_node;
			population[i].cost = 0.0;

			greed_search(instance, 1, 0, population[i].chromosome, population[i].startingNode, &population[i].cost);
		}
	}
}

void crossover(Instance* instance, Population* population, int parent1, int parent2, int* chromosome, int* childSize, int startingNode, double* cost)
{
	for (int i = 0; i < instance->nnodes; i++)
	{
		chromosome[i] = -1;
	}

	int* covered = calloc(instance->nnodes, sizeof(int));
	int ind = (int)instance->nnodes / 2;
	int succ_node = population[parent1].chromosome[startingNode];
	//in chromosome[sn] we save the succ of sn
	int sn = startingNode;

	for (int i = 0; i < ind; i++)
	{
		chromosome[sn] = succ_node;
		*cost = *cost + instance->dist[sn * instance->nnodes + succ_node];
		covered[succ_node] = 1;
		sn = succ_node;
		*childSize += 1;
		succ_node = population[parent1].chromosome[sn];
	}
	succ_node = sn;
	for (int i = 0; i < instance->nnodes; i++)
	{
		//check the successor of sn wrt parent2
		int tmp = population[parent2].chromosome[succ_node];

		if (covered[tmp] == 1)
		{
			succ_node = population[parent2].chromosome[tmp];
		}
		else
		{
			chromosome[sn] = tmp;
			*cost = *cost + instance->dist[sn * instance->nnodes + tmp];
			*childSize += 1;
			covered[tmp] = 1;
			sn = tmp;
			succ_node = population[parent2].chromosome[tmp];
			if (tmp == startingNode)
			{
				break;
			}
		}
	}
	if (chromosome[sn] == -1)
	{
		chromosome[sn] = startingNode;
		*childSize += 1;
		*cost = *cost + instance->dist[sn * instance->nnodes + startingNode];
	}
	free(covered);
}

void calc_fitness(Instance* instance, Population* population, int childIndex)
{
	int pred = population[childIndex].chromosome[0];
	population[childIndex].cost = 0;
	for (int i = 0; i < instance->nnodes; i++)
	{
		int node = population[childIndex].chromosome[i];
		population[childIndex].cost += instance->dist[pred * instance->nnodes + node];
	}
	population[childIndex].cost += instance->dist[pred * instance->nnodes + population[childIndex].chromosome[0]];
}

void repair(Population* child, Instance* instance, int cSize, int nnodes)
{
	for (int i = 0; i < cSize; i++)
	{
		if (child[i].childSize < nnodes)
		{
			extra_mileage(instance, child[i].chromosome, &child[i].childSize, &child[i].cost);

		}
		int tmp = 1;
		//tmp = two_opt(instance, child[i].chromosome, &child[i].cost, NULL, 0);
		while (tmp == 1)
		{
			tmp = two_opt(instance, child[i].chromosome, &child[i].cost, NULL, 0);
		}
		//printf("Child cost: %f\n", child[i].cost);
	}
}

int cmpfunc(const void* a, const void* b)
{
	const Population* right = a;
	const Population* left = b;
	return (right->cost - left->cost);
}

void kill_weak(Instance* instance, Population* population, const int pSize, Population* child, const int cSize)
{
	int tSize = pSize + cSize;
	Population* total = calloc(tSize, sizeof(Population));
	int count = 0;

	//saving both population and generated childs of the current generation into a temp collection
	for (int i = 0; i < pSize; i++) 
	{
		total[count++] = population[i];
	}
	for (int i = 0; i < cSize; i++) 
	{
		total[count++] = child[i];
	}

	count = 0;

	qsort(total, tSize, sizeof(Population), cmpfunc);

	while (count < pSize)
	{
		memcpy(population[count].chromosome, total[count].chromosome, sizeof(int) * instance->nnodes);
		population[count].cost = total[count].cost;
		count++;
	}
	for (int i = 0; i < cSize; i++)
	{
		free(child[i].chromosome);
		child[i].cost = 0.0;
	}
	free(total);
}

void gen_stats(Population* population, int pSize, double* best_cost, double* average, int* champion_idx) 
{
	*best_cost = MAX_VAL;
	*average = 0;
	for (int i = 0; i < pSize; i++) 
	{
		double fitness = population[i].cost;
		*average += fitness;
		if (fitness <= *best_cost) 
		{
			*best_cost = fitness;
			*champion_idx = i;
		}
	}
	*average /= pSize;
}

void genetic(Instance* instance)
{
	clock_t start = clock();
	// Initializations
	int nnodes = instance->nnodes;
	int pSize = 1000;
	int cSize = (int)pSize / 2;
	double best_cost = MAX_VAL;


	Population* population = (Population*)calloc(pSize, sizeof(Population));
	Population* child = (Population*)calloc(cSize, sizeof(Population));

	// Generating the population
	generatePopulation(population, nnodes, pSize, instance); 

	int numIter = 0;
	while (1)
	{
		//printf("Starting iteration %d\n", numIter + 1);
		for (int i = 0; i < cSize; i++)
		{
			child[i].chromosome = (int*)calloc(nnodes, sizeof(int));
			child[i].childSize = 0;
			int parent1 = rand() % pSize;
			int parent2 = rand() % pSize; //rand() % pSize -1
			if (parent1 == parent2)
			{
				while (1)
				{
					parent2 = rand() % pSize;
					if (parent1 != parent2)
					{
						break;
					}
				}
			}
			child[i].startingNode = population[parent1].startingNode;
			crossover(instance, population, parent1, parent2, child[i].chromosome, &child[i].childSize, child[i].startingNode, &child[i].cost);
		}

		repair(child, instance, cSize, instance->nnodes);

		kill_weak(instance, population, pSize, child, cSize);
		if (instance->verbose >= 10)
		{
			printf("New best cost: %f\n", population[0].cost);
		}
		numIter += 1;

		clock_t end = clock();
		double diff = (double)(end - start) / CLOCKS_PER_SEC;
		//printf("Elapsed time: %f\n", diff);
		if (diff > instance->timeLimit)
		{
			break;
		}
	}
	
	instance->geneticCost = population[0].cost;
	if (instance->verbose >= 10)
	{
		printSolution(instance, population[0].chromosome);
	}
	for (int i = 0; i < pSize; i++)
	{
		free(population[i].chromosome);
	}
	free(population);
	free(child);
}