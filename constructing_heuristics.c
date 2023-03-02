#include "constructing_heuristics.h"
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

int min_dist_node(int node, int uncovered_nodes[], int size, Instance* instance, double* current_cost, int grasp_flag, int second_choice_flag)
{
	int nn = instance->nnodes;
	int min_node_index = -1;
	// variables used for GRASP
	double min_dist = MAX_VAL;
	int second_min_node_index = 0;
	double second_min_dist = 0;
	for (int i = 0; i < nn; i++)
	{
		// skip already covered nodes
		if (uncovered_nodes[i] != -1)
		{
			// skip case: node == uncovered_nodes[i]
			if (instance->dist[node * nn + uncovered_nodes[i]] != 0)
			{
				if (instance->dist[node * nn + uncovered_nodes[i]] < min_dist)
				{
					second_min_dist = min_dist;
					second_min_node_index = min_node_index;
					min_dist = instance->dist[node * nn + uncovered_nodes[i]];
					min_node_index = i;
				}
			}
		}
	}
	// if using GRASP:
	if (grasp_flag == 1 && second_choice_flag == 9)
	{
		// case in which there exists a second best choice
		if (second_min_node_index != -1)
		{
			*current_cost = *current_cost + second_min_dist;
			return second_min_node_index;
		}
		else
		{
			*current_cost = *current_cost + min_dist;
			return min_node_index;
		}
	}
	// else:
	else
	{
		*current_cost = *current_cost + min_dist;
		return min_node_index;
	}
}

void greed_search(Instance* instance, int grasp_flag, int two_opt_flag, int* successors, int startingNode, double* cost)
{
	if (instance->geneticFlag == 0)
	{
		// will hold the index returned by min_dist_node
		int next_node_index = 0;

		for (int j = 0; j < instance->nnodes; j++)
		{
			int* uncovered_nodes = (int*)calloc(instance->nnodes, sizeof(int));
			for (int i = 0; i < instance->nnodes; i++)
			{
				// filling the uncovered_nodes array
				uncovered_nodes[i] = i;
				// succ[i] = -1;
				instance->greedySucc[i] = -1;
			}

			int starting_node = j;
			if (instance->verbose > 50)
			{
				printf("\nGREEDY first node: %d\n", starting_node);
			}
			// cost value of the j-th cycle
			double current_cost = 0.0;

			// setting the first element of the current solution as the starting node
			// succ[starting_node] = starting_node;
			instance->greedySucc[starting_node] = starting_node;

			int size = instance->nnodes - 1;
			// if j==size then j==47, so uncovered_nodes[j]=uncovered_nodes[size] is useless
			uncovered_nodes[starting_node] = -1;

			// used for selecting second best solution
			int second_choice_flag = 0;

			while (size > 0)
			{
				next_node_index = min_dist_node(starting_node, uncovered_nodes, size, instance, &current_cost, grasp_flag, second_choice_flag);
				instance->greedySucc[starting_node] = uncovered_nodes[next_node_index];
				if (instance->verbose > 50) 
				{ 
					printf("GREEDY next node: %d\n", starting_node); 
				}
				// re-assigning starting node so in the next call of min_dist_node the starting_node value passed will be the current "next_node"
				starting_node = uncovered_nodes[next_node_index];
				uncovered_nodes[next_node_index] = -1;

				size = size - 1;
				if (second_choice_flag == 9)
				{
					second_choice_flag = 0;
				}
				else
				{
					second_choice_flag += 1;
				}
			}
			for (int i = 0; i < instance->nnodes; i++)
			{
				// fill the missing successor of the initial node
				if (instance->greedySucc[i] == -1)
				{
					instance->greedySucc[i] = j;
					current_cost += instance->dist[i * instance->nnodes + j];
					break;
				}
			}
			if (current_cost < instance->greedyCost)
			{
				instance->greedyCost = current_cost;
				if (instance->verbose >= 10) 
				{ 
					printf("New best cost: %f\n", instance->greedyCost); 
				}
				if (cost != NULL)
				{
					*(cost) = current_cost;
				}
				for (int i = 0; i < instance->nnodes; i++)
				{
					if (successors != NULL)
					{
						successors[i] = instance->greedySucc[i];
					}
				}
			}
			free(uncovered_nodes);
		}
		int size = 0;
		int flag = two_opt_flag;
		// apply 2-opt
		while (flag == 1)
		{
			flag = two_opt(instance, instance->greedySucc, &instance->greedyCost, NULL, 0);
		}
		if (cost != NULL)
		{
			*(cost) = instance->greedyCost;
		}
		if (instance->verbose >= 10)
		{
			printSolution(instance, instance->greedySucc);
		}
	}
	// case in which we want to apply the greedy algo once, i.e. to create the population of the genetic algo
	else
	{
		int nnodes = instance->nnodes;
		int starting_node = startingNode;
		int* uncoveredNodes = (int*)calloc(nnodes, sizeof(int));
		int nextNodeIndex = 0;

		for (int j = 0; j < nnodes; j++)
		{
			// Filling the uncovered_nodes array
			uncoveredNodes[j] = j;
			successors[j] = -1;
		}

		// Setting the first element of the current solution as the starting node
		successors[starting_node] = starting_node;
		// Size of the solution
		int size = nnodes - 1;
		// If j==size then j==nnodes, so uncovered_nodes[j]=uncovered_nodes[size] is useless
		uncoveredNodes[starting_node] = -1;
		int second_choice_flag = 0;
		while (size > 0)
		{
			nextNodeIndex = min_dist_node(starting_node, uncoveredNodes, size, instance, cost, grasp_flag, second_choice_flag);
			successors[starting_node] = uncoveredNodes[nextNodeIndex];

			// re-assigning starting node so in the next call of min_dist_node the starting_node value passed will be the current "next_node"
			starting_node = uncoveredNodes[nextNodeIndex];
			uncoveredNodes[nextNodeIndex] = -1;
			size = size - 1;

			if (second_choice_flag == 9)
			{
				second_choice_flag = 0;
			}
			else
			{
				second_choice_flag += 1;
			}
		}
		for (int k = 0; k < nnodes; k++)
		{
			if (successors[k] == -1)
			{
				successors[k] = startingNode;
				*cost += instance->dist[k * instance->nnodes + startingNode];
				break;
			}
		}
		free(uncoveredNodes);
	}
}

int max_dist_node(int node, int uncovered_nodes[], Instance* instance, double* current_cost)
{
	int nn = instance->nnodes;
	int max_node_index = 0;
	double max_dist = 0;
	for (int i = 0; i < nn; i++)
	{
		if (uncovered_nodes[i] != -1)
		{
			if (instance->dist[node * nn + uncovered_nodes[i]] == 0)
			{
				continue;
			}
			else
			{
				if (instance->dist[node * nn + uncovered_nodes[i]] > max_dist)
				{
					max_dist = instance->dist[node * nn + uncovered_nodes[i]];
					max_node_index = i;
				}
			}
		}
	}
	*current_cost = *current_cost + max_dist;
	return max_node_index;
}

emstruct min_mileage(Instance* instance, int succ[], int uncovered_nodes[], double* current_cost)
{
	// return index of the node in uncovered_nodes[] for which the extra mileage is minimum
	// in a generic iteration, compute extra mileage for all covered nodes and all uncovered nodes and pick the minimum e-m:
	// E-M(i,h) = c[i,h] + c[h, succ[i]] - c[i, succ[i]]
	double absolute_mileage = MAX_VAL;
	double relative_mileage = 0;

	int nn = instance->nnodes;
	double min_dist = MAX_VAL;
	int min_node_index = 0;
	int min_node_succ = 0;
	// for each node i in succ, find the node in uncovered nodes at minimum distance from i
	for (int i = 0; i < nn; i++)
	{
		if (succ[i] == -1)
		{
			continue;
		}
		else
		{
			for (int h = 0; h < nn; h++)
			{
				if (uncovered_nodes[h] != -1)
				{
					relative_mileage = instance->dist[i * nn + uncovered_nodes[h]] + instance->dist[uncovered_nodes[h] * nn + succ[i]] - instance->dist[i * nn + succ[i]];
					// printf("relative mileage for node %i in succ wrt node %i in unc_nodes is %f\n", i, uncovered_nodes[h], relative_mileage);
					if (relative_mileage < absolute_mileage)
					{
						absolute_mileage = relative_mileage;
						min_node_index = h;
						min_node_succ = i;
					}
				}
			}
		}
	}
	emstruct ems;
	ems.succ_node = min_node_succ;
	ems.index_uncovered_node = min_node_index;
	*current_cost = *current_cost + absolute_mileage;
	return ems;
}

void extra_mileage(Instance* instance, int* successors, int* size, double* cost)
{
	if (instance->geneticFlag == 0)
	{
		// in a generic iteration, compute extra mileage for all covered nodes and all uncovered nodes and pick the minimum e-m:
		// E-M(i,h) = c[i,h] + c[h, succ[i]] - c[i, succ[i]]

		int* uncovered_nodes = (int*)calloc(instance->nnodes, sizeof(int));

		int max_node_index = 0;

		int succ_counter = 0;

		for (int j = 0; j < instance->nnodes; j++)
		{
			for (int i = 0; i < instance->nnodes; i++)
			{
				// filling the uncovered_nodes array
				uncovered_nodes[i] = i;
				// flag value for succ's elements
				instance->emSucc[i] = -1;
			}
			int starting_node = j;
			if (instance->verbose > 50)
			{
				printf("EM first node: %d\n", starting_node);
			}
			// cost value of the j-th cycle, initially set to 0
			double current_cost = 0.0;

			uncovered_nodes[starting_node] = -1;

			// max_node_index will hold the index of the node in uncovered_nodes at maximum distance from starting_node
			max_node_index = max_dist_node(starting_node, uncovered_nodes, instance, &current_cost);

			instance->emSucc[starting_node] = uncovered_nodes[max_node_index];
			instance->emSucc[max_node_index] = starting_node;
			succ_counter = 2;
			uncovered_nodes[max_node_index] = -1;

			while (succ_counter < instance->nnodes)
			{
				emstruct ems = min_mileage(instance, instance->emSucc, uncovered_nodes, &current_cost);
				if (instance->verbose > 50)
				{
					printf("EM next node: %d\n", ems.succ_node);
				}
				int tmp = instance->emSucc[ems.succ_node];
				instance->emSucc[ems.succ_node] = uncovered_nodes[ems.index_uncovered_node];
				instance->emSucc[uncovered_nodes[ems.index_uncovered_node]] = tmp;

				succ_counter += 1;
				uncovered_nodes[ems.index_uncovered_node] = -1;
			}
			for (int i = 0; i < instance->nnodes; i++)
			{
				if (instance->emSucc[i] == -1)
				{
					instance->emSucc[i] = j;
					current_cost += instance->dist[i * instance->nnodes + j];
					break;
				}
			}
			if (current_cost < instance->emCost)
			{
				instance->emCost = current_cost;
				if (instance->verbose >= 10) 
				{ 
					printf("New best cost: %f\n", instance->emCost); 
				}
			}
			clock_t end = clock();
			double diff = (double)((end - instance->tstart) / CLOCKS_PER_SEC);
			if (diff > instance->timeLimit)
			{
				break;
			}
		}
		free(uncovered_nodes);
		if (instance->verbose >= 10)
		{
			printSolution(instance, instance->emSucc);
		}
	}

	else
	{
		int* uncovered_nodes = (int*)calloc(instance->nnodes, sizeof(int));

		int max_node_index = 0;
		int min_mileage_index = 0;
		int succ_counter = 0;

		for (int i = 0; i < instance->nnodes; i++)
		{
			uncovered_nodes[i] = i;
		}

		for (int i = 0; i < instance->nnodes; i++)
		{
			// if node i already in chromosome then uncovered_nodes[i] = -1
			int tmp = successors[i];
			if (tmp != -1)
			{
				uncovered_nodes[tmp] = -1;
			}
		}

		while (*size < instance->nnodes)
		{
			emstruct ems = min_mileage(instance, successors, uncovered_nodes, cost);

			int tmp = successors[ems.succ_node];
			successors[ems.succ_node] = uncovered_nodes[ems.index_uncovered_node];
			successors[uncovered_nodes[ems.index_uncovered_node]] = tmp;

			*size += 1;
			uncovered_nodes[ems.index_uncovered_node] = -1;
		}
		free(uncovered_nodes);
	}
}