#include <cplex.h>  
#include <float.h>
#include <time.h>

#include "utilities.h"
#include "tsp.h"
#include "constructing_heuristics.h"
#include "meta-heuristics.h"

#define VERBOSE 9
#define LOCAL_BRANCHING_MIN_IMPROVEMENT 0.015 //1.5%
#define LOCAL_BRANCHING_MAX_SMALL_IMPROVEMENT 3


void print_error(const char* err)
{
	printf("\n\n ERROR: %s \n\n", err);
	fflush(NULL);
	exit(1);
}

void updateUB(Instance* inst, int succ[], double* UB)
{
	double cost = 0.0;
	int* sol = (int*)calloc(inst->nnodes, sizeof(int));
	succ2sol(inst, inst->cplexSucc, sol);
	for (int i = 0; i < inst->nnodes-1; i++)
	{
		cost += inst->dist[sol[i]*inst->nnodes + sol[i+1]];
	}
	cost += inst->dist[sol[inst->nnodes-1] * inst->nnodes + sol[0]];
	if (cost <= *UB)
	{
		*UB = cost;
	}
}

void computeCost(Instance* inst, int succ[], double* UB)
{
	double cost = 0.0;
	
	for (int i = 0; i < inst->nnodes; i++)
	{
		cost += inst->dist[i* inst->nnodes + succ[i]];
	}
	*UB = cost;
}

float randomFloat()
{
	float ran = (float)rand() / (float)RAND_MAX;
	return ran;
}


int xpos(int i, int j, Instance* inst)                                              
{
	if (i == j) print_error(" i == j in xpos");
	if (i > j) return xpos(j, i, inst);
	int pos = i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
	return pos;
}


/***************************************************************************************************************************/
void build_model(Instance* inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{

	double zero = 0.0;
	char binary = 'B';

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// add binary var.s x(i,j) for i < j  

	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);  		// ... x(1,2), x(1,3) ....
			double obj = inst->dist[i*inst->nnodes +j]; // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(" wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst)) print_error(" wrong position for x var.s");
		}
	}

	// add the degree constraints 

	int* index = (int*)calloc(inst->nnodes, sizeof(int));
	double* value = (double*)calloc(inst->nnodes, sizeof(double));

	for (int h = 0; h < inst->nnodes; h++)  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h + 1);
		int nnz = 0;
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			index[nnz] = xpos(i, h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0])) print_error("CPXaddrows(): error 1");
	}

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	if (VERBOSE >= 100) 
	{
		int status = CPXwriteprob(env, lp, "model.lp", NULL);
		if (status)
		{
			printf("Status error: %d\n", status);
			exit(1);
		}
	}
}

//#define DEBUG    // comment out to avoid debugging 
//#define EPS 1e-5

// Bender's loop method:
// receives integer solution made of one or mode cycles
//

/*********************************************************************************************************************************/
void build_sol(const double* xstar, Instance* inst, int* succ, int* comp, int* ncomp) // build succ() and comp() wrt xstar()...
/*********************************************************************************************************************************/
{

#ifdef DEBUG // I check if degree of all nodes is two and if xstar does not contain fract. components
	int* degree = (int*)calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			int k = xpos(i, j, inst);
			if (fabs(xstar[k]) > EPS && fabs(xstar[k] - 1.0) > EPS) print_error(" wrong xstar in build_sol()");
			if (xstar[k] > 0.5)
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for (int i = 0; i < inst->nnodes; i++)
	{
		if (degree[i] != 2) print_error("wrong degree in build_sol()");
	}
	free(degree);
#endif

	* ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	for (int start = 0; start < inst->nnodes; start++)
	{
		if (comp[start] >= 0) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while (!done)  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for (int j = 0; j < inst->nnodes; j++)
			{
				if (i != j && xstar[xpos(i, j, inst)] > 0.5 && comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}
		succ[i] = start;  // last arc to close the cycle

		// go to the next component...
	}
}

void patching(Instance* inst, int succ[], int comp[], int ncomp)
{
	while (ncomp > 1)
	{
		double deltaCost = DBL_MAX;
		double costNew = 0.0;
		double costOld = 0.0;
		int tmpA = 0, tmpB = 0;
		int a = 0, sa = 0, b = 0, sb = 0;
		for (int i = 1; i <= ncomp; i++)
		{
			for (int j = 0; j < inst->nnodes; j++)
			{
				if (comp[j] == i)
				{
					for (int k = j + 1; k < inst->nnodes; k++)
					{
						if (comp[j] < comp[k])
						{
							costOld = inst->dist[j * inst->nnodes + succ[j]] + inst->dist[k * inst->nnodes + succ[k]];
							costNew = inst->dist[j * inst->nnodes + succ[k]] + inst->dist[k * inst->nnodes + succ[j]];

							if (costNew - costOld < deltaCost)
							{
								deltaCost = costNew - costOld;
								a = j;
								sa = succ[j];
								b = k;
								sb = succ[k];
							}
						}
					}
				}
			}
		}
		// update components
		int compA = comp[a];
		int tmp = sb;
		while (tmp != b)
		{
			comp[tmp] = compA;
			tmp = succ[tmp];
		}
		comp[tmp] = compA;
		// update successors
		succ[a] = sb;
		succ[b] = sa;
		// check if there are more than two components
		int prev = comp[0];
		int flag = 0;
		for (int i = 1; i < inst->nnodes; i++)
		{
			tmp = comp[i];
			if (tmp != prev)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			break;
		}
		
	}
}


// alternative to Benders'
/********************************************************************************************************/
static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
/********************************************************************************************************/
{
	// +++++++++
	// Threadsafe: no variable inside inst must be modified
	// rand() is not thread safe!!!!! Create a thread safe RNG yourself
	// +++++++++
	//printf("\n+++++++++++++++++++++++++++++++\n");
	clock_t start = clock();
	Instance* inst = (Instance*)userhandle;
	double* xstar = (double*)malloc(inst->ncols * sizeof(double));  // current candidate solution will go here: IT'S A INTEGER SOLUTION (collection of cycles)
	double objval = CPX_INFBOUND;
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval)) print_error("CPXcallbackgetcandidatepoint error");

	// get some random information at the node (as an example for the students)
	int mythread = -1;
	CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread); //thread starts from 0
	// if you want to know at which node of the branching tree this function has been called:
	int mynode = -1;
	CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
	// value of the incumbent (cost of best solution) when this function (my_callback) is called:
	double incumbent = CPX_INFBOUND;
	CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
	//if ( VERBOSE >= 100 ) printf(" ... callback at node %5d thread %2d incumbent %10.2lf, candidate value %10.2lf\n", .....);


	int nnz = 0;
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* index = (int*)calloc(inst->ncols, sizeof(int)); //  put in instance the number of columns
	double* coeff = (double*)calloc(inst->ncols, sizeof(double));
	int ncomp = 0;
	char sense = 'L';
	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(inst->nnodes, sizeof(char));
	int izero = 0;
	build_sol(xstar, inst, succ, comp, &ncomp);
	//if xstar is infeasible, find a violated cut and store it in the usual Cplex's data structure (rhs, sense, nnz, index and value)
	//as we did in build_sol
	if (ncomp <= 1)
	{
		free(xstar);
		free(succ);
		free(comp);
		free(index);
		free(coeff);
		free(cname[0]);
		free(cname);
		return 0;
	}
	for (int k = 1; k <= ncomp; k++)
	{
		//check if overall time is over
		clock_t end = clock();
		double diff = (double)(end - inst->tstart) / CLOCKS_PER_SEC;
		//If time is over, patch the solution
		if (diff > inst->timeLimit)
		{			
			if (inst->patchFlag == 0)
			{
				inst->patchFlag = 1;
				patching(inst, succ, comp, ncomp);
				int flag = 1;
				double cost = DBL_MAX;
				updateUB(inst, succ, &cost);
				while (flag == 1)
				{
					flag = two_opt(inst, succ, &cost, 0, 0);
				}
				for (int index = 0; index < inst->nnodes; index++)
				{
					inst->cplexSucc[index] = succ[index];
				}
				if (inst->verbose >= 10) 
				{
					printf("New best cost: %f\n", cost);
				}
				inst->cplexCost = cost;
			}
			break;
		}

		int nnz = 0;
		double rhs = -1;
		for (int l = 0; l < inst->ncols; l++)
		{
			index[l] = 0;
			coeff[l] = 0.0;
		}
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (comp[i] == k)
			{
				rhs += 1;
				for (int j = i + 1; j < inst->nnodes - 1; j++)
				{
					if (comp[j] == k)
					{
						index[nnz] = xpos(i, j, inst);
						coeff[nnz] = 1.0;
						nnz += 1;
					}
				}
			}
		}
		if (nnz >= 2) // means that the solution is infeasible and a violated cut has been found -> insert SEC's
		{
			/*if (inst->verbose >= 10) 
			{
				printf("the solution is infeasible and a violated cut has been found\n");
			}*/
			int izero = 0;
			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, coeff)) print_error("CPXcallbackrejectcandidate() error"); // reject the solution and adds one cut 
		}
	}
	
	free(xstar);
	free(succ);
	free(comp);
	free(index);
	free(coeff);
	free(cname[0]);
	free(cname);
	return 0;
}

int TSPCallback(Instance* inst)
{
	inst->tstart = clock();
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);

	double LB = -CPX_INFBOUND;
	double UB = 0.0;
	inst->greedySucc = (int*)calloc(inst->nnodes, sizeof(int));
	greed_search(inst, 0, 0, NULL, 0, &UB);
	
	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);
	
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timeLimit); //MODIFY THIS FOR THE TIME LIMIT
	CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);

	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;
	// install callbacks (BEFORE CALLING CPXmipopt):
	// install a "lazyconstraint" callback to cut infeasible integer sol.s (found e.g. by heuristics)
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) print_error("CPXcallbacksetfunc() error");
	// I need to call the callback only when CPLEX finds a candidate integer solution
	error = CPXmipopt(env, lp);
	if (error)
	{
		//printf("CPX error code %d\n", error);
		printf("Error is %d!\n", error);
		print_error("CPXmipopt() error");
	}

	// Check if time is over
	clock_t end = clock();
	double diff = (double)(end - inst->tstart) / CLOCKS_PER_SEC;
	if (diff > inst->timeLimit)
	{
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);
		return 0;
	}

	double* xstar = (double*)calloc(ncols, sizeof(double));
	error = CPXgetx(env, lp, xstar, 0, ncols - 1);
	if (error)
	{
		printf("Error is %d!!\n", error);
		print_error("CPXgetx() error");
	}
	
	CPXgetbestobjval(env, lp, &UB);
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int ncomp = 0;
	build_sol(xstar, inst, succ, comp, &ncomp);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	//update UB and inst->cplexSucc
	for (int i = 0; i < inst->nnodes; i++)
	{
		inst->cplexSucc[i] = succ[i];
	}
	updateUB(inst, succ, &UB);
	inst->cplexCost = UB;

	free(xstar);
	free(succ);
	free(comp);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int CPXPUBLIC ad_hoc(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
{
	clock_t start = clock();
	Instance* inst = (Instance*)userhandle;
	double* xstar = (double*)calloc(inst->ncols, sizeof(double));  // current candidate solution will go here: IT'S A INTEGER SOLUTION (collection of cycles)
	double objval = CPX_INFBOUND;
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval)) print_error("CPXcallbackgetcandidatepoint error");

	int nnz = 0;
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* index = (int*)calloc(inst->ncols, sizeof(int)); //  put in instance the number of columns
	double* coeff = (double*)calloc(inst->ncols, sizeof(double));
	int ncomp = 0;
	char sense = 'L';
	int izero = 0;
	build_sol(xstar, inst, succ, comp, &ncomp);
	//if xstar is infeasible, find a violated cut and store it in the usual Cplex's data structure (rhs, sense, nnz, index and value)
	//as we did in build_sol
	if (ncomp <= 1)
	{
		
		return 0;
	}
	
	double* xheu = (double*)calloc(inst->ncols, sizeof(double));
	int* ind = (int*)calloc(inst->ncols, sizeof(int));
	for (int i = 0; i < inst->ncols; i++)
	{
		ind[i] = i;
		xheu[i] = 0.0;
	}

	patching(inst, succ, comp, ncomp);
	int flag = 1;
	while (flag == 1)
	{
		flag = two_opt(inst, succ, &objval, 0, 0);
	}
	//Check if time is over
	clock_t end = clock();
	double diff = (double)((end - inst->tstart)) / CLOCKS_PER_SEC;
	if (diff > inst->timeLimit)
	{
		if (inst->patchFlag == 0)
		{
			inst->patchFlag = 1;
			patching(inst, succ, comp, ncomp); 
			flag = 1;
			double cost = DBL_MAX;
			updateUB(inst, succ, &cost);
			while (flag == 1)
			{
				flag = two_opt(inst, succ, &cost, 0, 0);
			}
			for (int cntr = 0; cntr < inst->nnodes; cntr++)
			{
				inst->best_succ[cntr] = succ[cntr];
			}
			inst->cplexCost = cost;
			if (inst->verbose >= 10)
			{
				printf("New best cost: %f\n", cost);
			}
		}
		free(ind);
		free(xstar);
		free(succ);
		free(comp);
		free(index);
		free(coeff);
		free(xheu);
		return 0;
	}
	for (int i = 0; i < inst->nnodes; i++)
	{
		xheu[xpos(i, succ[i], inst)] = 1.0;
	}
	double UB = DBL_MAX;
	updateUB(inst, succ, &UB);

	if (CPXcallbackpostheursoln(context, inst->ncols, ind, xheu, UB, CPXCALLBACKSOLUTION_NOCHECK)) print_error("CPXcallbackpostheursoln() error");

	free(ind);
	free(xstar);
	free(succ);
	free(comp);
	free(index);
	free(coeff);
	free(xheu);
	return 0;
}

int TSPadHoc(Instance* inst)
{
	inst->tstart = clock();
	inst->starting_node = 0;
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);

	double LB = -CPX_INFBOUND;
	double UB = 0.0;

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);

	inst->greedySucc = (int*)calloc(inst->nnodes, sizeof(int));
	greed_search(inst, 0, 0, NULL, 0, &UB);

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timeLimit); //MODIFY THIS FOR THE TIME LIMIT
	CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);

	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;
	// install callbacks (BEFORE CALLING CPXmipopt):
	// install a "lazyconstraint" callback to cut infeasible integer sol.s (found e.g. by heuristics)
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (CPXcallbacksetfunc(env, lp, contextid, ad_hoc, inst)) print_error("CPXcallbacksetfunc() error");

	error = CPXmipopt(env, lp);
	if (error)
	{
		//printf("CPX error code %d\n", error);
		printf("Error is %d!\n", error);
		print_error("CPXmipopt() error");
	}

	 //Check if time is over
	clock_t end = clock();
	double diff = (double)(end - inst->tstart) / CLOCKS_PER_SEC;
	if (diff > inst->timeLimit)
	{
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);
		return 0;
	}

	double* xstar = (double*)calloc(ncols, sizeof(double));
	error = CPXgetx(env, lp, xstar, 0, ncols - 1);
	if (error)
	{
		printf("Error is %d!!\n", error);
		print_error("CPXgetx() error");
	}

	CPXgetbestobjval(env, lp, &UB);
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int ncomp = 0;
	build_sol(xstar, inst, succ, comp, &ncomp);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	//update UB and inst->cplexSucc
	for (int i = 0; i < inst->nnodes; i++)
	{
		inst->cplexSucc[i] = succ[i];
	}
	inst->cplexCost = UB;
	
	free(xstar);
	free(succ);
	free(comp);
}

int TSPHardFixing(Instance* inst)
{
	inst->tstart = clock();
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timeLimit);

	//use greedy to obtain a starting solution, which will be saved in inst->best_sol, and from there apply hard fixing
	// STOPPING CONDITION: On UB and time
	double UB = 0.0;
	double prev_UB = CPX_INFBOUND;
	double LB = -CPX_INFBOUND;
	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;
	inst->greedySucc = (int*)calloc(inst->nnodes, sizeof(int)); // will be freed by free_memory() in main;
	greed_search(inst, 0, 0, NULL, 0, &UB);
	
	CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);
	CPXsetintparam(env, lp, CPX_PARAM_NODELIM, 1000);
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) print_error("CPXcallbacksetfunc() error");

	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* sol = (int*)calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[i] = inst->greedySucc[i];
	}
	succ2sol(inst, succ, sol);

	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int ncomp = 0;
	char lu = 'L';
	double bd_one = 1.0;
	double bd_zero = 0.0;
	while (1)
	{
		// fix lower bound
		//use xpos to get the column position
		for (int i = 0; i < inst->nnodes; i++)
		{
			float probability = (float)rand()/RAND_MAX;

			if (probability < 0.2)
			{
				if (i == inst->nnodes - 1)
				{
					int position = xpos(sol[i], sol[0], inst);
					CPXchgbds(env, lp, 1, &position, &lu, &bd_one);
				}
				else
				{
					int position = xpos(sol[i], sol[i+1], inst);
					CPXchgbds(env, lp, 1, &position, &lu, &bd_one);
				}
			}
		}
		// solve
		error = CPXmipopt(env, lp);
		if (error)
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error");
		}
		// Check if time is over
		clock_t end = clock();
		double diff = (double)(end - inst->tstart) / CLOCKS_PER_SEC;
		
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timeLimit - diff);
		if (diff > inst->timeLimit)
		{
			free(succ);
			free(comp);
			free(sol);
			break;
		}
		// if time is not over, update UB
		double* xstar = (double*)calloc(ncols, sizeof(double));
		error = CPXgetx(env, lp, xstar, 0, ncols - 1);
		if (error)
		{
			printf("Error is %d!!\n", error);
			print_error("CPXgetx() error");
			
		}

		build_sol(xstar, inst, succ, comp, &ncomp);
			
		CPXgetbestobjval(env, lp, &UB);
		double handCompCost = 0.0;
		computeCost(inst, succ, &handCompCost);

		if (UB < prev_UB)
		{
			prev_UB = UB;
			
			for (int i = 0; i < inst->nnodes; i++)
			{
				inst->cplexSucc[0] = succ[0];
			}
			succ2sol(inst, succ, sol);
			inst->cplexCost = UB;
			if (inst->verbose >= 10)
			{
				printf("New best cost: %f\n", UB);
			}
		}
		free(xstar);

		for (int i = 0; i < inst->nnodes; i++)
		{
			for (int j = i + 1; j < inst->nnodes; j++)
			{
				int position = xpos(i, j, inst);
				CPXchgbds(env, lp, 1, &position, &lu, &bd_zero);
			}
		}
	}
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
}


void TSPLocalBranching(Instance* inst) 
{
	inst->tstart = clock();
	if (VERBOSE >= 10) printf("\nStart Local branching matheuristic\n");

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");
	//CPXsetdblparam(env, CPXPARAM_TimeLimit, inst->timeLimit);

	build_model(inst, env, lp);

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);
	//CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timeLimit);

	//first we need to find a feasible solution
	long ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;
	int* indexes = (int*)calloc(inst->ncols, sizeof(int));
	double* values = (double*)calloc(inst->ncols, sizeof(double));
	double* xh = (double*)calloc(inst->ncols, sizeof(double));

	double UB = 0.0;
	inst->greedySucc = (int*)calloc(inst->nnodes, sizeof(int));
	greed_search(inst, 0, 0, NULL, 0, &UB);
	CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);
	
	for (int i = 0; i < inst->nnodes; i++) 
	{
		xh[xpos(i,inst->greedySucc[i], inst)] = 1.0;  
	}

	//Setting up indices array. Setting up it here avoids on setting it up everytime the 2-opt callback need it. One time initialization and that's all.
	int* ind = (int*)calloc(inst->ncols, sizeof(int));
	int k = 0;
	for (int i = 0; i < inst->ncols; i++)
	{
		ind[i] = i;
	}

	int beg = 0;
	int level = CPX_MIPSTART_NOCHECK;
	if (CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, ind, xh, &level, NULL)) print_error("CPXaddmipstarts() error");


	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) print_error("CPXcallbacksetfunc() error");

	double K[] = { 10, 20, 30, 40 };
	int k_index = 0;
	double actual_val = 0.0;
	double actual_best = UB;

	int small_imp = 0;
	char sense = 'G'; // >=
	int matbeg = 0;
	int num_iter = 0; //counter of iteration
	double bd_one = 1.0;
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* sol = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[i] = inst->greedySucc[i];
	}
	succ2sol(inst, inst->greedySucc, sol);
	char** cname = (char**)calloc(1, sizeof(char));
	cname[0] = (char*)calloc(100, sizeof(char));
	sprintf(cname[0], "Local branching constraint");

	while (1) 
	{
		if (k_index >= (sizeof(K) / sizeof(*K))) { break; }

		clock_t end = clock();
		double diff = (double)(end - inst->tstart) / CLOCKS_PER_SEC;
		if (diff > inst->timeLimit) { break; } //stop

		double time_remain = inst->timeLimit - diff;

		CPXsetdblparam(env, CPXPARAM_TimeLimit, time_remain);

		int k = 0;
		for (int i = 0; i < inst->nnodes; i++) 
		{
			if (xh[xpos(i,succ[i], inst)] > 0.5)
			{
				indexes[k] = xpos(i,succ[i], inst);
				values[k] = 1.0;
				k++;
			}
		}

		double rhs = inst->nnodes - K[k_index];
		if (CPXaddrows(env, lp, 0, 1, k, &rhs, &sense, &matbeg, indexes, values, NULL, &cname[0])) print_error("CPXaddrows() error");
		//solve the model
		error = CPXmipopt(env, lp);
		if (error)
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error");
		}
		CPXgetx(env, lp, xh, 0, ncols - 1);

		build_sol(xh, inst, succ, comp, &ncomp);

		CPXgetobjval(env, lp, &actual_val);

		double improved_cost = 1 - actual_val / actual_best;

		if (actual_val < actual_best) 
		{
			if (improved_cost < LOCAL_BRANCHING_MIN_IMPROVEMENT) 
			{
				small_imp++;
			}
			else //se ho capito cosa volevi fare, penso che l'else vada qui, ma se mi sono sbagliato mettilo pure dov'era prima
			{
				small_imp = 0;
			}
			if (small_imp % LOCAL_BRANCHING_MAX_SMALL_IMPROVEMENT == 0 && k_index < (sizeof(K) / sizeof(*K)) - 1) 
			{
				k_index++;
			}
			
			//Update solution
			actual_best = actual_val;
			for (int i = 0; i < inst->nnodes; i++)
			{
				inst->cplexSucc[i] = succ[i];
			}
			succ2sol(inst, succ, sol);
			inst->cplexCost = actual_val;
			if (inst->verbose >= 10)
			{
				printf("New best cost: %f\n", inst->cplexCost);
			}
		}

		// Remove the added soft-fixing constraints
		int numrows = CPXgetnumrows(env, lp);
		error = CPXdelrows(env, lp, numrows - 1, numrows - 1);
		if (error)
		{
			printf("CPXdelrows error code %d\n", error);
			print_error("CPXdelrows() error");
		}
		num_iter++;
		for (int i = 0; i < ncols; i++) // "clean" indexes and values
		{
			indexes[i] = 0;
			values[i] = 0.0;
		}
	}
	
	free(cname[0]);
	free(cname);
	free(comp);
	free(succ);
	free(sol);
	free(xh);
	free(values);
	free(indexes);
	free(ind);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
}

/**************************************************************************************************************************/
int TSPopt(Instance* inst)
/**************************************************************************************************************************/
{
	inst->tstart = clock();
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);
	CPXsetdblparam(env, CPX_PARAM_TILIM, 3600.0);
	

	//CPXsetintparam(env, CPX_PARAM_THREADS, 1); 	// single thread -> just for debugging
	double UB = 0.0;
	inst->greedySucc = (int*)calloc(inst->nnodes, sizeof(int)); // will be freed by free_memory() in main;
	greed_search(inst, 0, 0, NULL, 0, &UB);
	
	CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);
	CPXsetintparam(env, lp, CPX_PARAM_NODELIM, 100);
	error = CPXmipopt(env, lp);
	if (error)
	{
		printf("CPX error code %d\n", error);
		print_error("CPXmipopt() error");
	}

	// use the optimal solution found by CPLEX
	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;
	double* xstar = (double*)calloc(ncols, sizeof(double));
	if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error");
	// Print xstar
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			//if (xstar[xpos(i, j, inst)] > 0.5) printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
		}
	}

	//Benders' loop method
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* index = (int*)calloc(ncols, sizeof(int));
	double* coeff = (double*)calloc(ncols, sizeof(double));
	int ncomp = 0;
	char sense = 'L';
	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(inst->nnodes, sizeof(char));
	int izero = 0;
	int Sname = 0;
	double LB = -CPX_INFBOUND;
	
	int iteration = 0;
	//TODO: UPDATE UB
	while (LB < 0.999*UB)
	{
		//printf("LB at iteration %d is: %f\n", iteration, LB);
		iteration++;
		//update CPLEX time limit
		clock_t end = clock();
		double diff = (double)(end - inst->tstart) / CLOCKS_PER_SEC;
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timeLimit - diff); 
		// build solution -> build succ, comp and get ncomp
		build_sol(xstar, inst, succ, comp, &ncomp); 
		computeCost(inst, succ, &UB);
		//printf("Iteration %d; UB: %f, LB: %f\n", iteration, UB, LB);
		if (ncomp <= 1)
		{
			error = CPXgetobjval(env, lp, &UB);
			if (error)
			{
				printf("Error: %d\n", error);
				print_error("CPXgetobj val error");
			}
			break;
		}
		else
		{
			for (int k = 1; k <= ncomp; k++)
			{
				Sname += 1;
				int nnz = 0;
				double rhs = -1;
				for (int l = 0; l < ncols; l++)
				{
					index[l] = 0;
					coeff[l] = 0.0;
				}
				for (int i = 0; i < inst->nnodes; i++)
				{
					if (comp[i] == k)
					{
						rhs += 1;
						for (int j = i + 1; j < inst->nnodes - 1; j++)
						{
							if (comp[j] == k)
							{
								index[nnz] = xpos(i, j, inst);
								coeff[nnz] = 1.0;
								nnz += 1;
								sprintf(cname[0], "S%d(%d,%d)",Sname, i + 1, j + 1);
							}
						}
					}
				}
				if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, coeff, NULL, &cname[0])) print_error("CPXaddrows(): error 1");
			}

			//check if time limit is over
			clock_t end = clock();
			double diff = (double)(end - inst->tstart) / CLOCKS_PER_SEC;
			//printf("Diff - timelimit:%f - %f \n", diff, inst->timeLimit);
			if (diff - inst->timeLimit > 0.0) // time is over
			{
				patching(inst, succ, comp, ncomp);
				int flag = 1;
				while (flag == 1)
				{
					flag = two_opt(inst, succ, &UB, 0, 0);
				}
				
				for (int index = 0; index < inst->nnodes; index++)
				{
					inst->cplexSucc[index] = succ[index];
				}
				inst->cplexCost = UB;
				if (inst->verbose >= 10)
				{
					printf("New best cost: %f\n", inst->cplexCost);
				}
				break;
			}
		}
		// Obtain the new cplex solution
		error = CPXmipopt(env, lp);
		LB = max(CPXgetobjval(env, lp, &LB), LB);
		
		// Obtain again xstar (like an update)
		error = CPXgetx(env, lp, xstar, 0, ncols - 1);
		if (error)
		{
			printf("CPX error code %d\n", error);
			print_error("CPXgetx() error");
		}
	}

	free(xstar);
	free(succ);
	free(comp);
	free(index);
	free(coeff);
	free(cname[0]);
	free(cname);

	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0; // or an appropriate nonzero error code
}
/***************************************************************************************************************************/


