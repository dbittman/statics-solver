#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <string.h>
#include <math.h>
#include "grid.h"

double jacobi_solve(struct grid *grid);
double solve(struct grid *grid, int method)
{
	printf("Started SOLVER - gridsz %d, method %d\n", grid->len, method);
	double res = jacobi_solve(grid);
	printf("SOLVER completed: %lf\n", res);
	return res;
}

